/* GPL-3.0-or-later. Conservative sectional aerosol replay on existing OpenFOAM
 * mesh and time-interpolated native mass face fluxes. Air fields are read-only.
 * Wall closure is an explicitly approximate small-Stokes drift/resistance model,
 * not an implementation claimed identical to any published named model. */
#include "fvCFD.H"
#include "wallFvPatch.H"
#include "inletOutletFvPatchFields.H"
#include <cmath>

scalar conductance(const scalar v, const scalar resistance)
{
    const scalar x=v*resistance;
    if (mag(x)<1e-6) return (1+x/2+sqr(x)/12)/resistance;
    if (x < -700) return 0;
    return v/(-std::expm1(-x));
}

int main(int argc,char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    if(Pstream::parRun())FatalErrorInFunction<<"Serial replay only"<<exit(FatalError);
    IOdictionary cfg(IOobject("ceilingProperties",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE));
    const scalarList frameTimes(cfg.lookup("frameTimes"));
    const wordList frameDirs(cfg.lookup("frameDirectories"));
    const wordList names(cfg.lookup("fieldNames"));
    const scalarList ds(cfg.lookup("diameters"));
    const scalar sticking=cfg.getOrDefault<scalar>("ceilingSticking",1);
    const scalar otherSticking=cfg.getOrDefault<scalar>("otherSticking",1);
    const scalar ScT=cfg.getOrDefault<scalar>("SchmidtT",.7);
    const scalar inertia=cfg.getOrDefault<scalar>("inertiaFactor",1);
    const scalar turbophoresis=cfg.getOrDefault<scalar>("turbophoresisFactor",0);
    const scalar impaction=cfg.getOrDefault<scalar>("impactionFactor",1);
    const word impactionMode=cfg.getOrDefault<word>("impactionMode","resistance");
    const scalar thermophoresis=cfg.getOrDefault<scalar>("thermophoresisFactor",1);
    const scalar wallDT=cfg.getOrDefault<scalar>("wallTemperatureOffsetK",0);
    const scalar updateEvery=cfg.getOrDefault<scalar>("coefficientUpdateSeconds",1);
    const scalar emit=cfg.getOrDefault<scalar>("emissionKgPerSecPerBin",1e-9);
    const scalar stop=cfg.getOrDefault<scalar>("sourceStopTime",120);
    const label nquad=cfg.getOrDefault<label>("wallQuadrature",72);
    const bool legacy=cfg.getOrDefault<bool>("legacyWallModel",false);
    const bool initialSaved=cfg.getOrDefault<bool>("useSavedInitialAerosols",false);
    const scalar rhoOil=850,kB=1.380649e-23,RR=6.0221417930e23*1.38065e-23*1000/28.965;
    if(frameTimes.size()<2 || frameDirs.size()!=frameTimes.size())FatalErrorInFunction<<"Need >=2 matching frames"<<exit(FatalError);
    label bracket=-1;
    autoPtr<volVectorField> ulo,uhi;
    autoPtr<volScalarField> tlo,thi,plo,phiP,klo,khi,nlo,nhi;
    autoPtr<surfaceScalarField> flo,fhi;
    auto load=[&](label i){
        auto readV=[&](const word& n,label j){return new volVectorField(IOobject(n,"frames/"+frameDirs[j],mesh,IOobject::MUST_READ,IOobject::NO_WRITE,false),mesh);};
        auto readS=[&](const word& n,label j){return new volScalarField(IOobject(n,"frames/"+frameDirs[j],mesh,IOobject::MUST_READ,IOobject::NO_WRITE,false),mesh);};
        ulo.reset(readV("U",i));uhi.reset(readV("U",i+1));
        tlo.reset(readS("T",i));thi.reset(readS("T",i+1));plo.reset(readS("p",i));phiP.reset(readS("p",i+1));
        klo.reset(readS("k",i));khi.reset(readS("k",i+1));nlo.reset(readS("nut",i));nhi.reset(readS("nut",i+1));
        flo.reset(new surfaceScalarField(IOobject("phi","frames/"+frameDirs[i],mesh,IOobject::MUST_READ,IOobject::NO_WRITE,false),mesh));
        fhi.reset(new surfaceScalarField(IOobject("phi","frames/"+frameDirs[i+1],mesh,IOobject::MUST_READ,IOobject::NO_WRITE,false),mesh));bracket=i;
    };
    load(0);
    volVectorField U(IOobject("Ureplay",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),ulo());
    volScalarField T(IOobject("Treplay",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),tlo());
    volScalarField p(IOobject("preplay",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),plo());
    volScalarField k(IOobject("kreplay",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),klo());
    volScalarField nut(IOobject("nutreplay",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),nlo());
    volScalarField rho(IOobject("rhoReplay",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),p/(dimensionedScalar("R",sqr(dimVelocity)/dimTemperature,RR)*T));
    surfaceScalarField phi(IOobject("phiReplay",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),flo());
    volScalarField shape(IOobject("sourceShape",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
    PtrList<volScalarField> Y(names.size()),D(names.size()),sink(names.size());
    PtrList<volVectorField> drift(names.size());
    PtrList<surfaceScalarField> J(names.size()),M(names.size());
    List<List<scalarField>> vdep(names.size());
    List<scalarField> ceiling(mesh.boundary().size());
    scalar ceilingArea=0;
    forAll(mesh.boundary(),pi){const auto& pp=mesh.boundary()[pi];ceiling[pi].setSize(pp.size(),0);forAll(pp,fi){if(isA<wallFvPatch>(pp)&&mag(pp.Cf()[fi].z()-2.2)<1e-7&&pp.nf()()[fi].z()>.99){ceiling[pi][fi]=1;ceilingArea+=pp.magSf()[fi];}}}
    scalarList released(names.size(),0),hood(names.size(),0),left(names.size(),0),right(names.size(),0),deposited(names.size(),0),ceilingMass(names.size(),0),initial(names.size(),0);
    forAll(names,bi){
        Y.set(bi,new volScalarField(IOobject(names[bi],runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh));
        if(!initialSaved)Y[bi]=dimensionedScalar("zero",dimless,0);
        initial[bi]=gSum(rho.primitiveField()*Y[bi].primitiveField()*mesh.V());mesh.setFluxRequired(names[bi]);
        D.set(bi,new volScalarField(IOobject("diffusion"+names[bi],runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedScalar("d",dimDensity*dimArea/dimTime,0)));
        sink.set(bi,new volScalarField(IOobject("sink"+names[bi],runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedScalar("s",dimDensity/dimTime,0)));
        drift.set(bi,new volVectorField(IOobject("drift"+names[bi],runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedVector("d",dimVelocity,Zero)));
        J.set(bi,new surfaceScalarField(IOobject("Jdep_"+names[bi],runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),mesh,dimensionedScalar("j",dimMass/dimArea/dimTime,0)));
        M.set(bi,new surfaceScalarField(IOobject("Mdep_"+names[bi],runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),mesh,dimensionedScalar("m",dimMass/dimArea,0)));
        vdep[bi].setSize(mesh.boundary().size());forAll(mesh.boundary(),pi)vdep[bi][pi].setSize(mesh.boundary()[pi].size(),0);
    }
    mkDir(runTime.path()/"diagnostics");
    OFstream budget(runTime.path()/"diagnostics/ceiling_budget.csv");budget.precision(14);
    budget<<"time_s,diameter_um,initial_kg,released_kg,stock_kg,hood_out_kg,left_out_kg,right_out_kg,ceiling_kg,other_wall_kg,mass_error_kg,ceiling_rate_kg_s,Qhood_m3_s,air_continuity_L1_kg_s,max_cell_Stokes,ceiling_yplus_max,ceiling_delta_max_m,min_Y"<<endl;
    scalar nextCoeff=-GREAT,maxSt=0,yplusCeil=0,deltaCeil=0;
    scalarField coefficientRho(rho.primitiveField());
    rho.oldTime();forAll(Y,bi)Y[bi].oldTime();
    auto exportFaces=[&](){
        OFstream out(runTime.path()/"diagnostics"/("ceilingFaces_"+runTime.timeName()+".csv"));out.precision(14);
        out<<"face_id,patch_id,patch_face,x_m,y_m,z_m,area_m2,owner,diameter_um,J_kg_m2_s,M_kg_m2"<<endl;
        forAll(mesh.boundary(),pi){const auto& pp=mesh.boundary()[pi];forAll(pp,fi)if(ceiling[pi][fi])forAll(names,bi){const vector q=pp.Cf()[fi];out<<pp.start()+fi<<','<<pi<<','<<fi<<','<<q.x()<<','<<q.y()<<','<<q.z()<<','<<pp.magSf()[fi]<<','<<pp.faceCells()[fi]<<','<<ds[bi]*1e6<<','<<J[bi].boundaryField()[pi][fi]<<','<<M[bi].boundaryField()[pi][fi]<<endl;}}
    };
    Info<<"Ceiling area "<<ceilingArea<<" m2; replay uses immutable frame files"<<endl;
    while(runTime.value()<runTime.endTime().value()-1e-9){
        const scalar old=runTime.value();runTime.setDeltaT(min(runTime.deltaTValue(),runTime.endTime().value()-old));++runTime;
        const scalar dt=runTime.deltaTValue(),now=runTime.value();
        label idx=0;while(idx+2<frameTimes.size()&&now>frameTimes[idx+1]+1e-9)idx++;
        if(idx!=bracket)load(idx);
        const scalar width=frameTimes[idx+1]-frameTimes[idx];
        const scalar w=max(scalar(0),min(scalar(1),(now-frameTimes[idx])/width));
        U=(1-w)*ulo()+w*uhi();T=(1-w)*tlo()+w*thi();p=(1-w)*plo()+w*phiP();k=(1-w)*klo()+w*khi();nut=(1-w)*nlo()+w*nhi();phi=(1-w)*flo()+w*fhi();
        rho=p/(dimensionedScalar("R",sqr(dimVelocity)/dimTemperature,RR)*T);
        const scalar sourceFrac=max(scalar(0),min(now,stop)-old)/dt;
        if(now>=nextCoeff-1e-9){
            nextCoeff=now+updateEvery;coefficientRho=rho.primitiveField();maxSt=0;yplusCeil=0;deltaCeil=0;
            volVectorField gradT(fvc::grad(T));
            volTensorField gradU(fvc::grad(U));
            volVectorField accel((uhi()-ulo())/dimensionedScalar("width",dimTime,width)+(U&gradU));
            scalarField mu(mesh.nCells()),kgas(mesh.nCells()),nu(mesh.nCells());
            forAll(mu,ci){mu[ci]=1.458e-6*Foam::pow(T[ci],1.5)/(T[ci]+110.4);kgas[ci]=mu[ci]*1006/.7;nu[ci]=mu[ci]/rho[ci];}
            forAll(names,bi){
                const scalar d=ds[bi];scalarField DB(mesh.nCells()),tau(mesh.nCells()),kth(mesh.nCells());
                sink[bi]=dimensionedScalar("zero",dimDensity/dimTime,0);
                forAll(DB,ci){
                    const scalar lam=mu[ci]/p[ci]*Foam::sqrt(constant::mathematical::pi*RR*T[ci]/2),Kn=2*lam/d;
                    const scalar cc=1+Kn*(1.257+.4*Foam::exp(-1.1/Kn));DB[ci]=kB*T[ci]*cc/(3*constant::mathematical::pi*mu[ci]*d);tau[ci]=rhoOil*sqr(d)*cc/(18*mu[ci]);
                    const scalar kr=kgas[ci]/.17;kth[ci]=2*1.17*cc*(kr+2.18*Kn)/((1+3*1.14*Kn)*(1+2*kr+4.36*Kn));
                    const scalar st=tau[ci]*mag(gradU[ci]);maxSt=max(maxSt,st);
                    drift[bi][ci]=tau[ci]*(1-rho[ci]/rhoOil)*vector(0,0,-9.81)-thermophoresis*kth[ci]*nu[ci]/T[ci]*gradT[ci]-inertia*tau[ci]*accel[ci];
                    D[bi][ci]=rho[ci]*(DB[ci]+max(nut[ci],scalar(0))/ScT);
                }
                forAll(mesh.boundary(),pi){const auto& pp=mesh.boundary()[pi];forAll(pp,fi){
                    const label ci=pp.faceCells()[fi];drift[bi].boundaryFieldRef()[pi][fi]=drift[bi][ci];D[bi].boundaryFieldRef()[pi][fi]=D[bi][ci];
                    if(!isA<wallFvPatch>(pp))continue;
                    const scalar dist=1/pp.deltaCoeffs()[fi];const vector nf=pp.Sf()[fi]/pp.magSf()[fi];
                    const vector Ut=U[ci]-(U[ci]&nf)*nf;
                    const scalar us=Foam::sqrt((nu[ci]+max(nut.boundaryField()[pi][fi],scalar(0)))*mag(Ut)/dist);
                    if(ceiling[pi][fi]){yplusCeil=max(yplusCeil,dist*us/nu[ci]);deltaCeil=max(deltaCeil,dist);}
                    scalar res=0,ya=0,invDa=1/DB[ci];
                    for(label j=0;j<nquad;j++){
                        const scalar yb=dist*Foam::pow(10.,-9.+9.*(j+1)/nquad),damp=1-Foam::exp(-yb*us/(26*nu[ci]));
                        const scalar invDb=1/(DB[ci]+.41*us*yb*sqr(damp)/ScT);
                        res+=(yb-ya)*.5*(invDa+invDb);ya=yb;invDa=invDb;
                    }
                    const scalar tw=T.boundaryField()[pi][fi]+(ceiling[pi][fi]?wallDT:0);
                    const scalar grav=tau[ci]*(1-rho[ci]/rhoOil)*(vector(0,0,-9.81)&nf);
                    const scalar therm=-thermophoresis*kth[ci]*nu[ci]/T[ci]*(tw-T[ci])/dist;
                    const scalar inert=-inertia*tau[ci]*(accel[ci]&nf);
                    // Approximate wall-normal turbulence variance declines over delta.
                    // This isotropic low-Stokes closure requires sensitivity analysis.
                    const scalar turbo=turbophoresis*tau[ci]*(2./3.)*max(k[ci],scalar(0))/dist;
                    // Wood-type turbulent inertia transfer, NOT a full trajectory model.
                    // Its vertical-pipe calibration is an important model uncertainty.
                    const scalar tauPlus=tau[ci]*sqr(us)/nu[ci];
                    const scalar vimp=impaction*us*min(scalar(.13),4.5e-4*sqr(tauPlus));
                    const scalar wallSpeed=(impactionMode=="parallel")?conductance(grav+therm+inert+turbo,res)+vimp:conductance(grav+therm+inert+turbo,1/(1/res+vimp));
                    const scalar vd=(ceiling[pi][fi]?sticking:otherSticking)*(legacy?(1/res+max(grav+therm,scalar(0))):wallSpeed);
                    vdep[bi][pi][fi]=vd;sink[bi][ci]+=rho[ci]*vd*pp.magSf()[fi]/mesh.V()[ci];
                }}
            }
        }
        scalar Q=0;forAll(mesh.boundary(),pi)if(mesh.boundary()[pi].name()=="hood")Q=gSum(phi.boundaryField()[pi]/rho.boundaryField()[pi]);
        const scalar airResidual=gSum(mag((fvc::ddt(rho)+fvc::div(phi))().primitiveField())*mesh.V());
        forAll(names,bi){
            volScalarField& y=Y[bi];
            surfaceScalarField aerosolPhi("aerosolPhi",phi+fvc::interpolate(rho)*(fvc::interpolate(drift[bi])&mesh.Sf()));
            forAll(mesh.boundary(),pi)if(isA<wallFvPatch>(mesh.boundary()[pi]))aerosolPhi.boundaryFieldRef()[pi]=0;
            y.correctBoundaryConditions();
            volScalarField effectiveSink("effectiveSink",sink[bi]);
            forAll(effectiveSink,ci)effectiveSink[ci]*=rho[ci]/coefficientRho[ci];
            #ifdef CEILING_CLEARANCE_INTEGRAL
            // Preserve the physical initial inventory in oldTime; start the
            // resolvent solve from zero so solver normalization is meaningful
            // even when the pseudo-time scale is very large.
            y.oldTime();y.primitiveFieldRef()=0;y.correctBoundaryConditions();
            #endif
            fvScalarMatrix eq(fvm::ddt(rho,y)+fvm::div(aerosolPhi,y,"div(aerosolPhi,Y)")-fvm::laplacian(D[bi],y)+fvm::Sp(effectiveSink,y)==dimensionedScalar("rate",dimMass/dimTime,emit*sourceFrac)*shape);
            eq.solve();y.correctBoundaryConditions();const surfaceScalarField flux(eq.flux());
            released[bi]+=emit*sourceFrac*dt;
            scalar ceilRate=0;
            forAll(mesh.boundary(),pi){const auto& pp=mesh.boundary()[pi];if(isA<wallFvPatch>(pp)){
                forAll(pp,fi){const label ci=pp.faceCells()[fi];const scalar jj=rho[ci]*y[ci]*vdep[bi][pi][fi];J[bi].boundaryFieldRef()[pi][fi]=jj;M[bi].boundaryFieldRef()[pi][fi]+=dt*jj;const scalar dm=dt*jj*pp.magSf()[fi];deposited[bi]+=dm;if(ceiling[pi][fi]){ceilingMass[bi]+=dm;ceilRate+=dm/dt;}}
            }else{const scalar out=dt*gSum(flux.boundaryField()[pi]);if(pp.name()=="hood")hood[bi]+=out;else if(pp.name()=="left")left[bi]+=out;else if(pp.name()=="right")right[bi]+=out;}}
            const scalar stock=gSum(rho.primitiveField()*y.primitiveField()*mesh.V());
            const scalar error=stock+hood[bi]+left[bi]+right[bi]+deposited[bi]-released[bi]-initial[bi];
            budget<<now<<','<<ds[bi]*1e6<<','<<initial[bi]<<','<<released[bi]<<','<<stock<<','<<hood[bi]<<','<<left[bi]<<','<<right[bi]<<','<<ceilingMass[bi]<<','<<deposited[bi]-ceilingMass[bi]<<','<<error<<','<<ceilRate<<','<<Q<<','<<airResidual<<','<<maxSt<<','<<yplusCeil<<','<<deltaCeil<<','<<gMin(y.primitiveField())<<endl;
        }
        if(runTime.writeTime()){runTime.write();exportFaces();Info<<"Saved ceiling deposition at "<<now<<" s"<<endl;}
    }
    runTime.writeNow();exportFaces();Info<<"End replay"<<endl;return 0;
}
