// Conservative Eulerian drift-diffusion for small-inertia, dilute oil droplets.
// Brownian diffusion (Cunningham), Stokes-Cunningham settling, assumed
// thermophoresis. Wall sink uses a resolved-coefficient 1D resistance integral,
// with van-Driest damped eddy diffusivity; perfect sticking is an assumption.
{
    const volScalarField& temp=thermo.T();
    const volScalarField mu(thermo.mu());
    const volScalarField kgas(thermo.kappa());
    const volScalarField nut(turbulence->nut());
    const volScalarField tk(turbulence->k());
    const volVectorField gradT(fvc::grad(temp));
    const scalar dt=runTime.deltaTValue();
    Info<< "Air maxima: U " << gMax(mag(U.primitiveField())) << " T " << gMax(temp.primitiveField()) << endl;
    const scalar kB=1.380649e-23;
    scalar heatAdv=0,heatWall=0,heatOpen=0,heatKin=0,netAirFlux=0,Qhood=0,walkVol=0,yplusMax=0;
    scalar airEnergy=gSum(rho.primitiveField()*(temp.primitiveField()-298.15)*1006*mesh.V());
    const scalar kinEnergy=gSum(rho.primitiveField()*.5*magSqr(U.primitiveField())*mesh.V());
    const scalar pIntegral=gSum(p.primitiveField()*mesh.V());
    const scalar airMass=gSum(rho.primitiveField()*mesh.V());
    const scalar gravityPower=gSum(rho.primitiveField()*(U.primitiveField()&g.value())*mesh.V());
    const volScalarField alphaEff(turbulence->alphaEff());
    forAll(mesh.boundary(),pi)
    {
        const auto& pp=mesh.boundary()[pi];
        if(pp.name()=="hood")Qhood=gSum(phi.boundaryField()[pi]/rho.boundaryField()[pi]);
        if(isA<wallFvPatch>(pp))
            heatWall+=gSum(-alphaEff.boundaryField()[pi]*thermo.he().boundaryField()[pi].snGrad()*pp.magSf());
        else {
            heatAdv+=gSum(phi.boundaryField()[pi]*1006*(temp.boundaryField()[pi]-298.15));
            heatKin+=gSum(phi.boundaryField()[pi]*.5*magSqr(U.boundaryField()[pi]));
            heatOpen+=gSum(-alphaEff.boundaryField()[pi]*thermo.he().boundaryField()[pi].snGrad()*pp.magSf());
        }
        netAirFlux+=gSum(phi.boundaryField()[pi]);
    }
    forAll(mesh.V(),ci) if(mesh.C()[ci].y()>.6 && mesh.C()[ci].z()>.9 && mesh.C()[ci].z()<1.8)walkVol+=mesh.V()[ci];
    reduce(walkVol,sumOp<scalar>());
    const scalar airMassError=(airMass-previousAirMass)/dt+netAirFlux;
    const scalar energyError=heatFraction*heatEach*burners.size()+gravityPower+(pIntegral-previousPIntegral)/dt-(airEnergy-previousAirEnergy+kinEnergy-previousKinEnergy)/dt-heatAdv-heatKin-heatWall-heatOpen;
    // Reuse the same wall turbulence profile across all five particle sizes.
    // This is algebraically identical to the per-bin computation below.
    static const scalarField wallFractions([](){scalarField a(72);for(label j=0;j<72;j++)a[j]=Foam::pow(10.,-9.+9.*(j+1)/72.);return a;}());
    List<scalarField> wallTurbulentDiff(mesh.boundary().size());
    forAll(mesh.boundary(),pi)
    {
        const auto& pp=mesh.boundary()[pi];
        if(!isA<wallFvPatch>(pp))continue;
        wallTurbulentDiff[pi].setSize(pp.size()*72);
        forAll(pp,fi)
        {
            const label ci=pp.faceCells()[fi];
            const scalar dist=1/pp.deltaCoeffs()[fi];
            const scalar nu=mu[ci]/rho[ci];
            const scalar ustar=Foam::pow(.09,.25)*Foam::sqrt(max(tk[ci],scalar(0)));
            for(label j=0;j<72;j++)
            {
                const scalar y=dist*wallFractions[j];
                const scalar damp=1-Foam::exp(-y*ustar/(26*nu));
                wallTurbulentDiff[pi][fi*72+j]=.41*ustar*y*sqr(damp)/ScT;
            }
        }
    }
    forAll(aerosolY,bi)
    {
        volScalarField& Y=aerosolY[bi];
        const scalar d=diameters[bi];
        volScalarField DB(IOobject("DB",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedScalar("D",dimArea/dimTime,0));
        volVectorField drift(IOobject("drift",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedVector("v",dimVelocity,Zero));
        volScalarField sink(IOobject("sink",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedScalar("rate",dimDensity/dimTime,0));
        forAll(Y,ci)
        {
            const scalar lam=mu[ci]/p[ci]*Foam::sqrt(constant::mathematical::pi*287.05*temp[ci]/2);
            const scalar Kn=2*lam/d;
            const scalar cc=1+Kn*(1.257+.4*Foam::exp(-1.1/Kn));
            DB[ci]=kB*temp[ci]*cc/(3*constant::mathematical::pi*mu[ci]*d);
            const scalar kr=kgas[ci]/particleK;
            const scalar kth=(thermophoresisModel=="Talbot") ? 2*1.17*cc*(kr+2.18*Kn)/((1+3*1.14*Kn)*(1+2*kr+2*2.18*Kn)) : thermoCoeff;
            drift[ci]=(particleRho-rho[ci])*sqr(d)*cc/(18*mu[ci])*g.value()
                      -kth*mu[ci]/(rho[ci]*temp[ci])*gradT[ci];
        }
        // Explicitly copy cell properties at patches, no fictitious zero slip there.
        forAll(mesh.boundary(),pi)
        {
            const auto& pp=mesh.boundary()[pi];
            forAll(pp,fi)
            {
                const label ci=pp.faceCells()[fi];
                DB.boundaryFieldRef()[pi][fi]=DB[ci];
                drift.boundaryFieldRef()[pi][fi]=drift[ci];
                if(!isA<wallFvPatch>(pp))continue;
                const scalar dist=1/pp.deltaCoeffs()[fi];
                const scalar nu=mu[ci]/rho[ci];
                const scalar ustar=Foam::pow(.09,.25)*Foam::sqrt(max(tk[ci],scalar(0)));
                yplusMax=max(yplusMax,dist*ustar/nu);
                // Integral dy / [DB + kappa*u_tau*y*(1-Foam::exp(-y+/26))^2/ScT]
                // 72 logarithmic intervals resolve the Brownian sublayer.
                scalar resistance=0,ya=0,invDa=1/DB[ci];
                for(label j=0;j<72;j++)
                {
                    const scalar yb=dist*wallFractions[j];
                    const scalar invDb=1/(DB[ci]+wallTurbulentDiff[pi][fi*72+j]);
                    resistance+=(yb-ya)*.5*(invDa+invDb);ya=yb;invDa=invDb;
                }
                const vector nf=pp.Sf()[fi]/pp.magSf()[fi];
                const scalar vdep=sticking*(1/resistance+max(scalar(0),drift[ci]&nf));
                sink[ci]+=rho[ci]*vdep*pp.magSf()[fi]/mesh.V()[ci];
            }
        }
        surfaceScalarField aerosolPhi("aerosolPhi",phi+fvc::interpolate(rho)*(fvc::interpolate(drift)&mesh.Sf()));
        // Wall transport is accounted for exactly once by the implicit sink.
        forAll(mesh.boundary(),pi)
            if(isA<wallFvPatch>(mesh.boundary()[pi]))aerosolPhi.boundaryFieldRef()[pi]=0;
        // inletOutlet switches according to particle flux, including settling.
        Y.correctBoundaryConditions();
        volScalarField diffusion("diffusion",rho*(DB+nut/ScT));
        fvScalarMatrix YE(fvm::ddt(rho,Y)+fvm::div(aerosolPhi,Y,"div(aerosolPhi,Y)")
                         -fvm::laplacian(diffusion,Y)+fvm::Sp(sink,Y)==sourceFraction*aerosolSource);
        YE.solve();
        Y.correctBoundaryConditions();
        // Use the equation's numerical boundary flux, not a separate approximation.
        const surfaceScalarField flux(YE.flux());
        released[bi]+=emissionRate*sourceFraction*dt;
        deposited[bi]+=dt*gSum(sink.primitiveField()*Y.primitiveField()*mesh.V());
        forAll(mesh.boundary(),pi)
        {
            const scalar out=dt*gSum(flux.boundaryField()[pi]);
            const word& n=mesh.boundary()[pi].name();
            if(n=="hood")hoodOut[bi]+=out;
            else if(n=="left")leftOut[bi]+=out;
            else if(n=="right")rightOut[bi]+=out;
        }
        const scalar stock=gSum(rho.primitiveField()*Y.primitiveField()*mesh.V());
        const scalar error=stock+hoodOut[bi]+leftOut[bi]+rightOut[bi]+deposited[bi]-released[bi]-startMass[bi];
        scalar walkStock=0;
        forAll(Y,ci) if(mesh.C()[ci].y()>.6 && mesh.C()[ci].z()>.9 && mesh.C()[ci].z()<1.8)walkStock+=rho[ci]*Y[ci]*mesh.V()[ci];
        reduce(walkStock,sumOp<scalar>());
        if(Pstream::master())
            aerosolCsv()<<runTime.value()<<','<<d*1e6<<','<<stock<<','<<released[bi]<<','<<hoodOut[bi]<<','<<leftOut[bi]<<','<<rightOut[bi]<<','<<deposited[bi]<<','<<error<<','<<walkStock/(walkVol*emissionRate)<<','<<gMax(rho.primitiveField()*Y.primitiveField())/emissionRate<<','<<gMin(Y.primitiveField())<<','<<Qhood<<','<<gMax(temp.primitiveField())<<','<<airMass<<','<<heatFraction*heatEach*burners.size()<<','<<heatAdv<<','<<heatWall<<','<<airEnergy<<','<<yplusMax<<','<<gMax(mag(U.primitiveField()))<<','<<airMassError<<','<<energyError<<','<<hoodPressureDrop<<','<<heatOpen<<endl;
    }
    aerosolState.set("released",released);aerosolState.set("hoodOut",hoodOut);
    aerosolState.set("leftOut",leftOut);aerosolState.set("rightOut",rightOut);
    aerosolState.set("deposited",deposited);aerosolState.set("initialMass",startMass);
    aerosolState.set("hoodPressureDrop",hoodPressureDrop);
    previousAirEnergy=airEnergy;previousKinEnergy=kinEnergy;previousPIntegral=pIntegral;previousAirMass=airMass;
}
