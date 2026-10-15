    const volScalarField& temp=thermo.T();
    const volScalarField mu(thermo.mu());
    const volScalarField nut(turbulence->nut());
    const volScalarField tk(turbulence->k());
    const volVectorField gradT(fvc::grad(temp));

    Info<< "Air maxima: U " << gMax(mag(U.primitiveField())) << " T " << gMax(temp.primitiveField()) << endl;
    const scalar kB=1.380649e-23;
    scalar heatAdv=0,heatWall=0,heatKin=0,Qhood=0,walkVol=0,yplusMax=0;
    scalar worstMassError=0;
    const scalar dU=Foam::sqrt(gSum(mesh.V()*magSqr(U.primitiveField()-oldU))/(gSum(mesh.V()*magSqr(U.primitiveField()))+VSMALL));
    const scalar dT=Foam::sqrt(gSum(mesh.V()*sqr(temp.primitiveField()-oldT))/(gSum(mesh.V()*sqr(temp.primitiveField()-298.15))+VSMALL));
    const scalar gravityPower=gSum(rho.primitiveField()*(U.primitiveField()&g.value())*mesh.V());
    scalar airEnergy=gSum(rho.primitiveField()*(temp.primitiveField()-298.15)*1006*mesh.V());
    const volScalarField alphaEff(turbulence->alphaEff());
    forAll(mesh.boundary(),pi)
    {
        const auto& pp=mesh.boundary()[pi];
        if(pp.name()=="hood")Qhood=gSum(phi.boundaryField()[pi]/rho.boundaryField()[pi]);
        heatWall+=gSum(-alphaEff.boundaryField()[pi]*thermo.he().boundaryField()[pi].snGrad()*pp.magSf());
        if(!isA<wallFvPatch>(pp)) {
          heatAdv+=gSum(phi.boundaryField()[pi]*1006*(temp.boundaryField()[pi]-298.15));
          heatKin+=gSum(phi.boundaryField()[pi]*.5*magSqr(U.boundaryField()[pi]));
        }
    }
    forAll(mesh.V(),ci) if(mesh.C()[ci].y()>.6 && mesh.C()[ci].z()>.9 && mesh.C()[ci].z()<1.8)walkVol+=mesh.V()[ci];
    reduce(walkVol,sumOp<scalar>());
    forAll(aerosolY,bi)
    {
        volScalarField& Y=aerosolY[bi];
        const scalar d=diameters[bi];
        volScalarField DB(IOobject("DB",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedScalar("D",dimArea/dimTime,0));
        volVectorField drift(IOobject("drift",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedVector("v",dimVelocity,Zero));
        volScalarField sink(IOobject("sink",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedScalar("rate",dimDensity/dimTime,0));
        forAll(Y,ci)
        {
            const scalar lam=66e-9*(temp[ci]/298.15)*(101325/p[ci]);
            const scalar Kn=2*lam/d;
            const scalar cc=1+Kn*(1.257+.4*Foam::exp(-1.1/Kn));
            DB[ci]=kB*temp[ci]*cc/(3*constant::mathematical::pi*mu[ci]*d);
            drift[ci]=(particleRho-rho[ci])*sqr(d)*cc/(18*mu[ci])*g.value()
                      -thermoCoeff*mu[ci]/(rho[ci]*temp[ci])*gradT[ci];
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
                    const scalar yb=dist*Foam::pow(10.,-9.+9.*(j+1)/72.);
                    const scalar damp=1-Foam::exp(-yb*ustar/(26*nu));
                    const scalar invDb=1/(DB[ci]+.41*ustar*yb*sqr(damp)/ScT);
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
        fvScalarMatrix YE(fvm::div(aerosolPhi,Y,"div(aerosolPhi,Y)")
                         -fvm::laplacian(diffusion,Y)+fvm::Sp(sink,Y)==aerosolSource);
        YE.relax(.7);
        YE.solve();
        Y.correctBoundaryConditions();
        // Use the equation's numerical boundary flux, not a separate approximation.
        const surfaceScalarField flux(YE.flux());

        scalar hood=0,left=0,right=0;
        forAll(mesh.boundary(),pi)
        {
            const scalar f=gSum(flux.boundaryField()[pi]);
            const word& n=mesh.boundary()[pi].name();
            if(n=="hood")hood=f;
            if(n=="left")left=f;
            if(n=="right")right=f;
        }
        const scalar dep=gSum(sink.primitiveField()*Y.primitiveField()*mesh.V());
        const scalar stock=gSum(rho.primitiveField()*Y.primitiveField()*mesh.V());
        const scalar err=(hood+left+right+dep-emissionRate)/emissionRate;
        worstMassError=max(worstMassError,mag(err));
        scalar walkStock=0;
        forAll(Y,ci) if(mesh.C()[ci].y()>.6 && mesh.C()[ci].z()>.9 && mesh.C()[ci].z()<1.8)walkStock+=rho[ci]*Y[ci]*mesh.V()[ci];
        reduce(walkStock,sumOp<scalar>());
        const scalar eErr=(heatEach*burners.size()+gravityPower-heatAdv-heatKin-heatWall)/(heatEach*burners.size()+SMALL);
        steadyCsv<<runTime.value()<<','<<d*1e6<<','<<stock/emissionRate<<','<<hood/emissionRate<<','<<left/emissionRate<<','<<right/emissionRate<<','<<dep/emissionRate<<','<<err<<','<<walkStock/(walkVol*emissionRate)<<','<<gMax(rho.primitiveField()*Y.primitiveField())/emissionRate<<','<<gMin(Y.primitiveField())<<','<<Qhood<<','<<gMax(temp.primitiveField())<<','<<gMax(mag(U.primitiveField()))<<','<<heatAdv<<','<<heatWall<<','<<heatKin<<','<<gravityPower<<','<<eErr<<','<<dU<<','<<dT<<endl;
    }
    if(runTime.value()>200 && dU<1e-5 && dT<1e-5 && worstMassError<1e-5)steadyCount++;else steadyCount=0;
    Info<<"Steady changes U "<<dU<<" T "<<dT<<" aerosol balance "<<worstMassError<<" count "<<steadyCount<<endl;
