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
    // Aerosols are rigorously dilute/one-way with respect to this air solve.
    // Solve them separately in ceilingReplayFoam on these saved native fluxes.
    // Existing Yp fields are placeholders, NOT calculated aerosol output here.
    if(Pstream::master())
        aerosolCsv()<<runTime.value()<<",0.1,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,"<<Qhood<<','<<gMax(temp.primitiveField())<<','<<airMass<<','<<heatFraction*heatEach*burners.size()<<','<<heatAdv<<','<<heatWall<<','<<airEnergy<<",nan,"<<gMax(mag(U.primitiveField()))<<','<<airMassError<<','<<energyError<<','<<hoodPressureDrop<<','<<heatOpen<<endl;
    aerosolState.set("released",released);aerosolState.set("hoodOut",hoodOut);
    aerosolState.set("leftOut",leftOut);aerosolState.set("rightOut",rightOut);
    aerosolState.set("deposited",deposited);aerosolState.set("initialMass",startMass);
    aerosolState.set("hoodPressureDrop",hoodPressureDrop);
    previousAirEnergy=airEnergy;previousKinEnergy=kinEnergy;previousPIntegral=pIntegral;previousAirMass=airMass;
}
