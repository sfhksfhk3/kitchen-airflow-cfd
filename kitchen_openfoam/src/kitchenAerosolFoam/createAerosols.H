// Additional dilute sectional aerosol model, GPL-3.0-or-later.
// All five Y fields are separate monodisperse unit-response experiments.
IOdictionary aerosolDict(IOobject("aerosolProperties",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE));
const bool projectInitialFlux=aerosolDict.getOrDefault<bool>("projectInitialMassFlux",false);
volScalarField ambientHydrostaticPressure(IOobject("ph_rgh",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedScalar("pAmbient",dimPressure,101325));
const scalarList diameters(aerosolDict.lookup("diameters"));
const wordList aerosolNames(aerosolDict.lookup("fieldNames"));
const scalar particleRho=aerosolDict.get<scalar>("particleDensity");
const scalar ScT=aerosolDict.get<scalar>("SchmidtT");
const scalar thermoCoeff=aerosolDict.get<scalar>("thermophoreticCoefficient");
const word thermophoresisModel=aerosolDict.getOrDefault<word>("thermophoresisModel","Talbot");
const scalar particleK=aerosolDict.getOrDefault<scalar>("particleThermalConductivity",.17);
const scalar targetQ=aerosolDict.getOrDefault<scalar>("hoodTargetQ",.5);
const scalar pressureGain=aerosolDict.getOrDefault<scalar>("pressureControllerGain",20);
const label hoodPatch=mesh.boundaryMesh().findPatchID("hood");
const scalar sticking=aerosolDict.get<scalar>("stickingProbability");
const scalar emissionRate=aerosolDict.get<scalar>("emissionKgPerSecPerBin");
const scalar sourceStop=aerosolDict.get<scalar>("sourceStopTime");
const scalar heatStop=aerosolDict.get<scalar>("heatStopTime");
const scalar heatEach=aerosolDict.get<scalar>("heatEachW");
const List<vector> burners(aerosolDict.lookup("burners"));
const scalar sigXY=aerosolDict.get<scalar>("sourceSigmaXY");
const scalar sigZ=aerosolDict.get<scalar>("sourceSigmaZ");
volScalarField sourceShape(IOobject("sourceShape",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),mesh,dimensionedScalar("zero",dimless/dimVolume,0));
forAll(burners,bi)
{
    scalarField w(mesh.nCells(),Zero);
    forAll(w,ci)
    {
        const vector d=mesh.C()[ci]-burners[bi];
        if(d.z()>0 && d.z()<.25 && mesh.C()[ci].y()<.6)
            w[ci]=Foam::exp(-.5*(sqr(d.x()/sigXY)+sqr(d.y()/sigXY)+sqr((d.z()-.07)/sigZ)));
    }
    const scalar sw=gSum(w*mesh.V());
    if(sw<=VSMALL) FatalErrorInFunction << "Empty burner source" << exit(FatalError);
    sourceShape.primitiveFieldRef()+=w/(sw*burners.size());
}
volScalarField heatSource(IOobject("heatSource",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),sourceShape*dimensionedScalar("P",dimPower,heatEach*burners.size()));
volScalarField aerosolSource(IOobject("aerosolSource",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),sourceShape*dimensionedScalar("rate",dimMass/dimTime,emissionRate));
PtrList<volScalarField> aerosolY(aerosolNames.size());
forAll(aerosolY,i)
{
    aerosolY.set(i,new volScalarField(IOobject(aerosolNames[i],runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh));
    mesh.setFluxRequired(aerosolNames[i]);
}
// State persists across restarts so cumulative mass accounting remains meaningful.
IOdictionary aerosolState(IOobject("aerosolState",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE));
scalarList released(aerosolY.size(),0),hoodOut(aerosolY.size(),0),leftOut(aerosolY.size(),0),rightOut(aerosolY.size(),0),deposited(aerosolY.size(),0),startMass(aerosolY.size(),0);
forAll(aerosolY,i) startMass[i]=gSum(rho.primitiveField()*aerosolY[i].primitiveField()*mesh.V());
if(aerosolState.found("released"))
{
    aerosolState.lookup("released")>>released;
    aerosolState.lookup("hoodOut")>>hoodOut;
    aerosolState.lookup("leftOut")>>leftOut;
    aerosolState.lookup("rightOut")>>rightOut;
    aerosolState.lookup("deposited")>>deposited;
    aerosolState.lookup("initialMass")>>startMass;
}
scalar previousAirEnergy=gSum(rho.primitiveField()*(thermo.T().primitiveField()-298.15)*1006*mesh.V());
scalar previousKinEnergy=gSum(rho.primitiveField()*.5*magSqr(U.primitiveField())*mesh.V());
scalar previousPIntegral=gSum(p.primitiveField()*mesh.V());
scalar previousAirMass=gSum(rho.primitiveField()*mesh.V());
scalar hoodPressureDrop=aerosolState.getOrDefault<scalar>("hoodPressureDrop",aerosolDict.getOrDefault<scalar>("initialHoodPressureDrop",2.));
autoPtr<OFstream> aerosolCsv;
if(Pstream::master())
{
    mkDir(runTime.path()/"diagnostics");
    aerosolCsv.reset(new OFstream(runTime.path()/"diagnostics"/("aerosol_"+runTime.timeName()+".csv")));
    aerosolCsv().precision(12);
    aerosolCsv()<<"time_s,diameter_um,stock_kg,released_kg,hood_out_kg,left_out_kg,right_out_kg,deposited_kg,mass_error_kg,walkway_c_over_R_s_m3,cmax_over_R_s_m3,min_Y,Qhood_m3_s,Tmax_K,air_mass_kg,heat_input_W,heat_advective_out_W,heat_wall_out_W,air_sensible_J,yplus_proxy_max,Umax_m_s,air_mass_residual_kg_s,energy_balance_W,hood_pressure_drop_Pa,heat_open_diffusive_W"<<endl;
}
Info<< "Normalized source integral "<<fvc::domainIntegrate(sourceShape).value()<<"; heat input "<<fvc::domainIntegrate(heatSource).value()<<" W"<<endl;
const scalar projectionStartTime=aerosolState.getOrDefault<scalar>("projectionStartTime",runTime.value());
if(projectInitialFlux && !aerosolState.getOrDefault<bool>("fluxProjectionCompleted",false))
{
    // Regridding a cell velocity does not automatically map its conservative
    // face flux. Project once before rhoEqn so the initial density stays finite.
    wordList bc(mesh.boundary().size(),"zeroGradient");
    forAll(mesh.boundary(),pi)if(!isA<wallFvPatch>(mesh.boundary()[pi]))bc[pi]="fixedValue";
    volScalarField correction(IOobject("fluxCorrection",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedScalar("zero",phi.dimensions()/dimLength,0),bc);
    mesh.setFluxRequired(correction.name());
    const surfaceScalarField initialPhi("initialPhi",phi);
    for(label n=0;n<12;n++)
    {
        fvScalarMatrix pe(fvm::laplacian(correction)==fvc::div(initialPhi));
        pe.solve();
        if(n==11)phi=initialPhi-pe.flux();
    }
    U+=fvc::reconstruct((phi-initialPhi)/fvc::interpolate(rho));
    U.correctBoundaryConditions();
    K=.5*magSqr(U);
    Info<<"Initial projected mass divergence max "<<gMax(mag(fvc::div(phi)().primitiveField()))<<endl;
    aerosolState.set("fluxProjectionCompleted",true);
    aerosolState.set("projectionStartTime",projectionStartTime);
}
