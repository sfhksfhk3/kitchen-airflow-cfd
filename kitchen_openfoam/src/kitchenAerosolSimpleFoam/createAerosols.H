// Additional dilute sectional aerosol model, GPL-3.0-or-later.
// All five Y fields are separate monodisperse unit-response experiments.
IOdictionary aerosolDict(IOobject("aerosolProperties",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE));
const scalarList diameters(aerosolDict.lookup("diameters"));
const wordList aerosolNames(aerosolDict.lookup("fieldNames"));
const scalar particleRho=aerosolDict.get<scalar>("particleDensity");
const scalar ScT=aerosolDict.get<scalar>("SchmidtT");
const scalar thermoCoeff=aerosolDict.get<scalar>("thermophoreticCoefficient");
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

label steadyCount=0;
OFstream steadyCsv(runTime.path()/"steady_diagnostics.csv");
steadyCsv.precision(12);
steadyCsv<<"iteration,diameter_um,stock_over_R_s,hood_fraction,left_fraction,right_fraction,deposition_fraction,mass_relative_error,walkway_c_over_R_s_m3,cmax_over_R_s_m3,min_Y,Qhood_m3_s,Tmax_K,Umax_m_s,heat_advective_W,heat_diffusive_W,heat_kinetic_W,gravity_power_W,thermal_relative_error,relative_U_change,relative_T_change"<<endl;
