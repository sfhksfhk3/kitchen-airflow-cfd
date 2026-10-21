/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    buoyantPimpleFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent flow of compressible fluids
    for ventilation and heat-transfer, with optional mesh motion
    and mesh topology changes.

    Turbulence is modelled using a run-time selectable compressible RAS or
    LES model.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallFvPatch.H"
#include "inletOutletFvPatchFields.H"
#include "dynamicFvMesh.H"
#include "rhoThermo.H"
#include "turbulentFluidThermoModel.H"
#include "radiationModel.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent fluid flow"
        " of compressible fluids, including radiation,"
        " with optional mesh motion and mesh topology changes."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createRhoUfIfPresent.H"
    #include "createAerosols.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (true)
    {
        // Time::run() normally accepts an endpoint within half a time step.
        // Take the remaining fraction explicitly so exported end fields have
        // the actual requested physical time, including after restarts.
        const scalar remainingBeforeRun=runTime.endTime().value()-runTime.value();
        if(runTime.stopAt()==Time::saEndTime && remainingBeforeRun>1e-9 && remainingBeforeRun<runTime.deltaTValue())
            runTime.setDeltaT(remainingBeforeRun,false);
        if(!runTime.run())break;
        #include "readDyMControls.H"

        // Store divrhoU from the previous mesh
        // so that it can be mapped and used in correctPhi
        // to ensure the corrected phi has the same divergence
        autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU.reset
            (
                new volScalarField
                (
                    "divrhoU",
                    fvc::div(fvc::absolute(phi, rho, U))
                )
            );
        }

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        // Anchor output ticks to absolute physical time after arbitrary restarts.
        const scalar nextOutput=10*(label((runTime.value()+1e-8)/10)+1);
        const scalar untilOutput=nextOutput-runTime.value();
        if(untilOutput>SMALL && untilOutput<runTime.deltaTValue())runTime.setDeltaT(untilOutput,false);
        const scalar untilEnd=runTime.endTime().value()-runTime.value();
        if(untilEnd>SMALL && untilEnd<runTime.deltaTValue())runTime.setDeltaT(untilEnd,false);
        if(projectInitialFlux && runTime.value()<projectionStartTime+.1)runTime.setDeltaT(min(runTime.deltaTValue(),scalar(.002)),false);

        ++runTime;
        const scalar tOld=runTime.value()-runTime.deltaTValue();
        const scalar sourceFraction=max(scalar(0),min(runTime.value(),sourceStop)-tOld)/runTime.deltaTValue();
        const scalar heatFraction=max(scalar(0),min(runTime.value(),heatStop)-tOld)/runTime.deltaTValue();
        // Uniform plenum suction leaves local outlet velocity unconstrained.
        // Slow feedback enforces the explicitly assumed network volume flow.
        const scalar measuredQ=gSum(phi.boundaryField()[hoodPatch]/rho.boundaryField()[hoodPatch]);
        if (!prescribedHoodFlow)
        {
            hoodPressureDrop=min(scalar(1300),max(scalar(.02),hoodPressureDrop+pressureGain*runTime.deltaTValue()*(targetQ-measuredQ)));
            p_rgh.boundaryFieldRef()[hoodPatch] == scalar(101325-hoodPressureDrop);
        }

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Store momentum to set rhoUf for introduced faces.
                autoPtr<volVectorField> rhoU;
                if (rhoUf.valid())
                {
                    rhoU.reset(new volVectorField("rhoU", rho*U));
                }

                // Do any mesh changes
                mesh.update();

                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & rhoUf();

                        #include "correctPhi.H"

                        // Make the fluxes relative to the mesh-motion
                        fvc::makeRelative(phi, rho, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            if (pimple.firstIter() && !pimple.SIMPLErho())
            {
                #include "rhoEqn.H"
            }

            #include "UEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        rho = thermo.rho();
        if(prescribedHoodFlow)
        {
            hoodPressureDrop=101325-gSum(p.boundaryField()[hoodPatch]*mesh.magSf().boundaryField()[hoodPatch])/gSum(mesh.magSf().boundaryField()[hoodPatch]);
        }

        #include "aerosolEqns.H"

        if(mag(runTime.value()-10*label(runTime.value()/10+.5))<1e-7)runTime.writeNow();
        else runTime.write();

        runTime.printExecutionTime(Info);
    }

    runTime.writeNow();
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
