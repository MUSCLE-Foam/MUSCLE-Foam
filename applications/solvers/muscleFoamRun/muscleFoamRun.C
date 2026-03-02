/*---------------------------------------------------------------------------*\

         __  __ _   _ ___  ___ _    ___         ___
        |  \/  | | | / __|/ __| |  | __|  ___  | __|__  __ _ _ __
        | |\/| | |_| \__ \ (__| |__| _|  |___| | _/ _ \/ _` | '  \
        |_|  |_|\___/|___/\___|____|___|       |_|\___/\__,_|_|_|_|

    MUSCLE-Foam: an OpenFOAM-based multiphysics solver for reactor physics
    Copyright (C) 2026 NRG PALLAS

\-----------------------------------------------------------------------------/

License
    This file is part of MUSCLE-Foam. It is published under the terms of the
    GNU General Public License. See <http://www.gnu.org/licenses/>.

Application
    muscleFoamRun

Description
    MUSCLE-Foam is a modular neutronics solver for OpenFOAM-13. It allows for
    coupled simulations of standard OpenFOAM solvers with the addition of
    neutronics.

    Currently, the following neutronics solvers are supported
        - Diffusion (see also diffusionSolver.H/.C)
        - SP3 (see also SP3Solver.H/.C)
        - SN (see also SNSolver.H/.C)

    MUSCLE-Foam also allows for the optional transport of delayed neutron and
    decay heat precursor groups by molecular diffusion, and turbulent diffusion
    and advection in the case of fuel/coolant motion.

Usage
    \b muscleFoamRun [OPTION]

      - \par -solver <name>
        Solver name

      - \par -neutronicsSolver <name>
        Neutronics solver name

      - \par -libs '(\"lib1.so\" ... \"libN.so\")'
        Specify the additional libraries loaded

    Example usage:
      - To run a \c rhoPimpleFoam case by specifying the solver on the
        command line:
        \verbatim
            muscleFoamRun -solver fluid -neutronicsSolver diffusion
        \endverbatim

      - To update and run a \c rhoPimpleFoam case add the following entries to
        the controlDict:
        \verbatim
            application     muscleFoamRun;

            solver          fluid;

            neutronicsSolver diffusion;
        \endverbatim
        then execute \c muscleFoamRun

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "solver.H"
#include "neutronicsSolver.H"
#include "pimpleControl.H"
#include "pimpleSingleRegionControl.H"
#include "setDeltaT.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "solver",
        "name",
        "Solver name"
    );

    argList::addOption
    (
        "neutronicsSolver",
        "name",
        "Neutronics solver name"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    // Read external cycles from controlDict
    int externalCycles
    (
        runTime.controlDict().lookupOrDefault("externalCycles", 1)
    );

    // Read the solverName from the optional solver entry in controlDict
    word solverName
    (
        runTime.controlDict().lookupOrDefault("solver", word::null)
    );

    // Read the neutronicsSolverName from the optional solver entry in controlDict
    word neutronicsSolverName
    (
        runTime.controlDict().lookupOrDefault("neutronicsSolver", word::null)
    );

    // Optionally reset the solver name from the -solver command-line argument
    args.optionReadIfPresent("solver", solverName);

    // Optionally reset the neutronicsSolver name from the -solver command-line argument
    args.optionReadIfPresent("neutronicsSolver", neutronicsSolverName);

    // Check the solverName has been set
    if (solverName == word::null)
    {
        args.printUsage();

        FatalErrorIn(args.executable())
            << "solver not specified in the controlDict or on the command-line"
            << exit(FatalError);
    }
    else
    {
        // Load the solver library
        solver::load(solverName);
    }

    // Check the neutronicsSolverName has been set
    if (neutronicsSolverName == word::null)
    {
        args.printUsage();

        FatalErrorIn(args.executable())
            << "neutronics solver not specified in the controlDict or on the command-line"
            << exit(FatalError);
    }

    // Create the default single region mesh
    #include "createMesh.H"

    // Instantiate the selected solver
    autoPtr<solver> solverPtr(solver::New(solverName, mesh));
    solver& solver = solverPtr();

    // Instantiate the selected neutronics solver
    autoPtr<neutronicsSolver> neutronicsSolverPtr(neutronicsSolver::New(neutronicsSolverName, mesh));
    neutronicsSolver& neutronicsSolver = neutronicsSolverPtr();

    // Create the outer PIMPLE loop and control structure
    pimpleSingleRegionControl pimple(solver.pimple);

    // Create the neutronics 'NPIMPLE' loop
    pimpleControl npimple(mesh, "NPIMPLE");

    // Set the initial time-step
    setDeltaT(runTime, solver);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Starting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        solver.preSolve();

        // Adjust the time-step according to the solver maxDeltaT
        adjustDeltaT(runTime, solver);

        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // External cycle loop
        for (int counter = 0; counter < externalCycles; counter++)
        {
            if (externalCycles > 1)
            {
                Info<< "External cycle " << counter + 1 << nl << endl;
            }

            // PIMPLE corrector loop
            while (pimple.loop())
            {
                if (solver.pimple.flow())
                {
                    solver.moveMesh();
                    solver.motionCorrector();
                }

                if (solver.pimple.models())
                {
                    solver.fvModels().correct();
                }

                solver.prePredictor();

                if (solver.pimple.predictTransport())
                {
                    if (solver.pimple.flow())
                    {
                        solver.momentumTransportPredictor();
                    }

                    if (solver.pimple.thermophysics())
                    {
                        solver.thermophysicalTransportPredictor();
                    }
                }

                if (solver.pimple.flow())
                {
                    solver.momentumPredictor();
                }

                if (solver.pimple.thermophysics())
                {
                    solver.thermophysicalPredictor();
                }

                if (solver.pimple.flow())
                {
                    solver.pressureCorrector();
                }

                if (solver.pimple.correctTransport())
                {
                    if (solver.pimple.flow())
                    {
                        solver.momentumTransportCorrector();
                    }

                    if (solver.pimple.thermophysics())
                    {
                        solver.thermophysicalTransportCorrector();
                    }
                }
            }

            solver.postSolve();

            // NPIMPLE loop
            while (npimple.loop())
            {
                neutronicsSolver.solve();
                neutronicsSolver.precursorTransport();
                neutronicsSolver.criticalityCalculation();
                neutronicsSolver.updatePower();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //