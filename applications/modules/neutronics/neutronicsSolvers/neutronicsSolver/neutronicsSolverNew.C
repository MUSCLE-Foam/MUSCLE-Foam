#include "neutronicsSolver.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::neutronicsSolver> Foam::neutronicsSolver::New
(
    const word& neutronicsSolverName,
    fvMesh& mesh
)
{
    Info<< "Selecting neutronics solver " << neutronicsSolverName << endl;

    if (!fvMeshConstructorTablePtr_)
    {
        FatalErrorInFunction
            << "neutronics solvers table is empty"
            << exit(FatalError);
    }

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(neutronicsSolverName);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown neutronics solver type "
            << neutronicsSolverName << nl << nl
                << "Valid neutronics solvers are :" << endl
                << fvMeshConstructorTablePtr_->sortedToc()
                << exit(FatalError);
    }

    autoPtr<neutronicsSolver> neutronicsSolverPtr(cstrIter()(mesh));

    // Ensure fvModels and fvConstraints are constructed
    // before time is incremented
    neutronicsSolverPtr->fvModels();
    neutronicsSolverPtr->fvConstraints();

    return neutronicsSolverPtr;
}


// ************************************************************************* //