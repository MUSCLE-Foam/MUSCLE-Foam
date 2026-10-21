// ************************************************************************* //

#include "reactivityInsertion.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(reactivityInsertion, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        reactivityInsertion,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::reactivityInsertion::reactivityInsertion
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    rho_(readScalar(dict.lookup("rho")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::reactivityInsertion::~reactivityInsertion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::reactivityInsertion::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    return true;
}


bool Foam::functionObjects::reactivityInsertion::execute()
{
    PKESolver& PKE =
        mesh_.lookupObjectRef<PKESolver>("neutronicsSolver");

    if (mesh_.time().value() >= 0.1)
    {
        PKE.reactivity() = rho_; // 0.1 * PKE.betaTotal();
    }

    return true;
}


bool Foam::functionObjects::reactivityInsertion::write()
{
    return true;
}


// ************************************************************************* //