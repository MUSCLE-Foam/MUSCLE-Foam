#include "neutronicsHeatSource.H"
#include "basicThermo.H"
#include "fvModels.H"
#include "fvMatrix.H"
#include "Scale.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(neutronicsHeatSource, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        neutronicsHeatSource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::neutronicsHeatSource::neutronicsHeatSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    gamma_
    (
        dimPower/(dimVolume*dimTemperature),
        dict.lookupOrDefault<scalar>("gamma", 0)
    ),
    Tref_
    (
        dimTemperature,
        dict.lookupOrDefault<scalar>("Tref", 300)
    ),
    zone_(mesh, dict),
    heNames_()
{
    if (mesh.foundObject<basicThermo>(physicalProperties::typeName))
    {
        const basicThermo& thermo =
            mesh().lookupObject<basicThermo>(physicalProperties::typeName);

        heNames_.append(thermo.he().name());
    }
    else
    {
        forAll(mesh.objectRegistry::toc<basicThermo>(), i)
        {
            heNames_.append
            (
                mesh.lookupObject<basicThermo>
                (
                    mesh.objectRegistry::toc<basicThermo>()[i]
                ).he().name()
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::neutronicsHeatSource::~neutronicsHeatSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::neutronicsHeatSource::addSupFields() const
{
    return heNames_;
}


void Foam::fv::neutronicsHeatSource::addSup
(
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    const labelUList& cells = zone_.zone();

    scalarField& eqnSource = eqn.source();

    const volScalarField& powerDensity
    (
        mesh().lookupObject<volScalarField>("powerDensity")
    );

    const volScalarField& T(mesh().lookupObject<volScalarField>("T"));

    eqnSource -= powerDensity*mesh().V();

    forAll(cells, i)
    {
        eqnSource[cells[i]] += mesh().V()[cells[i]]
            *gamma_.value()*(T[cells[i]]-Tref_.value());
    }
}


void Foam::fv::neutronicsHeatSource::addSup
(
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    addSup(he, eqn);
}


void Foam::fv::neutronicsHeatSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    const labelList& cells = zone_.zone();

    scalarField& eqnSource = eqn.source();

    const volScalarField& powerDensity
    (
        mesh().lookupObject<volScalarField>("powerDensity")
    );

    const word phaseName(he.name().substr(2));

    const volScalarField& T
    (
        mesh().lookupObject<volScalarField>("T."+phaseName)
    );

    eqnSource -= alpha*powerDensity*mesh().V();

    forAll(cells, i)
    {
        eqnSource[cells[i]] += mesh().V()[cells[i]]
            *gamma_.value()*(T[cells[i]]-Tref_.value())
            *alpha[cells[i]];
    }
}


bool Foam::fv::neutronicsHeatSource::movePoints()
{
    return true;
}


void Foam::fv::neutronicsHeatSource::topoChange(const polyTopoChangeMap& map)
{}


void Foam::fv::neutronicsHeatSource::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::neutronicsHeatSource::distribute(const polyDistributionMap& map)
{}


bool Foam::fv::neutronicsHeatSource::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        gamma_.value() = dict.lookupOrDefault<scalar>("gamma", 0);
        Tref_.value() = dict.lookupOrDefault<scalar>("Tref", 300);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //