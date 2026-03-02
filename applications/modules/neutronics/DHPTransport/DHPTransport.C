#include "DHPTransport.H"
#include "neutronicsSolver.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DHPTransport, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::DHPTransport::readDecayHeatConstants(IOdictionary& dict)
{
    forAll(dec_, d)
    {
        dec_.set
        (
            d,
            new volScalarField
            (
                IOobject
                (
                    "dec"+Foam::name(d),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );

        lambda_[d].dimensions().reset(inv(dimTime));
        lambda_[d].value() = readScalar(dict.lookup("lambda"+Foam::name(d)));

        beta_[d] = readScalar(dict.lookup("beta"+Foam::name(d)));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DHPTransport::DHPTransport(neutronicsSolver& solver, IOdictionary& dict)
:
    mesh_(solver.mesh),

    neutronicsSolver_(solver),

    DHPGroups_(dict.lookup<label>("DHPGroups")),

    Sct_(dict.lookupOrDefault("Sct", 0.85)),
    Sc_(dict.lookupOrDefault("Sc", 1e10)),

    decayHeatSource_
    (
        IOobject
        (
            "decayHeatSource",
            mesh_.time().name(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimForce/(dimArea*dimTime), Zero)
    ),
    dec_(DHPGroups_),
    lambda_(DHPGroups_),
    beta_(DHPGroups_)
{
    readDecayHeatConstants(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::DHPTransport::~DHPTransport()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::DHPTransport::solve()
{
    volScalarField D = neutronicsSolver_.fluid.mu()
        / neutronicsSolver_.fluid.rho() / Sc_
        + neutronicsSolver_.fluid.nut() / Sct_;

    forAll(dec_, d)
    {
        fvScalarMatrix decEqn
        (
              fvm::ddt(dec_[d])
            + fvm::div(neutronicsSolver_.fluid.phi(), dec_[d])
            ==
              fvm::laplacian(D, dec_[d])
            - fvm::Sp(lambda_[d], dec_[d])
            + neutronicsSolver_.powerSource() * beta_[d]
        );

        decEqn.relax();
        decEqn.solve();

        dec_[d].correctBoundaryConditions();
    }

    updateSourceTerms();
}

void Foam::DHPTransport::normalize(scalar normalization)
{
    decayHeatSource_ /= normalization;
}

void Foam::DHPTransport::updateSourceTerms()
{
    decayHeatSource_ = Zero;

    forAll(dec_, d)
    {
        decayHeatSource_ += dec_[d] * lambda_[d];
    }
}


// ************************************************************************* //