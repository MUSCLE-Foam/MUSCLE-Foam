#include "DNPTransport.H"
#include "neutronicsSolver.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DNPTransport, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::DNPTransport::readDelayedNeutronConstants(IOdictionary& dict)
{
    forAll(prec_, p)
    {
        prec_.set
        (
            p,
            new volScalarField
            (
                IOobject
                (
                    "prec"+Foam::name(p),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );

        lambda_[p].dimensions().reset(inv(dimTime));
        lambda_[p].value() = readScalar(dict.lookup("lambda"+Foam::name(p)));

        beta_[p] = readScalar(dict.lookup("beta"+Foam::name(p)));

        betaTotal_ += beta_[p];
    }

    updateSourceTerms();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DNPTransport::DNPTransport(neutronicsSolver& solver, IOdictionary& dict)
:
    mesh_(solver.mesh),

    neutronicsSolver_(solver),

    DNPGroups_(dict.lookup<label>("DNPGroups")),

    Sct_(dict.lookupOrDefault("Sct", 0.85)),
    Sc_(dict.lookupOrDefault("Sc", 1e10)),

    prec_(DNPGroups_),
    lambda_(DNPGroups_),
    beta_(DNPGroups_),
    betaTotal_(0),

    DNS_(neutronicsSolver_.energyGroups()),
    chi_(neutronicsSolver_.energyGroups())
{
    forAll(chi_, energyGroup)
    {
        chi_.set
        (
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "chiDel"+Foam::name(energyGroup),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );
    }

    forAll(DNS_, energyGroup)
    {
        DNS_.set
        (
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "DNS"+Foam::name(energyGroup),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(inv(dimVolume*dimTime), Zero)
            )
        );
    }

    readDelayedNeutronConstants(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::DNPTransport::~DNPTransport()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::DNPTransport::solve()
{
    volScalarField D = neutronicsSolver_.fluid.mu()
        / neutronicsSolver_.fluid.rho() / Sc_
        + neutronicsSolver_.fluid.nut() / Sct_;

    forAll(prec_, p)
    {
        fvScalarMatrix precEqn
        (
              fvm::ddt(prec_[p])
            + fvm::div(neutronicsSolver_.fluid.phi(), prec_[p])
            ==
              fvm::laplacian(D, prec_[p])
            - fvm::Sp(lambda_[p], prec_[p])
            + (1.0 / neutronicsSolver_.k)
                * neutronicsSolver_.neutronFissionSource()
                * beta_[p]
        );

        precEqn.relax();
        precEqn.solve();

        prec_[p].correctBoundaryConditions();
    }

    updateSourceTerms();
}


void Foam::DNPTransport::normalize(scalar normalization)
{
    forAll(DNS_, e)
    {
        DNS_[e] /= normalization;
    }
}


void Foam::DNPTransport::updateSourceTerms()
{
    forAll(DNS_, e)
    {
        DNS_[e] = Zero;

        forAll(prec_, p)
        {
            DNS_[e] += chi_[e] * prec_[p] * lambda_[p];
        }
    }
}


// ************************************************************************* //