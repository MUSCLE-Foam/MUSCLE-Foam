#include "neutronicsSolver.H"
#include "localEulerDdtScheme.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neutronicsSolver, 0);
    defineRunTimeSelectionTable(neutronicsSolver, fvMesh);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::neutronicsSolver::writeData(Ostream&) const
{
    NotImplemented;
    return false;
}

bool Foam::neutronicsSolver::criticalityCalculationActive()
{
    return
    (
        criticalityCalculationActive_
        && (criticalityInactiveCounter_ >= criticalityInactive_)
    );
}

void Foam::neutronicsSolver::setDNPTransport()
{
    IOdictionary delayedNeutronConstants
    (
        IOobject
        (
            "delayedNeutronConstants",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    DNP_.set(new DNPTransport(*this, delayedNeutronConstants));
}

void Foam::neutronicsSolver::setDHPTransport()
{
    if (!decayHeat_)
    {
        FatalErrorInFunction
            << "Decay heat transport cannot be set"
            << exit(FatalError);
    }
    else
    {
        IOdictionary decayHeatConstants
        (
            IOobject
            (
                "decayHeatConstants",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        DHP_.set(new DHPTransport(*this, decayHeatConstants));
    }
}

void Foam::neutronicsSolver::setNeutronSources()
{
    forAll(FNS_, energyGroup)
    {
        FNS_.set
        (
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "FNS"+Foam::name(energyGroup),
                    mesh.time().name(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(inv(dimVolume*dimTime), Zero)
            )
        );
    }

    forAll(SNS_, energyGroup)
    {
        SNS_.set
        (
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "SNS"+Foam::name(energyGroup),
                    mesh.time().name(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(inv(dimVolume*dimTime), Zero)
            )
        );
    }
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::neutronicsSolver::read()
{
    return true;
}

void Foam::neutronicsSolver::normalize(scalar normalization)
{
    // Fluxes are normalized in the specific neutronics solver class

    forAll(FNS_, e)
    {
        FNS_[e] /= normalization;
        SNS_[e] /= normalization;
    }

    neutronFissionSource_ /= normalization;

    promptHeatSource_ /= normalization;

    DNP_->normalize(normalization);

    if (decayHeat_)
    {
        DHP_->normalize(normalization);
    }

    powerDensity_ /= normalization;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::neutronicsSolver::neutronicsSolver(fvMesh& mesh)
:
    regIOobject
    (
        IOobject
        (
            typeName,
            mesh.time().name(),
            mesh
        )
    ),

    mesh_(mesh),
    steady(mesh_.schemes().steady()),
    LTS(fv::localEulerDdt::enabled(mesh)),

    fvModelsPtr(nullptr),
    fvConstraintsPtr(nullptr),

    fluid_(mesh),

    DNP_(),
    decayHeat_(mesh.time().controlDict().lookupOrDefault("decayHeat", false)),
    DHP_(),
    energyGroups_(mesh.time().controlDict().lookup<label>("energyGroups", 6)),
    k_(mesh.time().controlDict().lookupOrDefault("keff", 1.0)),
    criticalityInactive_
    (
        mesh.time().controlDict().lookupOrDefault("criticalityInactive", 10)
    ),
    criticalityInactiveCounter_(0),
    criticalityCalculationActive_
    (
        mesh.time().controlDict().lookupOrDefault
        (
            "criticalityCalculationActive",
            false
        )
    ),

    nominalPower_
    (
        mesh.time().controlDict().lookupOrDefault("nominalPower", 1e9)
    ),

    neutronFissionSource_
    (
        IOobject
        (
            "neutronFissionSource",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(inv(dimVolume*dimTime), Zero)
    ),
    powerSource_
    (
        IOobject
        (
            "powerSource",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPower/dimVolume, Zero)
    ),
    promptHeatSource_
    (
        IOobject
        (
            "promptHeatSource",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPower/dimVolume, Zero)
    ),
    powerDensity_
    (
        IOobject
        (
            "powerDensity",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPower/dimVolume, Zero)
    ),
    rhoXS_
    (
        IOobject
        (
            "rhoXS",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    ),
    TXS_
    (
        IOobject
        (
            "TXS",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    ),

    FNS_(energyGroups_),
    SNS_(energyGroups_),

    mesh(mesh_),
    runTime(mesh_.time()),

    k(k_),

    fluid(fluid_)
{
    setDNPTransport();

    if (decayHeat_)
    {
        setDHPTransport();
    }

    setNeutronSources();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::neutronicsSolver::~neutronicsSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::fvModels& Foam::neutronicsSolver::fvModels() const
{
    if (!fvModelsPtr)
    {
        fvModelsPtr = &Foam::fvModels::New(mesh);
    }

    return *fvModelsPtr;
}


Foam::fvConstraints& Foam::neutronicsSolver::fvConstraints() const
{
    if (!fvConstraintsPtr)
    {
        fvConstraintsPtr = &Foam::fvConstraints::New(mesh);
    }

    return *fvConstraintsPtr;
}


void Foam::neutronicsSolver::precursorTransport()
{
    DNP_->solve();

    if (decayHeat_)
    {
        DHP_->solve();
    }
}


void Foam::neutronicsSolver::criticalityCalculation()
{
    scalar change = 0;

    if (criticalityCalculationActive())
    {
        change = k_;

        k_ = k_ *
            (
                fvc::domainIntegrate(neutronFissionSource_).value()
                / fvc::domainIntegrate(neutronFissionSource_.oldTime()).value()
            );

        k_ = min(k_,3.0);
        k_ = max(k_,0.3);

        change = 1e5 * (k_ - change) / k_;
    }

    Info << "Effective multiplication factor k = " << k_
             << ", change = " << change << " pcm"<< nl << endl;

    criticalityInactiveCounter_++;
}


void Foam::neutronicsSolver::updatePower()
{
    powerDensity_ = promptHeatSource_;

    if (DHP_.valid())
    {
        powerDensity_ += DHP_->decayHeatSource();
    }

    scalar totalPower = fvc::domainIntegrate(powerDensity_).value();

    if (criticalityCalculationActive_)
    {
        scalar normalization = max(totalPower, 1e-10) / nominalPower_;

        this->normalize(normalization);
    }
    else
    {
        Info << "Total power: " << totalPower << " W" << endl;
    }
}


// ************************************************************************* //