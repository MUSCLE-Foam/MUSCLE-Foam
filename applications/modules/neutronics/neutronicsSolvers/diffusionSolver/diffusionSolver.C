#include "diffusionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diffusionSolver, 0);
    addToRunTimeSelectionTable(neutronicsSolver, diffusionSolver, fvMesh);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::diffusionSolver::updateCrossSections()
{
    volScalarField K(fluid.rho()/rhoXS_);
    volScalarField logT(log(fluid.T()/TXS_));

    forAll(flux_, ei)
    {
        TR_[ei] = (TR0_[ei] + alphaTR_[ei] * logT) * K;
        M_[ei]  = (M0_[ei]  + alphaM_[ei]  * logT) * K;
        SP_[ei] = (SP0_[ei] + alphaSP_[ei] * logT) * K;
        FN_[ei] = (FN0_[ei] + alphaFN_[ei] * logT) * K;

        F_[ei] = FN_[ei] * (1.0 - DNP_->betaTotal()) * chiPrompt_[ei];

        forAll(S_[ei], ej)
        {
            S_[ei][ej] = (S0_[ei][ej] + alphaS_[ei][ej] * logT) * K;
        }
    }
}


void Foam::diffusionSolver::updateSourceTerms()
{
    neutronFissionSource_ = Zero;
    powerSource_ = Zero;

    forAll(SNS_, e)
    {
        SNS_[e] = Zero;
        FNS_[e] = Zero;
    }

    forAll(flux_, ei)
    {
        neutronFissionSource_ += FN_[ei] * flux_[ei];
        powerSource_ += SP_[ei] * flux_[ei];

        forAll(SNS_, ej)
        {
            if (ei != ej)
            {
                SNS_[ei] += S_[ej][ei] * flux_[ej];
                FNS_[ei] += FN_[ej] * (1.0 - DNP_->betaTotal())
                    * flux_[ej] * chiPrompt_[ei];
            }
        }
    }

    promptHeatSource_ = (1.0 - DNP_->betaTotal()) * powerSource_;

    DNP_->updateSourceTerms();
}


void Foam::diffusionSolver::normalize(scalar normalization)
{
    forAll(flux_, e)
    {
        flux_[e] /= normalization;
    }

    this->neutronicsSolver::normalize(normalization);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusionSolver::diffusionSolver(fvMesh& mesh)
:
    neutronicsSolver(mesh),

    flux_(energyGroups_),

    fluxAlbedo_
    (
        IOobject
        (
            "fluxAlbedo",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    ),

    IV_(energyGroups_),
    chiPrompt_(energyGroups_),

    TR_(energyGroups_),
    TR0_(energyGroups_),
    alphaTR_(energyGroups_),

    M_(energyGroups_),
    M0_(energyGroups_),
    alphaM_(energyGroups_),

    SP_(energyGroups_),
    SP0_(energyGroups_),
    alphaSP_(energyGroups_),

    FN_(energyGroups_),
    FN0_(energyGroups_),
    alphaFN_(energyGroups_),

    S_(energyGroups_),
    S0_(energyGroups_),
    alphaS_(energyGroups_),

    F_(energyGroups_),

    DAlbedo_
    (
        IOobject
        (
            "DAlbedo",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    )
{
    init();
    updateCrossSections();
    updateSourceTerms();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusionSolver::~diffusionSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::diffusionSolver::solve()
{
    updateCrossSections();

    forAll(flux_, e)
    {
        DAlbedo_ = 1.0 / (3.0 * TR_[e]);

        fvScalarMatrix fluxEqn
        (
              fvm::ddt(IV_[e], flux_[e])
            ==
              fvm::laplacian(DAlbedo_, flux_[e])
            + fvm::Sp(F_[e], flux_[e]) / k_
            - fvm::Sp(M_[e], flux_[e])
            + DNP_->DNS()[e]
            + SNS_[e]
            + FNS_[e] / k_
        );

        fluxEqn.relax();
        fluxEqn.solve();

        flux_[e].correctBoundaryConditions();
    }

    updateSourceTerms();
}


// ************************************************************************* //