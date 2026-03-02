#include "SP3Solver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SP3Solver, 0);
    addToRunTimeSelectionTable(neutronicsSolver, SP3Solver, fvMesh);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::SP3Solver::updateCrossSections()
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
            S2_[ei][ej] = (S20_[ei][ej] + alphaS2_[ei][ej] * logT) * K;
            S3_[ei][ej] = (S30_[ei][ej] + alphaS3_[ei][ej] * logT) * K;
        }
    }
}


void Foam::SP3Solver::updateSourceTerms()
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
        neutronFissionSource_ += FN_[ei] * (flux_[ei] - 2.0 * secondFlux_[ei]);
        powerSource_ += SP_[ei] * (flux_[ei] - 2.0 * secondFlux_[ei]);

        forAll(SNS_, ej)
        {
            if (ei != ej)
            {
                SNS_[ei] += S_[ej][ei] * (flux_[ej] - 2.0 * secondFlux_[ej]);
                FNS_[ei] += FN_[ej] * (1.0 - DNP_->betaTotal())
                    * (flux_[ej] - 2.0 * secondFlux_[ej]) * chiPrompt_[ei];
            }
        }
    }

    promptHeatSource_ = (1.0 - DNP_->betaTotal()) * powerSource_;

    DNP_->updateSourceTerms();
}


void Foam::SP3Solver::normalize(scalar normalization)
{
    forAll(flux_, e)
    {
        flux_[e] /= normalization;
    }

    forAll(secondFlux_, e)
    {
        secondFlux_[e] /= normalization;
    }

    this->neutronicsSolver::normalize(normalization);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SP3Solver::SP3Solver(fvMesh& mesh)
:
    neutronicsSolver(mesh),

    flux_(energyGroups_),
    secondFlux_(energyGroups_),

    fluxAlbedo_
    (
        IOobject
        (
            "fluxAlbedo",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
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

    S2_(energyGroups_),
    S20_(energyGroups_),
    alphaS2_(energyGroups_),

    S3_(energyGroups_),
    S30_(energyGroups_),
    alphaS3_(energyGroups_),

    F_(energyGroups_),

    DAlbedo_
    (
        IOobject
        (
            "DAlbedo",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
{
    init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SP3Solver::~SP3Solver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::SP3Solver::solve()
{
    updateCrossSections();

    forAll(flux_, e)
    {
        DAlbedo_ = 1.0 / (3.0 * TR_[e]);

        fluxAlbedo_ = secondFlux_[e];

        fvScalarMatrix fluxEqn
        (
              fvm::ddt(IV_[e], flux_[e])
            ==
              fvm::laplacian(DAlbedo_, flux_[e])
            + fvm::Sp(F_[e], flux_[e])                / k_
            - fvc::Sp(F_[e], secondFlux_[e])    * 2.0 / k_
            - fvm::Sp(M_[e], flux_[e])
            + fvc::Sp(M_[e], secondFlux_[e])    * 2.0
            + DNP_->DNS()[e]
            + SNS_[e]
            + FNS_[e]                                 / k_
            + fvc::ddt(IV_[e], secondFlux_[e])  * 2.0

        );

        fluxEqn.relax();
        fluxEqn.solve();

        flux_[e].correctBoundaryConditions();
    }

    forAll(secondFlux_, e)
    {
        DAlbedo_ = 9.0 / (35.0 * TR_[e]);

        fluxAlbedo_ = flux_[e];

        fvScalarMatrix secondFluxEqn
        (
               9.0/5.0 * fvm::ddt(IV_[e], secondFlux_[e])
            ==
              fvm::laplacian(DAlbedo_, secondFlux_[e])
            - fvm::Sp(TR_[e], secondFlux_[e])
            - 4.0/5.0 * fvm::Sp(M_[e], secondFlux_[e])
            + 2.0/5.0 * fvc::Sp(M_[e], flux_[e])
            + 4.0/5.0 * fvm::Sp(F_[e], secondFlux_[e]) / k_
            - 2.0/5.0 * fvc::Sp(F_[e], flux_[e]) / k_
            - 2.0/5.0 * DNP_->DNS()[e]
            - 2.0/5.0 * SNS_[e]
            - 2.0/5.0 * FNS_[e] / k_
            + 2.0/5.0 * fvc::ddt(IV_[e], flux_[e])
        );

        secondFluxEqn.relax();
        secondFluxEqn.solve();

        secondFlux_[e].correctBoundaryConditions();
    }

    updateSourceTerms();
}


// ************************************************************************* //