#include "SNSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SNSolver, 0);
    addToRunTimeSelectionTable(neutronicsSolver, SNSolver, fvMesh);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::SNSolver::updateCrossSections()
{
    volScalarField K(fluid.rho()/rhoXS_);
    volScalarField logT(log(fluid.T()/TXS_));

    forAll(flux_, ei)
    {
        TOT_[ei] = (TOT0_[ei] + alphaTOT_[ei] * logT) * K;
        TR_[ei] = (TR0_[ei] + alphaTR_[ei] * logT) * K;
        M_[ei]  = (M0_[ei]  + alphaM_[ei]  * logT) * K;
        SP_[ei] = (SP0_[ei] + alphaSP_[ei] * logT) * K;
        FN_[ei] = (FN0_[ei] + alphaFN_[ei] * logT) * K;

        F_[ei] = FN_[ei] * (1.0 - DNP_->betaTotal()) * chiPrompt_[ei];
    }

    forAll(S_, ord)
    {
        forAll(S_[ord], ei)
        {
            forAll(S_[ord][ei], ej)
            {
                S_[ord][ei][ej] =
                    (S0_[ord][ei][ej] + alphaS_[ord][ei][ej] * logT) * K;
            }
        }
    }
}


void Foam::SNSolver::updateSourceTerms()
{
    neutronFissionSource_ = Zero;
    powerSource_ = Zero;

    forAll(totalFlux_, ei)
    {
        neutronFissionSource_ += FN_[ei] * totalFlux_[ei];
        powerSource_ += SP_[ei] * totalFlux_[ei];
    }

    promptHeatSource_ = (1.0 - DNP_->betaTotal()) * powerSource_;

    DNP_->updateSourceTerms();
}


void Foam::SNSolver::normalize(scalar normalization)
{
    forAll(flux_, e)
    {
        totalFlux_[e] /= normalization;

        forAll(flux_[e], d)
        {
            flux_[e][d] /= normalization;
        }
    }

    this->neutronicsSolver::normalize(normalization);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SNSolver::SNSolver(fvMesh& mesh)
:
    neutronicsSolver(mesh),

    directionsDict_
    (
        IOobject
        (
            "neutronFlightDirectionsDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    nDirections_(directionsDict_.lookup<label>("nDirections")),
    directions_(nDirections_),
    directionWeights_(nDirections_),

    flux_(energyGroups_),
    totalFlux_(energyGroups_),

    scatteringOrders_
    (
        mesh.time().controlDict().lookup<label>("scatteringOrders")
    ),

    IV_(energyGroups_),
    chiPrompt_(energyGroups_),

    TOT_(energyGroups_),
    TOT0_(energyGroups_),
    alphaTOT_(energyGroups_),

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

    S_(scatteringOrders_),
    S0_(scatteringOrders_),
    alphaS_(scatteringOrders_),

    F_(energyGroups_),

    neutronSource_(energyGroups_),

    phif_
    (
        IOobject
        (
            "phif",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(sqr(dimLength), Zero)
    ),

    polyDir_(scatteringOrders_)
{
    init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SNSolver::~SNSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::SNSolver::solve()
{
    updateCrossSections();

    forAll(flux_, ei)
    {
        forAll(flux_[ei], di)
        {
            phif_ = directions_[di] & mesh.Sf();

            neutronSource_[ei][di] =
                (1.0/8.0) * DNP_->DNS()[ei] * directionWeights_[di];

            forAll(flux_, ej)
            {
                forAll(flux_[ej], dj)
                {
                    if (!(ej == ei && dj == di))
                    {
                        neutronSource_[ei][di] += (1.0/8.0) / k_ * FN_[ej]
                            * chiPrompt_[ei] * (1.0 - DNP_->betaTotal())
                            * flux_[ej][dj] * directionWeights_[dj];

                        forAll(polyDir_, ord)
                        {
                            neutronSource_[ei][di] += polyDir_[ord][di][dj]
                                * S_[ord][ej][ei] * flux_[ej][dj];
                        }
                    }
                }
            }

            fvScalarMatrix fluxEqn
            (
                  fvm::ddt(IV_[ei], flux_[ei][di])
                + fvm::div(phif_, flux_[ei][di])
                + fvm::Sp(TOT_[ei], flux_[ei][di])
                ==
                  fvm::Sp
                  (
                      (1.0/8.0) / k_ * FN_[ei] * chiPrompt_[ei]
                    * (1.0 - DNP_->betaTotal()) * directionWeights_[di],
                      flux_[ei][di]
                  )
                + neutronSource_[ei][di]
            );

            forAll(S_, ord)
            {
                scalar orderWeight = (2.0 * ord + 1.0) / 8.0;

                fluxEqn -= fvm::Sp
                (
                    S_[ord][ei][ei] * orderWeight * directionWeights_[di],
                    flux_[ei][di]
                );
            }

            fluxEqn.relax();
            fluxEqn.solve();

            flux_[ei][di].correctBoundaryConditions();
        }
    }

    forAll(totalFlux_, e)
    {
        totalFlux_[e] = Zero;

        forAll(directionWeights_, d)
        {
            totalFlux_[e] += (1.0/8.0) * directionWeights_[d] * flux_[e][d];
        }
    }

    updateSourceTerms();
}


// ************************************************************************* //