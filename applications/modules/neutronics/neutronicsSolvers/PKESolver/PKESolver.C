#include "PKESolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PKESolver, 0);
    addToRunTimeSelectionTable(neutronicsSolver, PKESolver, fvMesh);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::PKESolver::normalize(scalar normalization)
{
    this->neutronicsSolver::normalize(normalization);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PKESolver::PKESolver(fvMesh& mesh)
:
    neutronicsSolver(mesh),

    kineticsConstants_
    (
        IOobject
        (
            "kineticsConstants",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    L_("L", dimTime, readScalar(kineticsConstants_.lookup("L"))),
    reactivity_(readScalar(kineticsConstants_.lookup("reactivity"))),
    reactivityPlus_(0),
    psi_(readScalar(kineticsConstants_.lookup("psi"))),

    eta_(this->DNP_->DNPGroups, kineticsConstants_.lookupOrDefault("eta", 1.0)),
    movingPrecursors_
    (
        kineticsConstants_.lookupOrDefault("movingPrecursors", false)
    ),
    controlSystem_
    (
        kineticsConstants_.lookupOrDefault("controlSystem", false)
    ),

    tauc_
    (
        "tauc",
        dimTime,
        movingPrecursors_ ?
        readScalar(kineticsConstants_.lookup("tauc"))
        :
        -1
    ),
    taul_
    (
        "taul",
        dimTime,
        movingPrecursors_ ?
        readScalar(kineticsConstants_.lookup("taul"))
        :
        -1
    ),
    alphaT_
    (
        "alphaT",
        inv(dimTemperature),
        readScalar(kineticsConstants_.lookup("alphaT"))
    ),
    alphaRho_
    (
        "alphaRho",
        inv(dimDensity),
        readScalar(kineticsConstants_.lookup("alphaRho"))
    ),
    coreCells_(mesh.cellZones()["core"])
{
    // Compute core volume
    coreVolume_ = 0;

    const scalarField& V = mesh.V();

    forAll(coreCells_, i)
    {
        coreVolume_ += V[coreCells_[i]];
    }

    reduce(coreVolume_, sumOp<scalar>());

    // Moving precursor initialization

    scalar betaTotal = this->DNP_->betaTotal();

    if (movingPrecursors_)
    {
        reactivityPlus_ += betaTotal;

        forAll(eta_,i)
        {
            scalar lam = this->DNP_->lambda[i].value();
            scalar beta = this->DNP_->beta[i];
            reactivityPlus_ -= (beta*lam)/(lam+(1.0-Foam::exp(-lam*taul_.value()))/tauc_.value());

            eta_[i] = lam/(lam+(1.0-Foam::exp(-lam*taul_.value()))/tauc_.value());

            Info << "Initial eta" << i << ": " << eta_[i] << nl;
        }
    }

    // Initial temperature and density
    coreInitTemp_ = coreAverage(this->fluid.T());
    coreInitRho_ = coreAverage(this->fluid.rho());

    Info << "Initial core temperature: " << coreInitTemp_
        << " K. Initial density: " << coreInitRho_ << " kg/m3." << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PKESolver::~PKESolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PKESolver::solve()
{
    scalar deltaT = this->runTime.deltaTValue();
    scalar betaTotal = this->DNP_->betaTotal();
    scalar L = L_.value();

    scalar reactivity = reactivity_;

    if (controlSystem_)
    {
       reactivityPlus_ = betaTotal;

        forAll(eta_,i)
        {
            scalar lam = this->DNP_->lambda[i].value();
            scalar beta = this->DNP_->beta[i];
            reactivityPlus_ -= (beta*lam)/(lam+(1.0-Foam::exp(-lam*taul_.value()))/tauc_.value());
        }
    }

    if (movingPrecursors_)
    {
        reactivity += reactivityPlus_;
    }

    scalar coreTemp = coreAverage(this->fluid.T());
    scalar coreRho  = coreAverage(this->fluid.rho());

    scalar feedbackReactivity =
        alphaT_.value() * (coreTemp - coreInitTemp_)
      + alphaRho_.value() * (coreRho - coreInitRho_);

    reactivity += feedbackReactivity;

    Info << "Core temperature: " << coreTemp << " K" << nl
        << "Core density: " << coreRho << " kg/m3" << nl
        << "Feedback reactivity: " << feedbackReactivity << endl;


    // Implicitly compute flux at time n+1

    scalar sumNum = 0.0;
    scalar sumDenom = 0.0;

    forAll(eta_,i)
    {
        scalar lambda = this->DNP_->lambda[i].value();
        scalar beta = this->DNP_->beta[i];

        scalar denominator = 1.0 + lambda * deltaT;

        if (movingPrecursors_)
        {
            denominator += deltaT * (1.0-Foam::exp(-lambda*taul_.value()))/tauc_.value();
        }

        sumNum += (beta * eta_[i]) / denominator;
        sumDenom += (beta * lambda) / denominator;
    }

    // scalar reactivity2 = betaTotal - sumNum/psi_ - deltaT*sumDenom;

    // Info << "Reactivity: " << reactivity << endl;
    // Info << "Reactivity2: " << reactivity2 << endl;


    psi_ = (psi_ + (deltaT/L)*sumNum)
         / (1.0 - deltaT*(reactivity-betaTotal)/L - (Foam::sqr(deltaT)/L)*sumDenom);


    // Prevent non-physical negative values during aggressive transients
    psi_ = max(psi_, 0.0);

    // Advance precursor concentrations in time
    forAll(eta_,i)
    {
        scalar lam = this->DNP_->lambda[i].value();
        // scalar beta = this->DNP_->beta[i];

        if (movingPrecursors_)
        {
            eta_[i] = (lam * psi_ + eta_[i] / deltaT)
                    / (
                          1.0/deltaT + lam + 1.0/tauc_.value()
                        - Foam::exp(-lam*taul_.value())/tauc_.value()
                      );
        }
        else
        {
            eta_[i] = (lam * psi_ + eta_[i] / deltaT)
                    / (1.0/deltaT + lam);
        }
    }

    // Output flux and precursor concentrations
    Info << "Psi: " << psi_ << endl;

    forAll(eta_,i)
    {
        Info << "Eta" << i << ": " << eta_[i] << nl;
    }
    Info << endl;
}

void Foam::PKESolver::updatePower()
{
    scalar totalPower = this->nominalPower_ * psi_;

    this->powerDensity_ = Zero;

    for (label cellI : coreCells_)
    {
        this->powerDensity_[cellI] = totalPower / coreVolume_;
    }

    powerDensity_.correctBoundaryConditions();

    Info << "Total power: " << totalPower << " W" << endl;
}

Foam::scalar Foam::PKESolver::coreAverage(const volScalarField& field)
{
    const scalarField& V = this->mesh_.V();

    scalar integral = 0;

    forAll(coreCells_, i)
    {
        integral += V[coreCells_[i]] * field[coreCells_[i]];
    }

    reduce(integral, sumOp<scalar>());

    return integral/coreVolume_;
}

// ************************************************************************* //