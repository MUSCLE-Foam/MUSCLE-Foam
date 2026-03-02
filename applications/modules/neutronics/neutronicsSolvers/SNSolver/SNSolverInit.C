#include "SNSolver.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::SNSolver::init()
{
    forAll(directions_, d)
    {
        directions_[d] = vector
        (
            directionsDict_.lookup("direction"+Foam::name(d))
        );

        directions_[d] /= Foam::mag(directions_[d]);
    }

    forAll(directionWeights_, d)
    {
        directionWeights_[d] = readScalar
        (
            directionsDict_.lookup("directionWeight"+Foam::name(d))
        );
    }

    forAll(flux_, e)
    {
        flux_[e].setSize(nDirections_);

        forAll(flux_[e], d)
        {
            flux_[e].set
            (
                d,
                new volScalarField
                (
                    IOobject
                    (
                        "flux"+Foam::name(e)+"_"+Foam::name(d),
                        mesh_.time().name(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_
                )
            );
        }
    }

    forAll(totalFlux_, e)
    {
        totalFlux_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "totalFlux"+Foam::name(e),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(inv(sqr(dimLength)*dimTime), Zero)
            )
        );
    }

    forAll(IV_, e)
    {
        IV_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "IV"+Foam::name(e),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );
    }

    forAll(chiPrompt_, e)
    {
        chiPrompt_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "chiPrompt"+Foam::name(e),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );
    }

    forAll(TOT_, e)
    {
        TOT_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "TOT"+Foam::name(e),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(inv(dimLength), Zero)
            )
        );

        TOT0_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "TOT"+Foam::name(e)+"_0",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        alphaTOT_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "alphaTOT"+Foam::name(e),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(inv(dimLength), Zero)
            )
        );
    }

    forAll(TR_, e)
    {
        TR_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "TR"+Foam::name(e),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(inv(dimLength), Zero)
            )
        );

        TR0_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "TR"+Foam::name(e)+"_0",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        alphaTR_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "alphaTR"+Foam::name(e),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(inv(dimLength), Zero)
            )
        );
    }

    forAll(M_, e)
    {
        M_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "M"+Foam::name(e),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(inv(dimLength), Zero)
            )
        );

        M0_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "M"+Foam::name(e)+"_0",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        alphaM_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "alphaM"+Foam::name(e),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(inv(dimLength), Zero)
            )
        );
    }

    forAll(SP_, e)
    {
        SP_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "SP"+Foam::name(e),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimForce, Zero)
            )
        );

        SP0_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "SP"+Foam::name(e)+"_0",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        alphaSP_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "alphaSP"+Foam::name(e),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimForce, Zero)
            )
        );
    }

    forAll(FN_, e)
    {
        FN_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "FN"+Foam::name(e),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(inv(dimLength), Zero)
            )
        );

        FN0_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "FN"+Foam::name(e)+"_0",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        alphaFN_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "alphaFN"+Foam::name(e),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(inv(dimLength), Zero)
            )
        );
    }

    forAll(S_, ord)
    {
        S_[ord].setSize(energyGroups_);

        forAll(S_[ord], ei)
        {
            S_[ord][ei].setSize(energyGroups_);

            forAll(S_[ord][ei], ej)
            {
                S_[ord][ei].set
                (
                    ej,
                    new volScalarField
                    (
                        IOobject
                        (
                            "S"+Foam::name(ord)+Foam::name(ei)+Foam::name(ej),
                            mesh_.time().name(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh_,
                        dimensionedScalar(inv(dimLength), Zero)
                    )
                );
            }
        }

        S0_[ord].setSize(energyGroups_);

        forAll(S0_[ord], ei)
        {
            S0_[ord][ei].setSize(energyGroups_);

            forAll(S0_[ord][ei], ej)
            {
                S0_[ord][ei].set
                (
                    ej,
                    new volScalarField
                    (
                        IOobject
                        (
                            "S"+Foam::name(ord)+Foam::name(ei)
                               +Foam::name(ej)+"_0",
                            mesh_.time().name(),
                            mesh_,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh_
                    )
                );
            }
        }

        alphaS_[ord].setSize(energyGroups_);

        forAll(alphaS_[ord], ei)
        {
            alphaS_[ord][ei].setSize(energyGroups_);

            forAll(alphaS_[ord][ei], ej)
            {
                alphaS_[ord][ei].set
                (
                    ej,
                    new volScalarField
                    (
                        IOobject
                        (
                            "alphaS"+Foam::name(ord)+Foam::name(ei)
                                    +Foam::name(ej),
                            mesh_.time().name(),
                            mesh_,
                            IOobject::READ_IF_PRESENT,
                            IOobject::NO_WRITE
                        ),
                        mesh_,
                        dimensionedScalar(inv(dimLength), Zero)
                    )
                );
            }
        }
    }

    forAll(F_, e)
    {
        F_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "F"+Foam::name(e),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(inv(dimLength), Zero)
            )
        );
    }

    forAll(neutronSource_, e)
    {
        neutronSource_[e].setSize(nDirections_);

        forAll(neutronSource_[e], d)
        {
            neutronSource_[e].set
            (
                d,
                new volScalarField
                (
                    IOobject
                    (
                        "neutronSource"+Foam::name(e)+"_"+Foam::name(d),
                        mesh_.time().name(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar(inv(pow3(dimLength)*dimTime), Zero)
                )
            );
        }
    }

    forAll(polyDir_, ord)
    {
        polyDir_[ord].setSize(nDirections_);

        forAll(polyDir_[ord], di)
        {
            polyDir_[ord][di].setSize(nDirections_);

            forAll(polyDir_[ord][di], dj)
            {
                scalar prod = directions_[di] & directions_[dj];

                scalar orderWeight = (2.0 * ord + 1.0) / 8.0;

                // Legendre Polynomial of prod
                scalar P = 1.0;

                if (ord == 1)
                {
                    P = prod;
                }
                else if (ord > 1)
                {
                    scalar P00 = 1.0;
                    scalar P0  = prod;

                    for (label i = 2; i <= ord; i++)
                    {
                        P = ((2.0*i - 1.0) * prod * P0 - (i - 1.0) * P00) / i;
                        P00 = P0;
                        P0 = P;
                    }
                }

                polyDir_[ord][di][dj] =
                    P * directionWeights_[dj] * orderWeight;
            }
        }
    }
}


// ************************************************************************* //