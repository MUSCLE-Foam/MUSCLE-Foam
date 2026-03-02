#include "SP3Solver.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::SP3Solver::init()
{
    forAll(flux_, energyGroup)
    {
        flux_.set
        (
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "flux"+Foam::name(energyGroup),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }

    forAll(secondFlux_, energyGroup)
    {
        secondFlux_.set
        (
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "secondFlux"+Foam::name(energyGroup),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }

    forAll(IV_, energyGroup)
    {
        IV_.set
        (
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "IV"+Foam::name(energyGroup),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );
    }

    forAll(chiPrompt_, energyGroup)
    {
        chiPrompt_.set
        (
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "chiPrompt"+Foam::name(energyGroup),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );
    }

    forAll(TR_, energyGroup)
    {
        TR_.set
        (
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "TR"+Foam::name(energyGroup),
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
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "TR"+Foam::name(energyGroup)+"_0",
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
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "alphaTR"+Foam::name(energyGroup),
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

    forAll(M_, energyGroup)
    {
        M_.set
        (
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "M"+Foam::name(energyGroup),
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
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "M"+Foam::name(energyGroup)+"_0",
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
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "alphaM"+Foam::name(energyGroup),
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

    forAll(SP_, energyGroup)
    {
        SP_.set
        (
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "SP"+Foam::name(energyGroup),
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
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "SP"+Foam::name(energyGroup)+"_0",
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
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "alphaSP"+Foam::name(energyGroup),
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

    forAll(FN_, energyGroup)
    {
        FN_.set
        (
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "FN"+Foam::name(energyGroup),
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
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "FN"+Foam::name(energyGroup)+"_0",
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
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "alphaFN"+Foam::name(energyGroup),
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

    forAll(S_, ei)
    {
        S_[ei].setSize(energyGroups_);

        forAll(S_[ei], ej)
        {
            S_[ei].set
            (
                ej,
                new volScalarField
                (
                    IOobject
                    (
                        "S"+Foam::name(ei)+Foam::name(ej),
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

        S0_[ei].setSize(energyGroups_);

        forAll(S0_[ei], ej)
        {
            S0_[ei].set
            (
                ej,
                new volScalarField
                (
                    IOobject
                    (
                        "S"+Foam::name(ei)+Foam::name(ej)+"_0",
                        mesh_.time().name(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_
                )
            );
        }

        alphaS_[ei].setSize(energyGroups_);

        forAll(alphaS_[ei], ej)
        {
            alphaS_[ei].set
            (
                ej,
                new volScalarField
                (
                    IOobject
                    (
                        "alphaS"+Foam::name(ei)+Foam::name(ej),
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

    forAll(S2_, ei)
    {
        S2_[ei].setSize(energyGroups_);

        forAll(S2_[ei], ej)
        {
            S2_[ei].set
            (
                ej,
                new volScalarField
                (
                    IOobject
                    (
                        "S2"+Foam::name(ei)+Foam::name(ej),
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

        S20_[ei].setSize(energyGroups_);

        forAll(S20_[ei], ej)
        {
            S20_[ei].set
            (
                ej,
                new volScalarField
                (
                    IOobject
                    (
                        "S2"+Foam::name(ei)+Foam::name(ej)+"_0",
                        mesh_.time().name(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_
                )
            );
        }

        alphaS2_[ei].setSize(energyGroups_);

        forAll(alphaS2_[ei], ej)
        {
            alphaS2_[ei].set
            (
                ej,
                new volScalarField
                (
                    IOobject
                    (
                        "alphaS2"+Foam::name(ei)+Foam::name(ej),
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

    forAll(S3_, ei)
    {
        S3_[ei].setSize(energyGroups_);

        forAll(S3_[ei], ej)
        {
            S3_[ei].set
            (
                ej,
                new volScalarField
                (
                    IOobject
                    (
                        "S3"+Foam::name(ei)+Foam::name(ej),
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

        S30_[ei].setSize(energyGroups_);

        forAll(S30_[ei], ej)
        {
            S30_[ei].set
            (
                ej,
                new volScalarField
                (
                    IOobject
                    (
                        "S3"+Foam::name(ei)+Foam::name(ej)+"_0",
                        mesh_.time().name(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_
                )
            );
        }

        alphaS3_[ei].setSize(energyGroups_);

        forAll(alphaS3_[ei], ej)
        {
            alphaS3_[ei].set
            (
                ej,
                new volScalarField
                (
                    IOobject
                    (
                        "alphaS3"+Foam::name(ei)+Foam::name(ej),
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

    forAll(F_, energyGroup)
    {
        F_.set
        (
            energyGroup,
            new volScalarField
            (
                IOobject
                (
                    "F"+Foam::name(energyGroup),
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


// ************************************************************************* //