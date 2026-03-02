#include "fluidProperties.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidProperties, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fluidProperties::mixtureProperty
(
    word propName,
    dimensionSet dims
) const
{
    tmp<volScalarField> tMixtureProp
    (
        new volScalarField
        (
            IOobject
            (
                "mixtureProperty",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dims, Zero)
        )
    );

    volScalarField& mixtureProp = tMixtureProp.ref();

    forAll(phaseNames_, i)
    {
        const word phaseProp = propName + "." + phaseNames_[i];

        if (mesh_.foundObject<volScalarField>(phaseProp))
        {
            const volScalarField& alpha =
                mesh_.lookupObject<volScalarField>("alpha." + phaseNames_[i]);

            const volScalarField& prop =
                mesh_.lookupObject<volScalarField>(phaseProp);

            mixtureProp += alpha * prop;
        }
    }

    return tMixtureProp;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidProperties::fluidProperties(const fvMesh& mesh)
:
    mesh_(mesh),

    mixture_(!mesh_.foundObject<volScalarField>("rho")),

    phaseNames_()
{
    if (mixture_)
    {
        const word alphaPrefix = "alpha.";

        forAll(mesh_.objectRegistry::toc<volScalarField>(), i)
        {
            const word fieldName =
                mesh_.objectRegistry::toc<volScalarField>()[i];

            if
            (
                fieldName.find(alphaPrefix) != std::string::npos
                &&
                !
                (
                    fieldName.size() > 2
                    && fieldName.rfind("_0") == fieldName.size() - 2
                )
            )
            {
                phaseNames_.append(fieldName.substr(alphaPrefix.size()));
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidProperties::~fluidProperties()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fluidProperties::rho() const
{
    if (mesh_.foundObject<volScalarField>("rho"))
    {
        return tmp<volScalarField>
        (
            mesh_.lookupObject<volScalarField>("rho")
        );
    }
    else
    {
        return mixtureProperty("rho", dimDensity);
    }
}

Foam::tmp<Foam::volScalarField> Foam::fluidProperties::T() const
{
    if (mesh_.foundObject<volScalarField>("T"))
    {
        return tmp<volScalarField>
        (
            mesh_.lookupObject<volScalarField>("T")
        );
    }
    else
    {
        return mixtureProperty("T", dimTemperature);
    }
}

Foam::tmp<Foam::volScalarField> Foam::fluidProperties::mu() const
{
    if (mesh_.foundObject<volScalarField>("mu"))
    {
        return tmp<volScalarField>
        (
            mesh_.lookupObject<volScalarField>("mu")
        );
    }
    else
    {
        return mixtureProperty("mu", dimDynamicViscosity);
    }
}

Foam::tmp<Foam::volScalarField> Foam::fluidProperties::nut() const
{
    if (mesh_.foundObject<volScalarField>("nut"))
    {
        return tmp<volScalarField>
        (
            mesh_.lookupObject<volScalarField>("nut")
        );
    }
    else if
    (
        mixture_
        && mesh_.toc<volScalarField>(wordRe("nut.*")).size()
    )
    {
        return mixtureProperty("nut", dimKinematicViscosity);
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "nut",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimKinematicViscosity, Zero)
            )
        );
    }
}

Foam::tmp<Foam::surfaceScalarField> Foam::fluidProperties::phi() const
{
    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>("phi");

    if (phi.dimensions() == dimVolume / dimTime)
    {
        return tmp<surfaceScalarField>(phi);
    }
    else
    {
        tmp<surfaceScalarField> tPhi
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "tempPhi",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimVolume / dimTime, Zero)
            )
        );

        surfaceScalarField& vPhi = tPhi.ref();

        vPhi = phi / fvc::interpolate(rho());

        return tPhi;
    }
}

// ************************************************************************* //