#include "albedoFvPatchField.H"
#include "GeometricFields.H"
#include "fieldTypes.H"
#include "IOobject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "fieldMapper.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::albedoFvPatchField::albedoFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar,volMesh>& iF
)
:
    fvPatchField<scalar>(p, iF),
    gamma_(0.5),
    DAlbedo_
    (
        this->patch().template lookupPatchField<volScalarField,scalar>
        (
            "DAlbedo"
        )
    ),
    fluxAlbedo_
    (
        this->patch().template lookupPatchField<volScalarField,scalar>
        (
            "fluxAlbedo"
        )
    ),
    forSecondMoment_(false)
{}


Foam::albedoFvPatchField::albedoFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar,volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fvPatchField<scalar>(p, iF, dict, valueRequired),
    gamma_(readScalar(dict.lookup("gamma"))),
    DAlbedo_
    (
        this->patch().template lookupPatchField<volScalarField,scalar>
        (
            "DAlbedo"
        )
    ),
    fluxAlbedo_
    (
        this->patch().template lookupPatchField<volScalarField,scalar>
        (
            "fluxAlbedo"
        )
    ),
    forSecondMoment_
    (
        this->internalField().name().find("secondFlux") != std::string::npos
    )
{}


Foam::albedoFvPatchField::albedoFvPatchField
(
    const albedoFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar,volMesh>& iF,
    const fieldMapper& mapper,
    const bool mappingRequired
)
:
    fvPatchField<scalar>(ptf, p, iF, mapper, mappingRequired),
    gamma_(ptf.gamma_),
    DAlbedo_(ptf.DAlbedo_),
    fluxAlbedo_(ptf.fluxAlbedo_),
    forSecondMoment_(ptf.forSecondMoment_)
{}


Foam::albedoFvPatchField::albedoFvPatchField
(
    const albedoFvPatchField& ptf,
    const DimensionedField<scalar,volMesh>& iF
)
:
    fvPatchField<scalar>(ptf, iF),
    gamma_(ptf.gamma_),
    DAlbedo_(ptf.DAlbedo_),
    fluxAlbedo_(ptf.fluxAlbedo_),
    forSecondMoment_(ptf.forSecondMoment_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar>> Foam::albedoFvPatchField::snGrad() const
{
    return
    (
        - this->patchInternalField() * gamma_ / DAlbedo_
    );
}


void Foam::albedoFvPatchField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    if (forSecondMoment_)
    {
        Field<scalar>::operator=
        (
            this->patchInternalField()
            + gamma_ / DAlbedo_
            *
            (
                21.0/20.0 * this->patchInternalField()
                - 3.0/20.0 * fluxAlbedo_
            )
            / this->patch().deltaCoeffs()
        );
    }
    else
    {
        Field<scalar>::operator=
        (
            this->patchInternalField()
            + gamma_ / DAlbedo_
            *
            (
                this->patchInternalField()
                - 3.0/4.0 * fluxAlbedo_
            )
            / this->patch().deltaCoeffs()
        );
    }

    fvPatchField<scalar>::evaluate();
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::albedoFvPatchField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    if (forSecondMoment_)
    {
        return
        (
            1.0 + gamma_ / DAlbedo_ * 21.0/20.0
            / this->patch().deltaCoeffs()
        );
    }
    else
    {
        return
        (
            1.0 + gamma_ / DAlbedo_
            / this->patch().deltaCoeffs()
        );
    }
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::albedoFvPatchField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    if (forSecondMoment_)
    {
        return
        (
            - fluxAlbedo_ * gamma_ / DAlbedo_ * 3.0/20.0
            / this->patch().deltaCoeffs()
        );
    }
    else
    {
        return
        (
            - fluxAlbedo_ * gamma_ / DAlbedo_ * 3.0/4.0
            / this->patch().deltaCoeffs()
        );
    }
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::albedoFvPatchField::gradientInternalCoeffs() const
{
    if (forSecondMoment_)
    {
        return (- gamma_ / DAlbedo_ * 21.0/20.0);
    }
    else
    {
        return (- gamma_ / DAlbedo_);
    }
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::albedoFvPatchField::gradientBoundaryCoeffs() const
{
    if (forSecondMoment_)
    {
        return (fluxAlbedo_ * gamma_ / DAlbedo_ * 3.0/20.0);
    }
    else
    {
        return (fluxAlbedo_ * gamma_ / DAlbedo_ * 3.0/4.0);
    }
}


void Foam::albedoFvPatchField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeKeyword("gamma")
        << gamma_ << token::END_STATEMENT << nl;
    os.writeKeyword("forSecondMoment")
        << forSecondMoment_ << token::END_STATEMENT << nl;

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField(fvPatchScalarField, albedoFvPatchField);
}


// ************************************************************************* //