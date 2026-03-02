/*---------------------------------------------------------------------------*\

    MUSCLE-Foam: an OpenFOAM-based multiphysics solver for reactor physics
    Copyright (C) 2026 NRG PALLAS

\-----------------------------------------------------------------------------/

License
    This file is part of MUSCLE-Foam. It is published under the terms of the
    GNU General Public License. See <http://www.gnu.org/licenses/>.

Class
    Foam::albedoFvPatchField

Description
    Albedo boundary condition for neutron fluxes and flux moments.

    For the diffusion solver, the albedo boundary condition for the flux is
    given by:

    \f[
        \nabla \varphi_i = -\frac{\gamma}{D_i}\varphi_i
    \f]

    where

    \f[
        \gamma = \frac{1}{2}\frac{(1-\alpha)}{(1+\alpha)}
    \f]

    and
    \vartable
        \varphi_i | ith group neutron flux [m^-2 s^-1]
        \alpha  | albedo coefficient [-]
        D_i     | diffusion coefficient [m]
    \endvartable

    For the SP3 solver, the albedo boundary condition for the flux is given by:
    \f[
        \nabla \varphi_i
            = -\frac{\gamma}{D_i} (\varphi_i-\frac{3}{4}\varphi_{2,i})
    \f]

    and for the second moment of the flux:

    \f[
        \nabla\varphi_{2,i} = \frac{\gamma}{D_{2,i}}
            (\frac{3}{20}\varphi_i-\frac{21}{20}\varphi_{2,i})
    \f]

    where
    \vartable
        \varphi_{2,i} | 2nd moment of ith group neutron flux [m^-2 s^-1]
        D_{2,i} | diffusion coefficient of 2nd moment of ith flux [m]
    \endvartable

    Reference:
    \verbatim
        Fiorina, C., Hursin, M., & Pautz, A. (2017)
        Extension of the GeN-Foam neutronic solver to SP3 analysis and
        application to the CROCUS experimental reactor
        Annals of Nuclear Energy, 101, 419-428.
    \endverbatim

SourceFiles
    albedoFvPatchField.C

    Acknowledgement:
        These source files draw inspiration from GeN-Foam's
        albedoSP3FvPatchField. The work of the GeN-Foam team is gratefully
        acknowledged, and no infringement of copyright is intended by the
        publication of this file.

\*---------------------------------------------------------------------------*/

#ifndef albedoFvPatchField_H
#define albedoFvPatchField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class albedoFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class albedoFvPatchField
:
    public fvPatchField<scalar>
{
private:

    // Private data

        //- Coefficient gamma
        scalar gamma_;

        //- Reference to diffusion coefficient patch field
        const Field<scalar>& DAlbedo_;

        //- Reference to other flux (for SP3)
        const Field<scalar>& fluxAlbedo_;

        //- Checks if the boundary condition is set for a second flux moment
        Switch forSecondMoment_;

public:

    //- Runtime type information
    TypeName("albedo");


    // Constructors

        //- Construct from patch and internal field
        albedoFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        albedoFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping the given albedoFvPatchField
        //  onto a new patch
        albedoFvPatchField
        (
            const albedoFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&,
            const bool mappingRequired=true
        );

        //- Disallow copy without setting internal field reference
        albedoFvPatchField(const albedoFvPatchField&) = delete;

        //- Copy constructor setting internal field reference
        albedoFvPatchField
        (
            const albedoFvPatchField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new albedoFvPatchField(*this, iF)
            );
        }


    // Member Functions

        // Attributes

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return true;
            }

            //- Return false: this patch field is not altered by assignment
            virtual bool assignable() const
            {
                return false;
            }


        // Evaluation functions

            //- Return patch-normal gradient
            virtual tmp<Field<scalar>> snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<scalar>> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<scalar>> gradientBoundaryCoeffs() const;


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        virtual void operator=(const UList<scalar>&) {}

        virtual void operator=(const fvPatchField<scalar>&) {}
        virtual void operator+=(const fvPatchField<scalar>&) {}
        virtual void operator-=(const fvPatchField<scalar>&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}

        virtual void operator+=(const Field<scalar>&) {}
        virtual void operator-=(const Field<scalar>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const scalar&) {}
        virtual void operator+=(const scalar&) {}
        virtual void operator-=(const scalar&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //