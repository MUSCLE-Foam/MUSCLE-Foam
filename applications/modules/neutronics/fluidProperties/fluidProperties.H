/*---------------------------------------------------------------------------*\

    MUSCLE-Foam: an OpenFOAM-based multiphysics solver for reactor physics
    Copyright (C) 2026 NRG PALLAS

\-----------------------------------------------------------------------------/

License
    This file is part of MUSCLE-Foam. It is published under the terms of the
    GNU General Public License. See <http://www.gnu.org/licenses/>.

Class
    Foam::fluidProperties

Description
    Fluid properties for neutronics solver. For single-phase solvers, return
    single phase properties. For multiphase solvers, calculate and return
    multiphase mixture properties.

SourceFiles
    fluidProperties.C

\*---------------------------------------------------------------------------*/

#ifndef fluidProperties_H
#define fluidProperties_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fluidProperties Declaration
\*---------------------------------------------------------------------------*/

class fluidProperties
{
private:

    // Private Member Data

        //- Region mesh
        const fvMesh& mesh_;

        //- Is the fluid a mixture?
        bool mixture_;

        //- List of phases, empty if single-phase solver
        wordList phaseNames_;

    // Private member functions

        //- Compute mixture properties given a property name
        tmp<volScalarField> mixtureProperty(word, dimensionSet) const;

public:

    //- Runtime type information
    TypeName("fluidProperties");


    // Constructors

        //- Construct from region mesh
        fluidProperties(const fvMesh&);

        //- Disallow default bitwise copy construction
        fluidProperties(const fluidProperties&) = delete;

    //- Destructor
    virtual ~fluidProperties();


    // Member Functions

        //- Return the (mixture) density field
        tmp<volScalarField> rho() const;

        //- Return the (mixture) temperature field
        tmp<volScalarField> T() const;

        //- Return the (mixture) viscosity field
        tmp<volScalarField> mu() const;

        //- Return the (mixture) turbulent viscosity field, if any
        tmp<volScalarField> nut() const;

        //- Return the (mixture) volumetric flux
        tmp<surfaceScalarField> phi() const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fluidProperties&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //