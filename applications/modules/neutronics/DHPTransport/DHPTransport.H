/*---------------------------------------------------------------------------*\

    MUSCLE-Foam: an OpenFOAM-based multiphysics solver for reactor physics
    Copyright (C) 2026 NRG PALLAS

\-----------------------------------------------------------------------------/

License
    This file is part of MUSCLE-Foam. It is published under the terms of the
    GNU General Public License. See <http://www.gnu.org/licenses/>.

Class
    Foam::DHPTransport

Description
    Transport module for Decay Heat Precursor (DHP).

    A user-specified number of DHP groups are transported by diffusion as well
    as fuel/coolant motion (if any).

    Reference:
    \verbatim
        Cervi, E., Lorenzi, S., Cammi, A. & Luzzi, L. (2019)
        Development of a multiphysics model for the study of fuel
        compressibility effects in the Molten Salt Fast Reactor.
        Chemical Engineering Science, 193, 379-393.
    \endverbatim

SourceFiles
    DHPTransport.C

\*---------------------------------------------------------------------------*/

#ifndef DHPTransport_H
#define DHPTransport_H

#include "fvMesh.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class DHPTransport Declaration
\*---------------------------------------------------------------------------*/

class neutronicsSolver;

class DHPTransport
{
protected:

    // Protected Data

        //- Region mesh
        const fvMesh& mesh_;

private:

    // Private Member Data

        //- Pointer to neutronics solver
        neutronicsSolver& neutronicsSolver_;

        //- Number of decay heat precursor groups
        label DHPGroups_;

        //- Turbulent Schmidt number
        scalar Sct_;

        //- Molecular Schmidt number,
        //  defaults to a large number such that molecular diffusion
        //  is neglected when the molecular Schmidt number is not known
        scalar Sc_;

        //- Decay heat source
        volScalarField decayHeatSource_;

        //- Decay heat precursor density
        PtrList<volScalarField> dec_;

        //- Decay heat precursor constants
        List<dimensionedScalar> lambda_;

        //- Decay heat precursor fractions
        scalarList beta_;


    // Private Member Functions

        //- Read delayed neutron constants dictionary
        void readDecayHeatConstants(IOdictionary&);


public:

    // Public Data


    //- Runtime type information
    TypeName("DHPTransport");


    // Constructors

        //- Construct from region mesh and dictionary
        DHPTransport(neutronicsSolver&, IOdictionary&);

        //- Disallow default bitwise copy construction
        DHPTransport(const DHPTransport&) = delete;

    //- Destructor
    virtual ~DHPTransport();


    // Member Functions

        //- Solve the decay heat precursor transport equations
        void solve();

        //- Normalize decay heat source
        void normalize(scalar);

        //- Update decay heat source
        void updateSourceTerms();


        //- Return decay heat source
        inline volScalarField& decayHeatSource() { return decayHeatSource_; }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const DHPTransport&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //