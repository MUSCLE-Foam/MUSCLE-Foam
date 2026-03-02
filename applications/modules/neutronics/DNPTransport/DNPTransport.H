/*---------------------------------------------------------------------------*\

    MUSCLE-Foam: an OpenFOAM-based multiphysics solver for reactor physics
    Copyright (C) 2026 NRG PALLAS

\-----------------------------------------------------------------------------/

License
    This file is part of MUSCLE-Foam. It is published under the terms of the
    GNU General Public License. See <http://www.gnu.org/licenses/>.

Class
    Foam::DNPTransport

Description
    Transport module for Delayed Neutron Precursor (DNP).

    A user-specified number of DNP groups are transported by diffusion as well
    as fuel/coolant motion (if any).

    Reference:
    \verbatim
        Cervi, E., Lorenzi, S., Cammi, A. & Luzzi, L. (2019)
        Development of a multiphysics model for the study of fuel
        compressibility effects in the Molten Salt Fast Reactor.
        Chemical Engineering Science, 193, 379-393.
    \endverbatim

SourceFiles
    DNPTransport.C

\*---------------------------------------------------------------------------*/

#ifndef DNPTransport_H
#define DNPTransport_H

#include "fvMesh.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class DNPTransport Declaration
\*---------------------------------------------------------------------------*/

class neutronicsSolver;

class DNPTransport
{
protected:

    // Protected Data

        //- Region mesh
        const fvMesh& mesh_;

private:

    // Private Member Data

        //- Pointer to neutronics solver
        neutronicsSolver& neutronicsSolver_;

        //- Number of delayed neutron precursor groups
        label DNPGroups_;

        //- Turbulent Schmidt number
        scalar Sct_;

        //- Molecular Schmidt number,
        //  defaults to a large number such that molecular diffusion
        //  is neglected when the molecular Schmidt number is not known
        scalar Sc_;

        //- Delayed neutron precursor densities
        PtrList<volScalarField> prec_;

        //- Delayed neutron precursor constants
        List<dimensionedScalar> lambda_;

        //- Delayed neutron precursor fractions
        scalarList beta_;

        //- Total delayed neutron fraction
        scalar betaTotal_;

        //- Delayed neutron sources per energy group
        PtrList<volScalarField> DNS_;

        //- Delayed neutron yields per energy group
        PtrList<volScalarField> chi_;


    // Private Member Functions

        //- Read delayed neutron constants dictionary
        void readDelayedNeutronConstants(IOdictionary&);

public:

    // Public Data


    //- Runtime type information
    TypeName("DNPTransport");


    // Constructors

        //- Construct from region mesh
        DNPTransport(neutronicsSolver&, IOdictionary&);

        //- Disallow default bitwise copy construction
        DNPTransport(const DNPTransport&) = delete;

    //- Destructor
    virtual ~DNPTransport();


    // Member Functions

        //- Solve the delayed neutron precursor transport equations
        void solve();

        //- Normalize delayed neutron source
        void normalize(scalar);

        //- Update delayed neutron source
        void updateSourceTerms();

        //- Return reference to delayed neutron sources
        inline const PtrList<volScalarField>& DNS() { return DNS_; }

        //- Return total delayed neutron fraction
        inline scalar betaTotal() { return betaTotal_; }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const DNPTransport&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //