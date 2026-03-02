/*---------------------------------------------------------------------------*\

    MUSCLE-Foam: an OpenFOAM-based multiphysics solver for reactor physics
    Copyright (C) 2026 NRG PALLAS

\-----------------------------------------------------------------------------/

License
    This file is part of MUSCLE-Foam. It is published under the terms of the
    GNU General Public License. See <http://www.gnu.org/licenses/>.

Class
    Foam::SNSolver

Description
    Discrete ordinate (S_N) solver for neutronics.

    The solver discretizes neutron fluxes based on energy groups as well as on
    the basis of the direction of flight of neutrons.

    Reference:
    \verbatim
        Cervi, E. & Cammi, A. (2020)
        A coupled neutronics shock physics solver: implementation of an SN
        neutron transport module.
        arXiv preprint - arXiv:2001.02291.
    \endverbatim

SourceFiles
    SNSolver.C
    SNSolverInit.C

\*---------------------------------------------------------------------------*/

#ifndef SNSolver_H
#define SNSolver_H

#include "neutronicsSolver.H"
#include "vectorList.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class SNSolver Declaration
\*---------------------------------------------------------------------------*/

class SNSolver
:
    public neutronicsSolver
{
private:

    // Private Member Data

        //- Flight direction dictionary
        IOdictionary directionsDict_;

        //- Number of neutron flight directions
        label nDirections_;

        //- Neutron flight directions. Must have magnitude of 1.
        vectorList directions_;

        //- Neutron flight direction weights
        scalarList directionWeights_;

        //- Neutron fluxes (energy, direction)
        List<PtrList<volScalarField>> flux_;

        //- Total neutron fluxes (energy)
        PtrList<volScalarField> totalFlux_;

        //- Number of scattering orders
        label scatteringOrders_;

        //- Inverse neutron velocities
        PtrList<volScalarField> IV_;

        //- Prompt neutron yields
        PtrList<volScalarField> chiPrompt_;

        //- Total cross-sections
        PtrList<volScalarField> TOT_;
        PtrList<volScalarField> TOT0_;
        PtrList<volScalarField> alphaTOT_;

        //- Transport cross sections
        PtrList<volScalarField> TR_;
        PtrList<volScalarField> TR0_;
        PtrList<volScalarField> alphaTR_;

        //- Combined fission, capture and scattering cross sections
        PtrList<volScalarField> M_;
        PtrList<volScalarField> M0_;
        PtrList<volScalarField> alphaM_;

        //- Fission cross-section times energy per fission
        PtrList<volScalarField> SP_;
        PtrList<volScalarField> SP0_;
        PtrList<volScalarField> alphaSP_;

        //- Fission cross-section times neutrons per fission
        PtrList<volScalarField> FN_;
        PtrList<volScalarField> FN0_;
        PtrList<volScalarField> alphaFN_;

        //- Scattering cross sections
        //  (scattering order, source group, target group)
        List<List<PtrList<volScalarField>>> S_;
        List<List<PtrList<volScalarField>>> S0_;
        List<List<PtrList<volScalarField>>> alphaS_;

        //- Normalized prompt neutron sources
        PtrList<volScalarField> F_;

        //- Neutron source
        List<PtrList<volScalarField>> neutronSource_;

        //- Surface flux
        surfaceScalarField phif_;

        //- Scattering polynomials
        //  (scattering order, source direction, target direction)
        List<List<scalarList>> polyDir_;

    // Private Member Functions

        //- Initialize neutronics fields
        void init();

        //- Update cross-sections
        void updateCrossSections();

        //- Update source terms
        void updateSourceTerms();

        //- Normalize fluxes
        void normalize(scalar);


public:

    // Public Data


    //- Runtime type information
    TypeName("SN");


    // Constructors

        //- Construct from region mesh
        SNSolver(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        SNSolver(const SNSolver&) = delete;


    //- Destructor
    virtual ~SNSolver();


    // Member Functions

        //- Solve neutron SN equations
        void solve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SNSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //