/*---------------------------------------------------------------------------*\

    MUSCLE-Foam: an OpenFOAM-based multiphysics solver for reactor physics
    Copyright (C) 2026 NRG PALLAS

\-----------------------------------------------------------------------------/

License
    This file is part of MUSCLE-Foam. It is published under the terms of the
    GNU General Public License. See <http://www.gnu.org/licenses/>.

Class
    Foam::diffusionSolver

Description
    Neutronics solver based on neutron diffusion theory.

    The solver computes neutron fluxes for a user-specified number of energy
    groups. Cross-sections and thermal feedbacks should be determined in
    advance, e.g., using precursor Monte Carlo simulation.

    Reference:
    \verbatim
        Stacey, W.M. (2018)
        Nuclear Reactor Physics, Chapter 3.
        John Wiley & Sons.
    \endverbatim

SourceFiles
    diffusionSolver.C
    diffusionSolverInit.C

\*---------------------------------------------------------------------------*/

#ifndef diffusionSolver_H
#define diffusionSolver_H

#include "neutronicsSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class diffusionSolver Declaration
\*---------------------------------------------------------------------------*/

class diffusionSolver
:
    public neutronicsSolver
{
private:

    // Private Member Data

        //- Neutron fluxes
        PtrList<volScalarField> flux_;

        //- Albedo flux
        volScalarField fluxAlbedo_;

        //- Inverse neutron velocities
        PtrList<volScalarField> IV_;

        //- Prompt neutron yields
        PtrList<volScalarField> chiPrompt_;

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

        //- Isotropic scattering cross sections
        List<PtrList<volScalarField>> S_;
        List<PtrList<volScalarField>> S0_;
        List<PtrList<volScalarField>> alphaS_;

        //- Normalized prompt neutron sources
        PtrList<volScalarField> F_;

        //- Diffusion coefficient
        volScalarField DAlbedo_;

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
    TypeName("diffusion");


    // Constructors

        //- Construct from region mesh
        diffusionSolver(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        diffusionSolver(const diffusionSolver&) = delete;


    //- Destructor
    virtual ~diffusionSolver();


    // Member Functions

        //- Solve neutron diffusion equations
        void solve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const diffusionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //