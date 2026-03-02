/*---------------------------------------------------------------------------*\

    MUSCLE-Foam: an OpenFOAM-based multiphysics solver for reactor physics
    Copyright (C) 2026 NRG PALLAS

\-----------------------------------------------------------------------------/

License
    This file is part of MUSCLE-Foam. It is published under the terms of the
    GNU General Public License. See <http://www.gnu.org/licenses/>.

Class
    Foam::SP3Solver

Description
    Neutronics solver based on the Simplified P_N equations with N=3 and using
    diffusion-like within-group form.

    The solver computes the zeroth and second order neutron flux moments for a
    user-specified number of energy groups. Cross-sections and thermal feed-
    backs should be determined in advance, e.g., using precursor Monte Carlo
    simulation.

    Reference:
    \verbatim
        Cervi, E., Lorenzi, S., Cammi, A. & Luzzi, L. (2019)
        Development of an SP3 neutron transport solver for the analysis of the
        Molten Salt Fast Reactor.
        Nuclear Engineering and Design, 346, 209-219.
    \endverbatim

SourceFiles
    SP3Solver.C
    SP3SolverInit.C

\*---------------------------------------------------------------------------*/

#ifndef SP3Solver_H
#define SP3Solver_H

#include "neutronicsSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class SP3Solver Declaration
\*---------------------------------------------------------------------------*/

class SP3Solver
:
    public neutronicsSolver
{
private:

    // Private Member Data

        //- Neutron fluxes
        PtrList<volScalarField> flux_;

        //- Second moment neutron fluxes
        PtrList<volScalarField> secondFlux_;

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

        //- Anisotropic scattering cross sections
        List<PtrList<volScalarField>> S2_;
        List<PtrList<volScalarField>> S20_;
        List<PtrList<volScalarField>> alphaS2_;

        List<PtrList<volScalarField>> S3_;
        List<PtrList<volScalarField>> S30_;
        List<PtrList<volScalarField>> alphaS3_;

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
    TypeName("SP3");


    // Constructors

        //- Construct from region mesh
        SP3Solver(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        SP3Solver(const SP3Solver&) = delete;


    //- Destructor
    virtual ~SP3Solver();


    // Member Functions

        //- Solve neutron SP3 equations
        void solve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SP3Solver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //