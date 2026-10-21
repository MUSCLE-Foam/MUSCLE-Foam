/*---------------------------------------------------------------------------*\

    MUSCLE-Foam: an OpenFOAM-based multiphysics solver for reactor physics
    Copyright (C) 2026 NRG PALLAS

\-----------------------------------------------------------------------------/

License
    This file is part of MUSCLE-Foam. It is published under the terms of the
    GNU General Public License. See <http://www.gnu.org/licenses/>.

Class
    Foam::PKESolver

Description
    Neutronics solver based on neutron PKE theory.

    Reference:
    \verbatim
        Stacey, W.M. (2018)
        Nuclear Reactor Physics, Chapter 3.
        John Wiley & Sons.
    \endverbatim

SourceFiles
    PKESolver.C
    PKESolverInit.C

\*---------------------------------------------------------------------------*/

#ifndef PKESolver_H
#define PKESolver_H

#include "neutronicsSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class PKESolver Declaration
\*---------------------------------------------------------------------------*/

class PKESolver
:
    public neutronicsSolver
{
private:

    // Private Member Data

        //- Kinetics constants dictionary
        IOdictionary kineticsConstants_;

        //- Mean prompt generation time
        dimensionedScalar L_;

        //- System reactivity
        scalar reactivity_;

        //- Extra reactivity for moving fuel
        scalar reactivityPlus_;

        //- Normalized neutron flux
        scalar psi_;

        //- Normalized precursor concentrations
        scalarList eta_;

        //- Moving DNPs
        Switch movingPrecursors_;

        //- Reactivity control system
        Switch controlSystem_;

        //- Precursor in-core time
        dimensionedScalar tauc_;

        //- Precursor ex-core time
        dimensionedScalar taul_;

        //- Temperature feedback coefficient
        dimensionedScalar alphaT_;

        //- Density feedback coefficient
        dimensionedScalar alphaRho_;

        //- Cells in core cellZone
        const labelList& coreCells_;

        //- Core volume
        scalar coreVolume_;

        //- Core initial temperature
        scalar coreInitTemp_;

        //- Core initial density
        scalar coreInitRho_;

    // Private Member Functions

        //- Normalize fluxes
        void normalize(scalar);


public:

    // Public Data


    //- Runtime type information
    TypeName("PKE");


    // Constructors

        //- Construct from region mesh
        PKESolver(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        PKESolver(const PKESolver&) = delete;


    //- Destructor
    virtual ~PKESolver();


    // Member Functions

        //- Solve neutron PKE equations
        void solve();

        //- Update power for PKE solver
        void updatePower() override;

        //- Compute a core volume average
        scalar coreAverage(const volScalarField&);

        //- Access reactivity
        inline scalar& reactivity() { return reactivity_; }

        //- Access beta tot
        inline scalar betaTotal() { return this->DNP_->betaTotal(); }

        //- Access tauc
        inline dimensionedScalar& tauc() { return tauc_; }

        //- Access tauc
        inline dimensionedScalar& taul() { return taul_; }

        //- Access kineticsDict
        inline IOdictionary& kineticsConstants() { return kineticsConstants_; }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const PKESolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //