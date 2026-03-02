/*---------------------------------------------------------------------------*\

    MUSCLE-Foam: an OpenFOAM-based multiphysics solver for reactor physics
    Copyright (C) 2026 NRG PALLAS

\-----------------------------------------------------------------------------/

License
    This file is part of MUSCLE-Foam. It is published under the terms of the
    GNU General Public License. See <http://www.gnu.org/licenses/>.

Class
    Foam::neutronicsSolver

Description
    Abstract base class for run-time selectable neutronics solvers.

SourceFiles
    neutronicsSolver.C
    neutronicsSolverNew.C

\*---------------------------------------------------------------------------*/

#ifndef neutronicsSolver_H
#define neutronicsSolver_H

#include "fvMesh.H"
#include "pimpleNoLoopControl.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "regionName.H"
#include "DNPTransport.H"
#include "DHPTransport.H"
#include "IOdictionary.H"
#include "fluidProperties.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class neutronicsSolver Declaration
\*---------------------------------------------------------------------------*/

class neutronicsSolver
:
    public regIOobject
{
protected:

    // Protected Data

        //- Region mesh
        fvMesh& mesh_;

        //- Switch for steady-state
        //  Set automatically from the selected ddt scheme
        bool steady;

        //- Switch for local time step transient operation
        //  Set automatically from the selected ddt scheme
        bool LTS;

        //- Pointer to the demand driven fvModels MeshObject
        mutable Foam::fvModels* fvModelsPtr;

        //- Pointer to the demand driven fvConstraints MeshObject
        mutable Foam::fvConstraints* fvConstraintsPtr;

        //- (Mixture) fluid properties needed for neutronics solver
        fluidProperties fluid_;

        //- Pointer to delayed neutron transport object
        autoPtr<DNPTransport> DNP_;

        //- Switch for decay heat
        bool decayHeat_;

        //- Pointer to decay heat precursor transport (if any)
        autoPtr<DHPTransport> DHP_;

        //- Number of neutron energy groups
        label energyGroups_;

        //- Effective multiplication factor for steady-state calculations
        scalar k_;

        //- Number of timesteps for which k is kept constant at k=k0
        label criticalityInactive_;

        //- Timestep counter before criticality calculation
        label criticalityInactiveCounter_;

        //- Switch for criticality calculation
        bool criticalityCalculationActive_;

        //- Nominal power of the system
        scalar nominalPower_;

        //- Neutron fission source
        volScalarField neutronFissionSource_;

        //- Power source
        volScalarField powerSource_;

        //- Prompt heat source
        volScalarField promptHeatSource_;

        //- Total power density
        volScalarField powerDensity_;

        //- Density field used for cross-section calculation
        volScalarField rhoXS_;

        //- Temperature field used for cross-section calculation
        volScalarField TXS_;

        //- Neutron fission sources per energy group
        PtrList<volScalarField> FNS_;

        //- Scattering neutron sources per energy group
        PtrList<volScalarField> SNS_;


    // Protected Member Functions

        //- Dummy write for regIOobject
        virtual bool writeData(Ostream&) const;

        //- Bool for whether or not criticality should be calculated
        //  in the current timestep
        bool criticalityCalculationActive();

        //- Set delayed neutron transport
        void setDNPTransport();

        //- Set decay heat transport (if any)
        void setDHPTransport();

        //- Set neutron sources
        void setNeutronSources();

        //- Read controls
        virtual bool read();

        //- Normalize fluxes and power
        virtual void normalize(const scalar);


public:

    // Public Data

        //- Region mesh
        const fvMesh& mesh;

        //- Time
        const Time& runTime;

        //- Effective multiplication factor
        const scalar& k;

        //- (Mixture) fluid properties needed for neutronics solver
        const fluidProperties& fluid;

    //- Runtime type information
    TypeName("neutronicsSolver");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            neutronicsSolver,
            fvMesh,
            (fvMesh& mesh),
            (mesh)
        );


    // Constructors

        //- Construct from region mesh
        neutronicsSolver(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        neutronicsSolver(const neutronicsSolver&) = delete;


    // Selectors

        //- Select, construct and return the neutronicsSolver
        static autoPtr<neutronicsSolver> New
        (
            const word& neutronicsSolverName,
            fvMesh& mesh
        );


    //- Destructor
    virtual ~neutronicsSolver();


    // Member Functions

        //- Return the fvModels that are created on demand
        Foam::fvModels& fvModels() const;

        //- Return the fvConstraints that are created on demand
        Foam::fvConstraints& fvConstraints() const;

        //- Solve the neutron flux equations
        virtual void solve() = 0;

        //- Transport delayed neutron and decay heat precursors if needed
        void precursorTransport();

        //- Perform criticality calculation for steady-state calculations
        void criticalityCalculation();

        //- Update the power density
        void updatePower();

        //- Return number of energy groups
        inline label energyGroups() { return energyGroups_; }

        //- Return reference to mesh
        inline const fvMesh& fvMsh() { return mesh; }

        //- Reference to neutron fission source
        inline const volScalarField& neutronFissionSource() { return neutronFissionSource_; }

        //- Reference to power source
        inline const volScalarField& powerSource() { return powerSource_; }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const neutronicsSolver&) = delete;
};


template<>
inline const word& regionName(const neutronicsSolver& region)
{
    return regionName(region.mesh);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //