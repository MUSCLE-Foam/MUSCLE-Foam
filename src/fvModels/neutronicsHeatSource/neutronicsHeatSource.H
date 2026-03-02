/*---------------------------------------------------------------------------*\

    MUSCLE-Foam: an OpenFOAM-based multiphysics solver for reactor physics
    Copyright (C) 2026 NRG PALLAS

\-----------------------------------------------------------------------------/

License
    This file is part of MUSCLE-Foam. It is published under the terms of the
    GNU General Public License. See <http://www.gnu.org/licenses/>.

Class
    Foam::fv::neutronicsHeatSource

Description
    fvModel for the inclusion of a heat source based on the power density
    calculated by the neutronics solver.

    An optional heat sink can also be included to model a heat exchanger.
    The heat sink is active in the chosen selection region (e.g., a cellSet)
    and is calculated as follows:

    \f[
        Q = \gamma * (T - T_{ref})
    \f]

    where
    \vartable
        Q       | Optional heat sink [W/m^3]
        \gamma  | heat sink coefficient [W/(m^3*K)]
        T       | Temperature [K]
        T_{ref} | Reference temperature [K]
    \endvartable

Usage
    Example usage:
    \verbatim
    neutronicsHeatSource
    {
        type            neutronicsHeatSource;

        select          cellSet;
        cellSet         heatExchanger;

        // Optional entries

        gamma           1e6;
        Tref            900;
    }
    \endverbatim

SourceFiles
    neutronicsHeatSource.C

\*---------------------------------------------------------------------------*/

#ifndef neutronicsHeatSource_H
#define neutronicsHeatSource_H

#include "fvModel.H"
#include "dimensionedType.H"
#include "fvCellZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                Class neutronicsHeatSource Declaration
\*---------------------------------------------------------------------------*/

class neutronicsHeatSource
:
    public fvModel
{
    // Private data

        //- Optional heat sink heat transfer coefficient
        dimensionedScalar gamma_;

        //- Reference temperature for optional heat sink
        dimensionedScalar Tref_;

        //- Optional heat exchanger cellSet
        fvCellZone zone_;

        //- List of fields for which a source term is added
        wordList heNames_;

public:

    //- Runtime type information
    TypeName("neutronicsHeatSource");


    // Constructors

        //- Construct from dictionary
        neutronicsHeatSource
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~neutronicsHeatSource();


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Sources

            //- Source term to energy equation
            virtual void addSup
            (
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;

            //- Source term to compressible energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;

            //- Source term to phase energy equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //