/*---------------------------------------------------------------------------*\

Class
    Foam::functionObjects::reactivityInsertion

Description
    Reactivity insertion

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    reactivityInsertion.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_reactivityInsertion_H
#define functionObjects_reactivityInsertion_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"
#include "dictionary.H"
#include "PKESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class reactivityInsertion Declaration
\*---------------------------------------------------------------------------*/

class reactivityInsertion
:
    public fvMeshFunctionObject
{

    //- Private member data

        // Reactivity insertion
        scalar rho_;

public:

    //- Runtime type information
    TypeName("reactivityInsertion");


    // Constructors

        //- Construct from Time and dictionary
        reactivityInsertion
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        reactivityInsertion(const reactivityInsertion&) = delete;


    //- Destructor
    virtual ~reactivityInsertion();


    // Member Functions

        //- Read the reactivityInsertion data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Calculate the reactivityInsertion field
        virtual bool execute();

        //- Write the reactivityInsertion field
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reactivityInsertion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //