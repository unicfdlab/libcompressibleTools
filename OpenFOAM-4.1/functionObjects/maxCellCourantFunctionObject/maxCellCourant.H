/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::maxCellCourant

Group
    grpmaxCellCourantFunctionObjects

Description
    

SeeAlso
    

SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef maxCellCourant_H
#define maxCellCourant_H

#include "typeInfo.H"
#include "regionFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class maxCellCourant Declaration
\*---------------------------------------------------------------------------*/

class maxCellCourant
: public regionFunctionObject
{

private:

    //-
    word CourantType_;
    
    //-
    word rhoName_;
    
    //-
    word phiName_;

private:
    
    //- Disallow default bitwise copy construct
    maxCellCourant(const maxCellCourant&);

    //- Disallow default bitwise assignment
    void operator=(const maxCellCourant&);

public:

    //- Runtime type information
    TypeName("maxCellCourant");

    // Constructors

        //- Construct from Time and dictionary
        maxCellCourant
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from objectRegistry and dictionary
        maxCellCourant
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

    //- Destructor
        virtual ~maxCellCourant();


    // Member Functions

        //- Read the maxCellCourant data
        virtual bool read(const dictionary&);
        
        //-
        virtual bool execute();

        //- Write the maxCellCourant (write forces output to console)
        virtual bool write();
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
