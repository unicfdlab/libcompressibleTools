/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
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
    Foam::fv::setComponentValue

Description
    Set component, specified by direction to value dirValue


SourceFiles
    setComponentValue.C

\*---------------------------------------------------------------------------*/

#ifndef setComponentValue_H
#define setComponentValue_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class setComponentValue Declaration
\*---------------------------------------------------------------------------*/

class setComponentValue
:
    public cellSetOption
{

protected:

    // Protected data

        //- Direction
        label direction_;
        
        //-
        scalar dirValue_;
        
        //-
        word vectorName_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        setComponentValue(const setComponentValue&);

        //- Disallow default bitwise assignment
        void operator=(const setComponentValue&);


public:

    //- Runtime type information
    TypeName("setComponentValue");


    // Constructors

        //- Construct from components
        setComponentValue
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~setComponentValue()
    {}


    // Member Functions

        // Evaluate

            //- Correct the energy field
            virtual void correct(volVectorField& U);


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
