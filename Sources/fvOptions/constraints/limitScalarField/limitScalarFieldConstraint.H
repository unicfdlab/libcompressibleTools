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
    Foam::fv::limitScalarFieldConstraint

Description
    Limits scalar field values by specifying minimum boundary.

    Sources described by:

        limitScalarFieldConstraintCoeffs
        {
            mode            uniform;    // uniform or lookup

            // uniform option
            minFieldValue     constant 500; // fixed temperature with time 500[K]
            
            scaleCoeff       1;

            // lookup option
            fieldName        T;          // optional scalar field name
        }

Note:
    The 'uniform' option allows the use of a time-varying uniform temperature
    by means of the DataEntry type.

SourceFiles
    fvOption.C

\*---------------------------------------------------------------------------*/

#ifndef limitScalarFieldConstraint_H
#define limitScalarFieldConstraint_H

#include "cellSetOption.H"
#include "Enum.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class limitScalarFieldConstraint Declaration
\*---------------------------------------------------------------------------*/

class limitScalarFieldConstraint
:
    public cellSetOption
{

protected:

    // Protected data

        //- Uniform value
        autoPtr<Function1<scalar> > fieldValue_;
        
        //-
        scalar scaleCoeff_;

        //- Scalar field name
        word fieldName_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        limitScalarFieldConstraint(const limitScalarFieldConstraint&);

        //- Disallow default bitwise assignment
        void operator=(const limitScalarFieldConstraint&);


public:

    //- Runtime type information
    TypeName("limitScalarField");


    // Constructors

        //- Construct from components
        limitScalarFieldConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~limitScalarFieldConstraint()
    {}


    // Member Functions

        //- Constrain energy equation to fix the temperature
        virtual void constrain(fvMatrix<scalar>& eqn, const label fieldI);


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
