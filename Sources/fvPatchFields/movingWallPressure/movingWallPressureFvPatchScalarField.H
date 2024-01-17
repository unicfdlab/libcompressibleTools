/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::subsonicSupersonicPressureOutletFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides a free-stream condition for pressure where flow 
    can be subsonic or supersonic. If flux is supersonic, then zero gradient condition
    is applied, otherwise - total presure
    
    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            movingWallPressure;
    }
    \endverbatim

Note
    This condition is designed to operate with a zero gradient velocity condition

SeeAlso
    Foam::zeroGradientFvPatchField
    Foam::freestreamFvPatchField

SourceFiles
    movingWallPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef movingWallPressureFvPatchScalarFields_H
#define movingWallPressureFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class movingWallPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class movingWallPressureFvPatchScalarField
:
    public mixedFvPatchScalarField
{

private:

    //-
    word psiName_;
    
    //-
    word phiName_;
    
    //-
    word UName_;
    
    //-
    scalarField p0_;

public:

    //- Runtime type information
    TypeName("movingWallPressure");


    // Constructors

        //- Construct from patch and internal field
        movingWallPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        movingWallPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given movingWallPressureFvPatchScalarField onto
        //  a new patch
        movingWallPressureFvPatchScalarField
        (
            const movingWallPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        movingWallPressureFvPatchScalarField
        (
            const movingWallPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new movingWallPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        movingWallPressureFvPatchScalarField
        (
            const movingWallPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new movingWallPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
            
            //-
            virtual void write(Ostream& os) const;
            
            //-
            virtual void operator=(const fvPatchField<scalar>& pvf);


        // Mapping functions
        
	    //- Map (and resize as needed) from self given a mapping object
	    virtual void autoMap
	    (
		const fvPatchFieldMapper&
	    );

	    //- Reverse map the given fvPatchField onto this fvPatchField
	    virtual void rmap
	    (
		const fvPatchScalarField&,
		const labelList&
	    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
