/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::normalDirectionInletOutletFvPatchField

Group
    grpOutletBoundaryConditions

Description
    Variant of inletOutlet boundary condition with uniform inletValue depending on
    direction of reference vector.

Usage

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type                normalDirectionInletOutlet;
        phi                 phiJm;
        Uref                (-3 0 0);
        uniformInletValue   0;
        value               uniform 0;
    }
    \endverbatim

    The mode of operation is determined by the sign of the flux across the
    patch faces.

Note
    Sign conventions:
    - positive flux (out of domain): apply zero-gradient condition
    - negative flux (into of domain): apply the user-specified fixed value

See also
    Foam::inletOutletFvPatchField
    Foam::Function1Types

SourceFiles
    normalDirectionInletOutletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef normalDirectionInletOutletFvPatchField_H
#define normalDirectionInletOutletFvPatchField_H

#include "mixedFvPatchField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class normalDirectionInletOutletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class normalDirectionInletOutletFvPatchField
:
    public mixedFvPatchField<Type>
{

protected:

    // Protected data

        //- Name of flux field
        vector Uref_;
        
        //-
        scalar tol_;
        
        //-
        scalarField phiRef_;
        
        //- Value
        autoPtr<Function1<Type>> uniformInletValue_;


public:

    //- Runtime type information
    TypeName("normalDirectionInletOutlet");


    // Constructors

        //- Construct from patch and internal field
        normalDirectionInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        normalDirectionInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given normalDirectionInletOutletFvPatchField
        //  onto a new patch
        normalDirectionInletOutletFvPatchField
        (
            const normalDirectionInletOutletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        normalDirectionInletOutletFvPatchField
        (
            const normalDirectionInletOutletFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new normalDirectionInletOutletFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        normalDirectionInletOutletFvPatchField
        (
            const normalDirectionInletOutletFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new normalDirectionInletOutletFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Attributes

            //- Return true: this patch field is altered by assignment
            virtual bool assignable() const
            {
                return true;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<Type>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "normalDirectionInletOutletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
