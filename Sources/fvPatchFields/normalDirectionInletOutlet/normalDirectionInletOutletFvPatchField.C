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

\*---------------------------------------------------------------------------*/

#include "normalDirectionInletOutletFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::normalDirectionInletOutletFvPatchField<Type>::normalDirectionInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    Uref_(vector::zero),
    tol_(0.0),
    phiRef_(p.size(), 0.0)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::normalDirectionInletOutletFvPatchField<Type>::normalDirectionInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    Uref_(dict.lookup("Uref")),
    tol_(readScalar(dict.lookup("tol"))),
    phiRef_((p.Sf()/p.magSf()) & Uref_),
    uniformInletValue_(Function1<Type>::New("uniformInletValue", dict))
{
    this->refValue() =
        uniformInletValue_->value(this->db().time().timeOutputValue());

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(this->refValue());
    }

    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::normalDirectionInletOutletFvPatchField<Type>::normalDirectionInletOutletFvPatchField
(
    const normalDirectionInletOutletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(p, iF),  // Don't map
    Uref_(ptf.Uref_),
    tol_(ptf.tol_),
    phiRef_((p.Sf()/p.magSf()) & Uref_),
    uniformInletValue_(ptf.uniformInletValue_)
{
    this->patchType() = ptf.patchType();

    // Evaluate refValue since not mapped
    this->refValue() =
        uniformInletValue_->value(this->db().time().timeOutputValue());

    this->refGrad() = Zero;
    this->valueFraction() = 0.0;

    // Initialize the patch value to the refValue
    fvPatchField<Type>::operator=(this->refValue());

    this->map(ptf, mapper);
}


template<class Type>
Foam::normalDirectionInletOutletFvPatchField<Type>::normalDirectionInletOutletFvPatchField
(
    const normalDirectionInletOutletFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    Uref_(ptf.Uref_),
    tol_(ptf.tol_),
    phiRef_((this->patch().Sf() / this->patch().magSf()) & Uref_),
    uniformInletValue_(ptf.uniformInletValue_)
{}


template<class Type>
Foam::normalDirectionInletOutletFvPatchField<Type>::normalDirectionInletOutletFvPatchField
(
    const normalDirectionInletOutletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    Uref_(ptf.Uref_),
    tol_(ptf.tol_),
    phiRef_(this->patch().nf() & Uref_),
    uniformInletValue_(ptf.uniformInletValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::normalDirectionInletOutletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    
//    scalarField phinb = 
//    this->patch().nf() & 
//    this->patch().template
//        lookupPatchField<volVectorField,vector>("U");
    
    this->refValue() =
        uniformInletValue_->value(this->db().time().timeOutputValue());
    
    const Field<scalar>& phip =
        phiRef_;
    
    this->valueFraction() = neg(phip + tol_); //zero gradient for zero flux
//    this->valueFraction() *= neg(phinb);
    
    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::normalDirectionInletOutletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("Uref") << Uref_ << token::END_STATEMENT << nl;
    os.writeKeyword("tol") << tol_ << token::END_STATEMENT << nl;
    
    this->uniformInletValue_->writeData(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::normalDirectionInletOutletFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchField<Type>::autoMap(m);

    // Override
    this->refValue() =
        uniformInletValue_->value(this->db().time().timeOutputValue());
}


template<class Type>
void Foam::normalDirectionInletOutletFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<Type>::rmap(ptf, addr);

    // Override
    this->refValue() =
        uniformInletValue_->value(this->db().time().timeOutputValue());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::normalDirectionInletOutletFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator=
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );
}


// ************************************************************************* //
