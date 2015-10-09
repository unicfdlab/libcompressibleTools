/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "subsonicSupersonicPressureOutletFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subsonicSupersonicPressureOutletFvPatchScalarField::
subsonicSupersonicPressureOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    psiName_("psi"),
    phiName_("phi"),
    UName_("U"),
    p0_(p.size(), 0.0)
{}

Foam::subsonicSupersonicPressureOutletFvPatchScalarField::
subsonicSupersonicPressureOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict),
    psiName_(dict.lookup("psi")),
    phiName_(dict.lookup("phi")),
    UName_(dict.lookup("U")),
    p0_("p0", dict, p.size())
{
    if (!dict.found("value"))
    {
	fvPatchField<scalar>::operator=(p0_);
    }
}

Foam::subsonicSupersonicPressureOutletFvPatchScalarField::
subsonicSupersonicPressureOutletFvPatchScalarField
(
    const subsonicSupersonicPressureOutletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    psiName_(ptf.psiName_),
    phiName_(ptf.phiName_),
    UName_(ptf.UName_),
    p0_(ptf.p0_, mapper)
{}

Foam::subsonicSupersonicPressureOutletFvPatchScalarField::
subsonicSupersonicPressureOutletFvPatchScalarField
(
    const subsonicSupersonicPressureOutletFvPatchScalarField& wbppsf
)
:
    mixedFvPatchScalarField(wbppsf),
    psiName_(wbppsf.psiName_),
    phiName_(wbppsf.phiName_),
    UName_(wbppsf.UName_),
    p0_(wbppsf.p0_)
{}


Foam::subsonicSupersonicPressureOutletFvPatchScalarField::
subsonicSupersonicPressureOutletFvPatchScalarField
(
    const subsonicSupersonicPressureOutletFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wbppsf, iF),
    psiName_(wbppsf.psiName_),
    phiName_(wbppsf.phiName_),
    UName_(wbppsf.UName_),
    p0_(wbppsf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::subsonicSupersonicPressureOutletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    //-
    const basicThermo& thermo = db().lookupObject<basicThermo>("thermophysicalProperties");
    
    //
    const scalarField gammap = thermo.gamma()().boundaryField()[patch().index()];
    
    //
    const scalarField psip   = db().lookupObject<volScalarField>(psiName_).boundaryField()[patch().index()];
    
    //-
    const vectorField Up     = db().lookupObject<volVectorField>(UName_).boundaryField()[patch().index()];
    
    //
    const scalarField c = sqrt(gammap / psip);
    
    //-
    const scalarField phip   = db().lookupObject<surfaceScalarField>(phiName_).boundaryField()[patch().index()];
    
    forAll(this->valueFraction(), iFace)
    {
	if ( mag(Up[iFace]) < c[iFace] ) //subsonic flow
	{
	    this->valueFraction()[iFace] = 1.0;
	    this->refValue()[iFace] = 
		p0_[iFace] - this->operator[](iFace)*psip[iFace]*magSqr(Up[iFace])*0.5*(1.0 - pos(phip[iFace]));
	    this->refGrad()[iFace] = 0.0;
	}
	else
	{
	    this->valueFraction()[iFace] = 0.0;
	    this->refValue()[iFace] = p0_[iFace];
	    this->refGrad()[iFace] = 0.0;
	}
    }


//    if (phi.dimensions() == dimVelocity*dimArea)
//    {
//        phip = patch().Sf() & Up.freestreamValue();
//    }
//    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
//    {
//        const fvPatchField<scalar>& rhop =
//            patch().lookupPatchField<volScalarField, scalar>("rho");
//
//        phip = rhop*(patch().Sf() & Up.freestreamValue());
//    }
//    else
//    {
//        FatalErrorIn("subsonicSupersonicPressureOutletFvPatchScalarField::updateCoeffs()")
//            << "dimensions of phi are not correct"
//            << "\n    on patch " << this->patch().name()
//            << " of field " << this->dimensionedInternalField().name()
//            << " in file " << this->dimensionedInternalField().objectPath()
//            << exit(FatalError);
//    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::subsonicSupersonicPressureOutletFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    p0_.autoMap(m);
}

void Foam::subsonicSupersonicPressureOutletFvPatchScalarField::rmap
(
	const fvPatchScalarField& ptf,
	const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const subsonicSupersonicPressureOutletFvPatchScalarField& tiptf =
	refCast<const subsonicSupersonicPressureOutletFvPatchScalarField>(ptf);

    p0_.rmap(tiptf.p0_, addr);
}



void Foam::subsonicSupersonicPressureOutletFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("psi") << psiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
    p0_.writeEntry("p0", os);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


void Foam::subsonicSupersonicPressureOutletFvPatchScalarField::operator=
(
    const fvPatchField<scalar>& ptf
)
{
    mixedFvPatchScalarField::operator=
    (
	ptf
    );
    
//    p0_ = ptf.p0_;
//    psiName_ = ptf.psiName_;
//    phiName_ = ptf.phiName_;
//    UName_   = ptf.UName_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        subsonicSupersonicPressureOutletFvPatchScalarField
    );
}

// ************************************************************************* //
