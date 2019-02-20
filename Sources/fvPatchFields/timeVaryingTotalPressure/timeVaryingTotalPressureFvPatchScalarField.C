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

#include "timeVaryingTotalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidThermo.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeVaryingTotalPressureFvPatchScalarField::timeVaryingTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    rhoName_("none"),
    useGamma_(false),
    compressible_(false),
    p0_()
{}


Foam::timeVaryingTotalPressureFvPatchScalarField::timeVaryingTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "none")),
    useGamma_((dict.lookup("useGamma"))),
    compressible_((dict.lookup("compressible"))),
    p0_
    (
        Function1<scalar>::New ("p0", dict)
    )
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        const scalar ct = db().time().timeOutputValue();
        scalarField p0 (this->size(), p0_->value(ct));
        this->operator == (p0);
    }
    
}


Foam::timeVaryingTotalPressureFvPatchScalarField::timeVaryingTotalPressureFvPatchScalarField
(
    const timeVaryingTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    useGamma_(ptf.useGamma_),
    compressible_(ptf.compressible_),
    p0_(ptf.p0_().clone().ptr())
{}


Foam::timeVaryingTotalPressureFvPatchScalarField::timeVaryingTotalPressureFvPatchScalarField
(
    const timeVaryingTotalPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    useGamma_(tppsf.useGamma_),
    compressible_(tppsf.compressible_),
    p0_(tppsf.p0_().clone().ptr())
{}


Foam::timeVaryingTotalPressureFvPatchScalarField::timeVaryingTotalPressureFvPatchScalarField
(
    const timeVaryingTotalPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    useGamma_(tppsf.useGamma_),
    compressible_(tppsf.compressible_),
    p0_(tppsf.p0_().clone().ptr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeVaryingTotalPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}


void Foam::timeVaryingTotalPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);
}


void Foam::timeVaryingTotalPressureFvPatchScalarField::updateCoeffs
(
    const scalarField& p0p,
    const vectorField& Up
)
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    if (!compressible_ && rhoName_ == "none")
    {
        operator==(p0p - 0.5*(1.0 - pos(phip))*magSqr(Up));
    }
    else if (compressible_)
    {
        const fluidThermo& thermo = 
                this->patch().boundaryMesh().mesh().thisDb().lookupObject<fluidThermo>("thermophysicalProperties");

        const fvPatchField<scalar>& psip = 
            thermo.psi().boundaryField()[patch().index()];
        
        if (useGamma_)
        {

            const fvPatchField<scalar>& pp =
                thermo.p().boundaryField()[patch().index()];
        
            const fvPatchField<scalar>& Tp =
                thermo.T().boundaryField()[patch().index()];
        
            scalarField gammap
            (
                thermo.gamma(pp, Tp, patch().index())
            );
            
            scalarField gM1ByG = (gammap - 1.0)/gammap;

            operator==
            (
                p0p
               /pow
                (
                    (1.0 + 0.5*psip*gM1ByG*(1.0 - pos(phip))*magSqr(Up)),
                    1.0/gM1ByG
                )
            );
        }
        else
        {
            operator==(p0p/(1.0 + 0.5*psip*(1.0 - pos(phip))*magSqr(Up)));
        }
    }
    else if (!compressible_ && rhoName_ != "none")
    {
        const fvPatchField<scalar>& rho =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);
        operator==(p0p - 0.5*rho*(1.0 - pos(phip))*magSqr(Up));
    }
    else
    {
        FatalErrorIn
        (
            "timeVaryingTotalPressureFvPatchScalarField::updateCoeffs()"
        )   << " rho or compressibility set inconsistently, rho = " << rhoName_
            << ", compressibility = " << compressible_ << ".\n"
            << "    Set either rho or compressible_ or neither depending on the "
               "definition of total pressure." << nl
            << "    Set the unused variable(s) to 'none'.\n"
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }
    
    Info << "Max/min of pressure at patch is: " << gMax(*this) << "/" << gMin(*this) << endl;

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::timeVaryingTotalPressureFvPatchScalarField::updateCoeffs()
{

    const scalar ct = db().time().timeOutputValue();
    
    const scalarField p0
    (
	this->size(),
	p0_->value(ct)
    );
    
    updateCoeffs
    (
        p0,
        patch().lookupPatchField<volVectorField, vector>(UName())
    );
}


void Foam::timeVaryingTotalPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    os.writeKeyword("compressible") << compressible_ << token::END_STATEMENT << nl;
    os.writeKeyword("useGamma") << useGamma_ << token::END_STATEMENT << nl;
    p0_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        timeVaryingTotalPressureFvPatchScalarField
    );
}

// ************************************************************************* //
