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

#include "flowRateControlledWithPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateControlledWithPressureFvPatchScalarField::flowRateControlledWithPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    rhoName_("none"),
    psiName_("none"),
    gamma_(0.0),
    p0_(0.0),
    flowRate_(),
    flowController_()
{}


Foam::flowRateControlledWithPressureFvPatchScalarField::flowRateControlledWithPressureFvPatchScalarField
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
    psiName_(dict.lookupOrDefault<word>("psi", "none")),
    gamma_(readScalar(dict.lookup("gamma"))),
    p0_(readScalar(dict.lookup("p0"))),
    flowRate_
    (
	Function1<scalar>::New ("flowRate", dict)
    ),
    flowController_
    (
	Foam::fv::controllerModel::New
	(
	    "flowRateController" + p.name(),
	    dict.subDict("flowRateController"),
	    p.boundaryMesh().mesh().time()
	)
    )
{
    scalarField p0 (this->size(), p0_);
    this->operator == (p0);
}


Foam::flowRateControlledWithPressureFvPatchScalarField::flowRateControlledWithPressureFvPatchScalarField
(
    const flowRateControlledWithPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    p0_(ptf.p0_),
    flowRate_(ptf.flowRate_().clone().ptr()),
    flowController_(ptf.flowController_().clone().ptr())
{}


Foam::flowRateControlledWithPressureFvPatchScalarField::flowRateControlledWithPressureFvPatchScalarField
(
    const flowRateControlledWithPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    p0_(tppsf.p0_),
    flowRate_(tppsf.flowRate_().clone().ptr()),
    flowController_(tppsf.flowController_().clone().ptr())
{}


Foam::flowRateControlledWithPressureFvPatchScalarField::flowRateControlledWithPressureFvPatchScalarField
(
    const flowRateControlledWithPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    p0_(tppsf.p0_),
    flowRate_(tppsf.flowRate_().clone().ptr()),
    flowController_(tppsf.flowController_().clone().ptr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::flowRateControlledWithPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}


void Foam::flowRateControlledWithPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);
}


void Foam::flowRateControlledWithPressureFvPatchScalarField::updateCoeffs
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

    if (psiName_ == "none" && rhoName_ == "none")
    {
        operator==(p0p - 0.5*(1.0 - pos(phip))*magSqr(Up));
    }
    else if (rhoName_ == "none")
    {
        const fvPatchField<scalar>& psip =
            patch().lookupPatchField<volScalarField, scalar>(psiName_);

        if (gamma_ > 1.0)
        {
            scalar gM1ByG = (gamma_ - 1.0)/gamma_;

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
    else if (psiName_ == "none")
    {
        const fvPatchField<scalar>& rho =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);
        operator==(p0p - 0.5*rho*(1.0 - pos(phip))*magSqr(Up));
    }
    else
    {
        FatalErrorIn
        (
            "flowRateControlledWithPressureFvPatchScalarField::updateCoeffs()"
        )   << " rho or psi set inconsistently, rho = " << rhoName_
            << ", psi = " << psiName_ << ".\n"
            << "    Set either rho or psi or neither depending on the "
               "definition of total pressure." << nl
            << "    Set the unused variable(s) to 'none'.\n"
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::flowRateControlledWithPressureFvPatchScalarField::updateCoeffs()
{

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const scalar ct = db().time().timeOutputValue();
    
    const scalar refFlow = flowRate_->value(ct);
    
    const scalar actualFlow = -gSum(phip);
    
    const scalar p0Increment = 
        flowController_->newIncrement
        (
            flowController_->error(refFlow, actualFlow),
            ct
        );
    
    p0_ += p0Increment;
    
    Info << "At tune - " << ct << " actual flx = " << actualFlow << endl;
    Info << "At time = " << ct << " flow error = " << flowController_->error(refFlow, actualFlow) <<  " new pressure will be = " << p0_ << endl;
    
    const scalarField p0
    (
        this->size(),
        p0_
    );
    
    updateCoeffs
    (
        p0,
        patch().lookupPatchField<volVectorField, vector>(UName())
    );
}


void Foam::flowRateControlledWithPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    os.writeKeyword("psi") << psiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << nl;
    os.writeKeyword("p0") << p0_ << token::END_STATEMENT << nl;
    
    os.incrIndent();
    os.incrIndent();
    flowRate_->writeData(os);
    os << endl;
    os << "flowRateController" <<  endl;
    flowController_->writeData(os);
    os.decrIndent();
    os.decrIndent();
    
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        flowRateControlledWithPressureFvPatchScalarField
    );
}

// ************************************************************************* //
