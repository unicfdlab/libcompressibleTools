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

\*---------------------------------------------------------------------------*/

#include "varyingGammaTotalTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::varyingGammaTotalTemperatureFvPatchScalarField::varyingGammaTotalTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    T0_(p.size(), 0.0)
{}


Foam::varyingGammaTotalTemperatureFvPatchScalarField::varyingGammaTotalTemperatureFvPatchScalarField
(
    const varyingGammaTotalTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    T0_(ptf.T0_, mapper)
{}


Foam::varyingGammaTotalTemperatureFvPatchScalarField::varyingGammaTotalTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    T0_("T0", dict, p.size())
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
        fvPatchField<scalar>::operator=(T0_);
    }
}


Foam::varyingGammaTotalTemperatureFvPatchScalarField::varyingGammaTotalTemperatureFvPatchScalarField
(
    const varyingGammaTotalTemperatureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    T0_(tppsf.T0_)
{}


Foam::varyingGammaTotalTemperatureFvPatchScalarField::varyingGammaTotalTemperatureFvPatchScalarField
(
    const varyingGammaTotalTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    T0_(tppsf.T0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::varyingGammaTotalTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    T0_.autoMap(m);
}


void Foam::varyingGammaTotalTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const varyingGammaTotalTemperatureFvPatchScalarField& tiptf =
        refCast<const varyingGammaTotalTemperatureFvPatchScalarField>(ptf);

    T0_.rmap(tiptf.T0_, addr);
}


void Foam::varyingGammaTotalTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const fluidThermo& thermo = 
        this->patch().boundaryMesh().mesh().thisDb().lookupObject<fluidThermo>("thermophysicalProperties");
    
    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);
    
    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);
    
    const fvPatchField<scalar>& psip = 
        thermo.psi().boundaryField()[patch().index()];
    
    const fvPatchField<scalar>& pp =
        thermo.p().boundaryField()[patch().index()];
    
    const fvPatchField<scalar>& Tp =
        thermo.T().boundaryField()[patch().index()];
    
    scalarField gammap
    (
        thermo.gamma(pp, Tp, patch().index())
    );
    
    scalarField gM1ByG ((gammap - 1.0)/gammap);
    
    operator==
    (
        T0_/(1.0 + 0.5*psip*gM1ByG*(1.0 - pos(phip))*magSqr(Up))
    );
    
    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::varyingGammaTotalTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    T0_.writeEntry("T0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        varyingGammaTotalTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
