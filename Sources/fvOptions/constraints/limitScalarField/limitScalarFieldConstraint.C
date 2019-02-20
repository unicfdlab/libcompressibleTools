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

\*---------------------------------------------------------------------------*/

#include "limitScalarFieldConstraint.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(limitScalarFieldConstraint, 0);
        addToRunTimeSelectionTable
        (
            option,
            limitScalarFieldConstraint,
            dictionary
        );
    }
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitScalarFieldConstraint::limitScalarFieldConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    fieldValue_(nullptr),
    scaleCoeff_(1.0),
    fieldName_("p")
{

    this->read(dict);
    
    //can be applied only to one field
    fieldNames_.setSize (1, fieldName_);
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::limitScalarFieldConstraint::constrain
(
    fvMatrix<scalar>& eqn,
    const label
)
{
    const volScalarField& field = 
        mesh().lookupObject<volScalarField>(fieldName_);

    scalar minValue = fieldValue_->value(mesh().time().value());
    cells_.resize(0);
    
    forAll(field, cellI)
    {
        if (field[cellI]*scaleCoeff_  <= minValue)
        {
            cells_.append(cellI);
        }
    }
    
    scalarField minValues (cells_.size(), minValue);
    
    eqn.setValues(cells_, minValues);
}


bool Foam::fv::limitScalarFieldConstraint::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        fieldValue_.reset
        (
            Function1<scalar>::New("minFieldValue", coeffs_).ptr()
        );
        
        coeffs_.lookup("fieldName") >> fieldName_;
        
        coeffs_.lookup("scaleCoeff") >> scaleCoeff_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
