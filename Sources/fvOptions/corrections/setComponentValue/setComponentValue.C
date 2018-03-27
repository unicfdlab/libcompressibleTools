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

#include "setComponentValue.H"
#include "fvMesh.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(setComponentValue, 0);
    addToRunTimeSelectionTable
    (
        option,
        setComponentValue,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::setComponentValue::setComponentValue
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    direction_(readLabel(coeffs_.lookup("direction"))),
    dirValue_(readScalar(coeffs_.lookup("dirValue"))),
    vectorName_(coeffs_.lookup("vectorName"))
{
    // Set the field name to that of the energy field from which the temperature
    // is obtained

    //const psiThermo& thermo =
    //    mesh_.lookupObject<psiThermo>(psiThermo::dictName);
    
    fieldNames_.setSize(1, vectorName_);

    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::setComponentValue::correct(volVectorField& U)
{
    forAll(U, i)
    {
        U[i][direction_] = dirValue_;
    }
    
    // handle boundaries in the case of 'all'
    if (selectionMode_ == smAll)
    {
        volVectorField::Boundary& bf = U.boundaryFieldRef();
        
        forAll(bf, patchi)
        {
            fvPatchVectorField& Up = bf[patchi];
            
            forAll(Up, facei)
            {
                Up[facei][direction_] = dirValue_;
            }
        }
    }
}


bool Foam::fv::setComponentValue::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.readIfPresent("direction", direction_);
        coeffs_.readIfPresent("dirValue", dirValue_);
        coeffs_.readIfPresent("vectorName", vectorName_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
