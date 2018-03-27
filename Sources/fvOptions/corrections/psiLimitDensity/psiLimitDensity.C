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

#include "psiLimitDensity.H"
#include "fvMesh.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(psiLimitDensity, 0);
    addToRunTimeSelectionTable
    (
        option,
        psiLimitDensity,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::psiLimitDensity::psiLimitDensity
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    pmin_(readScalar(coeffs_.lookup("pmin"))),
    rhoMin_(readScalar(coeffs_.lookup("rhoMin"))),
    rhoName_(coeffs_.lookup("rho"))
{
    // Set the field name to that of the energy field from which the temperature
    // is obtained

    //const psiThermo& thermo =
    //    mesh_.lookupObject<psiThermo>(psiThermo::dictName);

    fieldNames_.setSize(1, rhoName_);

    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::psiLimitDensity::correct(volScalarField& rho)
{
    const fluidThermo& thermo =
        mesh_.lookupObject<fluidThermo>(fluidThermo::dictName);
    
    const volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();
    const volScalarField& rho0 = rho.oldTime();
    
    scalar pCell = pmin_;
    
    label cellI = -1;
    label rhol = 0.0;
    
    forAll(cells_, i)
    {
        cellI = cells_[i];
        if(rho[cellI] < rhoMin_)
        {
            pCell = max(pmin_, p[cellI]);
            rhol = psi[cellI] * pCell;
            rho[cellI] = max(rhol, rho0[cellI]);
        }
    }
    
    // handle boundaries in the case of 'all'
    if (selectionMode_ == smAll)
    {
        volScalarField::Boundary& bf = rho.boundaryFieldRef();
        scalar pFace = pmin_;
        
        forAll(bf, patchi)
        {
            fvPatchScalarField& rhop = bf[patchi];
            const scalarField& pp = p.boundaryField()[patchi];
            const scalarField& psip = psi.boundaryField()[patchi];
            const scalarField& rho0p = rho0.boundaryField()[patchi];
            
            forAll(rhop, facei)
            {
                if (rhop[facei] < rhoMin_)
                {
                    pFace = max(pmin_, pp[facei]);
                    rhol = pFace * psip[facei];
                    rhop[facei] = max(rhol, rho0p[facei]);
                }
            }
        }
    }
}


bool Foam::fv::psiLimitDensity::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.readIfPresent("pmin", pmin_);
        coeffs_.readIfPresent("rhoMin", rhoMin_);
        coeffs_.readIfPresent("rhoName", rhoName_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
