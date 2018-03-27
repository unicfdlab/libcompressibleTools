/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "maxCellCourant.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(maxCellCourant, 0);
    
    addToRunTimeSelectionTable(functionObject, maxCellCourant, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::maxCellCourant::maxCellCourant
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject
    (
        name,
        runTime,
        dict
    ),
    CourantType_(word::null),
    rhoName_("rho"),
    phiName_("phi")
{
    this->read(dict);
}


Foam::functionObjects::maxCellCourant::maxCellCourant
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    regionFunctionObject
    (
        name,
        obr,
        dict
    ),
    CourantType_(word::null),
    rhoName_("rho"),
    phiName_("phi")
{
    this->read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::maxCellCourant::~maxCellCourant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::maxCellCourant::read(const dictionary& dict)
{
    if(regionFunctionObject::read(dict))
    {
        dict.lookup("CourantType") >> CourantType_;
        
        if (dict.found("rho"))
        {
            dict.lookup("rho") >> rhoName_;
        }
        
        if (dict.found("phi"))
        {
            dict.lookup("phi") >> phiName_;
        }
    }
    
    return false;
}


bool Foam::functionObjects::maxCellCourant::write()
{
    if (time_.outputTime())
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);
        const Time& runTime = mesh.time();
        
        autoPtr<volScalarField> localCoPtr;
        localCoPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "Co",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("localCo", dimless,  1.0)
            )
        );
        
        const volScalarField& rho = mesh.thisDb().lookupObject<volScalarField>(rhoName_);
        const surfaceScalarField& phi = mesh.thisDb().lookupObject<surfaceScalarField>(phiName_);
        
        if (CourantType_ == "cellCourant")
        {
            localCoPtr().primitiveFieldRef() = 
            (
                fvc::surfaceSum(mag(phi))().primitiveField()
                /
                rho.primitiveField()
                /
                mesh.V().field()
            ) * 0.5 * runTime.deltaT().value();
        }
        else if ( CourantType_ == "faceCourant" )
        {
            const surfaceScalarField& rho_own = mesh.thisDb().lookupObject<surfaceScalarField>("rho_own");
            const surfaceScalarField& rho_nei = mesh.thisDb().lookupObject<surfaceScalarField>("rho_nei");
            const surfaceScalarField& alpha_own = mesh.thisDb().lookupObject<surfaceScalarField>("alpha_own");
            const surfaceScalarField& alpha_nei = mesh.thisDb().lookupObject<surfaceScalarField>("alpha_nei");
            const surfaceScalarField& uMagSf    = mesh.thisDb().lookupObject<surfaceScalarField>("uMagSf");
            
            surfaceScalarField magUf =
            mag
            (
                phi
                / 
                (
                    rho_own * alpha_own
                    +
                    rho_nei * alpha_nei
                ) / uMagSf
            );
            
            surfaceScalarField faceCo =
                mesh.surfaceInterpolation::deltaCoeffs()
                *
                magUf
                *
                runTime.deltaT() * 2.0;
            
            forAll(mesh.cells(), iCell)
            {
                scalar maxCoCell =  -0.01;
                
                const labelList& cellFaces = mesh.cells()[iCell];
                
                forAll(cellFaces, iFace)
                {
                    if (mesh.isInternalFace(cellFaces[iFace]))
                    {
                        if (faceCo[cellFaces[iFace]] > maxCoCell)
                        {
                            maxCoCell = faceCo[cellFaces[iFace]];
                        }
                   }
                }
                
                localCoPtr()[iCell] = maxCoCell;
            }
            
            forAll(localCoPtr().boundaryField(), iPatch)
            {
                forAll(localCoPtr().boundaryField()[iPatch], iFace)
                {
                    localCoPtr().boundaryFieldRef()[iPatch][iFace] = 
                        faceCo.boundaryField()[iPatch][iFace];
                }
            }
        }
        else
        {
            FatalErrorIn
            (
                "maxCellCourant.C:"
            )   << "Wrong type of Courant criterion: " << CourantType_
            << endl << " must be one of:" 
            << endl << "1) cellCourant"
            << endl << "2) faceCourant"
            << endl << abort(FatalError);
        }
        
        if (localCoPtr.valid())
        {
            localCoPtr().write();
        }
        return true;
    }
    return true;
}

bool Foam::functionObjects::maxCellCourant::execute()
{
    return true;
}





// ************************************************************************* //
