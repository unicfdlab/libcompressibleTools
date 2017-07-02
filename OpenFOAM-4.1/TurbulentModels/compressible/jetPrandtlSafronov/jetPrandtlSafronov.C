/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 OpenFOAM Foundation
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

#include "jetPrandtlSafronov.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
//#include "volPointInterpolation.H"
//#include "isoSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
jetPrandtlSafronov<BasicTurbulenceModel>::jetPrandtlSafronov
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kEpsilon<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    GasName_
    (
        "Gas"
    ),
    pName_
    (
        "p"
    ),
    UName_
    (
        "U"
    ),
    TName_
    (
        "T"
    ),
    hName_
    (
        "h"
    ),
    axialDir_(1, 0, 0),
    radialDir_(0, 1, 0),
    inletPatchName_("inlet"),
    pext_(101325),
    Text_(292.15),
    Cpext_(1005.0),
    DIn_(0.286),
    searchEngine_(U.mesh()),
    axisStart_(0, 0.001, 0),
    axisEnd_(15, 0.001, 0),
    axisSamplePtr_
    (
        new faceOnlySet
        (
            "nozzleAxis",
            U.mesh(),
            searchEngine_,
            "xyz",
            axisStart_,
            axisEnd_
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
    
    this->read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool jetPrandtlSafronov<BasicTurbulenceModel>::read()
{
    if (kEpsilon<BasicTurbulenceModel>::read())
    {
        if (this->coeffDict().found("Gas"))
        {
            this->coeffDict().lookup("Gas") >> GasName_;
        }

        if (this->coeffDict().found("p"))
        {
            this->coeffDict().lookup("p") >> pName_;
        }

        if (this->coeffDict().found("U"))
        {
            this->coeffDict().lookup("U") >> UName_;
        }

        if (this->coeffDict().found("T"))
        {
            this->coeffDict().lookup("T") >> TName_;
        }
        
        if (this->coeffDict().found("h"))
        {
            this->coeffDict().lookup("h") >> hName_;
        }

        if (this->coeffDict().found("axialDir"))
        {
            this->coeffDict().lookup("axialDir") >> axialDir_;
        }
        
        if (this->coeffDict().found("radialDir"))
        {
            this->coeffDict().lookup("radialDir") >> radialDir_;
        }
        
        if (this->coeffDict().found("inlet"))
        {
            this->coeffDict().lookup("inlet") >> inletPatchName_;
        }
        
        if (this->coeffDict().found("pext"))
        {
            this->coeffDict().lookup("pext") >> pext_;
        }

        if (this->coeffDict().found("Text"))
        {
            this->coeffDict().lookup("Text") >> Text_;
        }

        if (this->coeffDict().found("Cpext"))
        {
            this->coeffDict().lookup("Cpext") >> Cpext_;
        }

        if (this->coeffDict().found("DIn"))
        {
            this->coeffDict().lookup("DIn") >> DIn_;
        }
        
        if (this->coeffDict().found("axisStart"))
        {
            if (this->coeffDict().found("axisEnd"))
            {
                this->coeffDict().lookup("axisStart") >> axisStart_;
                this->coeffDict().lookup("axisEnd") >> axisEnd_;

                this->axisSamplePtr_.reset
                (
                    new faceOnlySet
                    (
                        "nozzleAxis",
                        this->epsilon_.mesh(),
                        searchEngine_,
                        "xyz",
                        axisStart_,
                        axisEnd_
                    )
                );
            }
        }

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
jetPrandtlSafronov<BasicTurbulenceModel>::kSource() const
{
    tmp<fvScalarMatrix> tkSp
    (
        kEpsilon<BasicTurbulenceModel>::kSource()
    );
    
    return tkSp;
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
jetPrandtlSafronov<BasicTurbulenceModel>::epsilonSource() const
{
    tmp<fvScalarMatrix> tepsilonSp
    (
       kEpsilon<BasicTurbulenceModel>::epsilonSource()
    );
    
    return tepsilonSp;
}

template<class BasicTurbulenceModel>
void
jetPrandtlSafronov<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    
    const fvMesh& mesh = this->epsilon_.mesh();
    
    const label inPatchID = mesh.boundaryMesh().findPatchID(inletPatchName_);
    
    if (inPatchID < 0)
    {
        Info << "Cant find inlet patch to estimate kinematic turbulent viscosity" << endl;
    }
    
    const volScalarField& YGas = mesh.objectRegistry::template
            lookupObject<volScalarField>(GasName_);
            
    const volScalarField& p = mesh.objectRegistry::template
            lookupObject<volScalarField>(pName_);
            
    const volVectorField& U = mesh.objectRegistry::template
            lookupObject<volVectorField>(UName_);

    const volScalarField& T = mesh.objectRegistry::template
            lookupObject<volScalarField>(TName_);

    const volScalarField& h = mesh.objectRegistry::template
            lookupObject<volScalarField>(hName_);
    
    volScalarField Cp ("Cp", this->transport_.Cp());
    volScalarField gamma ("gamma",  Cp / this->transport_.Cv());
    
    volScalarField c ("c", sqrt(gamma/this->transport_.psi()));
    
    volScalarField M ("M", mag(U)/c);
    
    scalar inPatchSumSf = gSum(mesh.magSf().boundaryField()[inPatchID]);
    
    Info << "inlet area = " << inPatchSumSf << endl;
    
    vector UIn = gSum(U.boundaryField()[inPatchID] * mesh.magSf().boundaryField()[inPatchID]) / 
        inPatchSumSf;
    
    scalar gammaIn = gSum(gamma.boundaryField()[inPatchID] * mesh.magSf().boundaryField()[inPatchID]) / 
        inPatchSumSf;
    
    scalar pIn = gSum(p.boundaryField()[inPatchID] * mesh.magSf().boundaryField()[inPatchID]) / 
        inPatchSumSf;
    
    scalar MIn = gSum(M.boundaryField()[inPatchID] * mesh.magSf().boundaryField()[inPatchID]) / 
        inPatchSumSf;
        
    scalar cIn = gSum(c.boundaryField()[inPatchID] * mesh.magSf().boundaryField()[inPatchID]) / 
        inPatchSumSf;
        
//    scalar TIn = gSum(T.boundaryField()[inPatchID] * mesh.magSf().boundaryField()[inPatchID]) / 
//        inPatchSumSf;
//        
//    scalar CpIn = gSum(Cp.boundaryField()[inPatchID] * mesh.magSf().boundaryField()[inPatchID]) / 
//        inPatchSumSf;

    scalar hIn = gSum(h.boundaryField()[inPatchID] * mesh.magSf().boundaryField()[inPatchID]) / 
        inPatchSumSf;
    
    Info << "gammaIn  = " << gammaIn << endl;
    Info << "MIn      = " << MIn << endl;
    
    scalar MInEquiv = sqrt
            (
                (2.0 / (gammaIn - 1.0))
                *
                (
                    (1.0 + (gammaIn - 1.0)*0.5*MIn*MIn)
                    *
                    (pow(pIn/pext_, (gammaIn-1.0)/gammaIn))
                    -
                    1
                )
            );
    
    Info << "MInEquiv = " << MInEquiv << endl;
    
    scalar UInEquiv = MInEquiv * cIn * sqrt (pow(pIn/pext_, (gammaIn-1.0)/gammaIn));
    
    scalar a1 = (gammaIn - 1.0) / 2.0;
    scalar a2 = (gammaIn + 1.0) / (gammaIn - 1.0) / 2.0;
    
    scalar DInEquiv = DIn_ * sqrt ( (MIn / MInEquiv) * pow(
        (1.0 + a1*MIn*MIn)/(1.0 + a1*MInEquiv*MInEquiv),
        a2));
    
    Info << "Equiv variables calculated" << endl;
    
    scalar epsH = Cpext_*Text_ / (hIn + 0.5*magSqr(UIn));
    
    scalar nutMaxRel = 0.015*pow(epsH, 0.2);
    
    scalar nut1Rel = 0.0055 + 0.00079*log(1.352*epsH); //log is ln in OpenFOAM
    
    scalar xT = DInEquiv*pow(0.22+4.985*gammaIn*MInEquiv*MInEquiv,0.45)*(0.64+0.36*epsH);
    
    volScalarField nutRel ("nutRel", gamma*0.0); //set equal to 0
    
    forAll(nutRel, iCell)
    {
        nutRel.primitiveFieldRef()[iCell] = min(nut1Rel*(mesh.C()[iCell] & axialDir_) / xT, nutMaxRel);
    }
    forAll(nutRel.boundaryField(), iPatch)
    {
        forAll(nutRel.boundaryField()[iPatch], iFace)
        {
            nutRel.boundaryFieldRef()[iPatch][iFace] = 
                min(nut1Rel*(mesh.Cf().boundaryField()[iPatch][iFace] & axialDir_) / xT, nutMaxRel);
        }
    }
    
    //serial version
    interpolationTable<scalar> gasOnAxis;
    gasOnAxis.outOfBounds(interpolationTable<scalar>::CLAMP);
    labelList cellsOnAxis (axisSamplePtr_().cells());
    label cellI = -1;
    forAll(cellsOnAxis, iCell)
    {
        cellI = cellsOnAxis[iCell];
        gasOnAxis.append
        (
            Tuple2<scalar,scalar>
            (
                mesh.C()[cellI] & axialDir_,
                YGas[cellI]
            )
        );
        
        //Info << "cellI = " << gasOnAxis[iCell] << endl;
    }
    
    volScalarField theta ("theta", nutRel*0.0);
    scalar r0Gas = 0.0;
    forAll(theta, iCell)
    {
        r0Gas = gasOnAxis(mesh.C()[iCell] & axialDir_);
        if (r0Gas > 1.0e-4)
        {
            theta.primitiveFieldRef()[iCell] =
            min
            (
                1.0,
                0.1*
                YGas[iCell]/r0Gas
            );
        }
        else
        {
            theta.primitiveFieldRef()[iCell] = 0.0;
        }
    }
    forAll(theta.boundaryField(), iPatch)
    {
        forAll(theta.boundaryField()[iPatch], iFace)
        {
            r0Gas = gasOnAxis(mesh.Cf().boundaryField()[iPatch][iFace] & axialDir_);
            if (r0Gas > 1.0e-4)
            {
                theta.boundaryFieldRef()[iPatch][iFace] = min
                (
                    1.0, 
                    0.1*YGas.boundaryField()[iPatch][iFace] / r0Gas
                );
            }
            else
            {
                theta.boundaryFieldRef()[iPatch][iFace] = 0.0;
            }
        }
    }
    
    dimensionedScalar Ua ("Ua", dimLength / dimTime, UInEquiv);
    dimensionedScalar Da ("Da", dimLength, DInEquiv);
    dimensionedScalar zeronut ("zeronut", dimLength*dimLength / dimTime, 0.0);
    
    Info << "Ua = " << Ua.value() << endl;
    Info << "Da = " << Da.value() << endl;
    
    //corect kinematic viscosity
    this->nut_ = max(nutRel*theta*Ua*Da - this->transport_.mu() / this->transport_.rho(), zeronut);
    
    if (this->nut_.time().outputTime())
    {
        nutRel.write();
        theta.write();
    }
    
    correctNut();

    Info<< "Max/min of turbulent kinematic viscosity: "
        << max(this->nut_()).value()
        << "/"
        << min(this->nut_()).value()
        << endl;
}

template<class BasicTurbulenceModel>
void
jetPrandtlSafronov<BasicTurbulenceModel>::correctNut()
{

//    volScalarField MtSqr = Mt_*Mt_;
//
//    this->nut_ = this->Cmu_*(1.0/(1.0 + alpha3_*MtSqr))*sqr(this->k_)/this->epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
    
    BasicTurbulenceModel::correctNut();
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
