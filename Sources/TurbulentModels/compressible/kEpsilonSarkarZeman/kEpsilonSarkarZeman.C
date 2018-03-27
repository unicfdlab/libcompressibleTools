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

#include "kEpsilonSarkarZeman.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kEpsilonSarkarZeman<BasicTurbulenceModel>::kEpsilonSarkarZeman
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

    alpha1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alpha1",
            this->coeffDict_,
            2.5
        )
    ),

    alpha2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alpha2",
            this->coeffDict_,
            2.0
        )
    ),

    alpha3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alpha3",
            this->coeffDict_,
            0.0
        )
    ),

    lambdaT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "lambdaT",
            this->coeffDict_,
            0.2
        )
    ),

    etab_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "etab",
            this->coeffDict_,
            0
        )
    ),

    Mt_
    (
        (this->transport_.Cp() / this->transport_.Cv())*0.0 //init as zero w/o dimensions
    ),
    
    phiNames_(NULL)
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
    
    this->read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kEpsilonSarkarZeman<BasicTurbulenceModel>::read()
{
    if (kEpsilon<BasicTurbulenceModel>::read())
    {
        alpha1_.readIfPresent(this->coeffDict());
        alpha2_.readIfPresent(this->coeffDict());
        alpha3_.readIfPresent(this->coeffDict());
        lambdaT_.readIfPresent(this->coeffDict());
        etab_.readIfPresent(this->coeffDict());

        if (this->coeffDict().found("phiNames") && phiNames_.empty())
        {
            phiNames_.set(new Pair<word>(word::null, word::null));
            this->coeffDict().lookup("phiNames") >> phiNames_();
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
kEpsilonSarkarZeman<BasicTurbulenceModel>::kSource() const
{
    tmp<fvScalarMatrix> tkSp
    (
        kEpsilon<BasicTurbulenceModel>::kSource()
    );
    
    return tkSp;
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
kEpsilonSarkarZeman<BasicTurbulenceModel>::epsilonSource() const
{
    tmp<fvScalarMatrix> tepsilonSp
    (
       kEpsilon<BasicTurbulenceModel>::epsilonSource()
    );
    
    return tepsilonSp;
}

template<class BasicTurbulenceModel>
void
kEpsilonSarkarZeman<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    //update Sarkar-Zeman fields
    Mt_ =  sqrt(2.0*this->k_)/ sqrt(this->transport_.Cp()/this->transport_.Cv()/this->transport_.psi());
    Mt_ = max(Mt_ - lambdaT_, dimensionedScalar("zero", dimless, 0));
    volScalarField MtSqr = Mt_*Mt_;

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));
    volScalarField& epsilon = this->epsilon_;
    volScalarField& k = this->k_;
    
    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();
    
    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );
    
    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );
    tgradU.clear();
    
    // Update epsilon and G at the wall
    epsilon.boundaryFieldRef().updateCoeffs();
    
    autoPtr<fvScalarMatrix> epsConv;
    autoPtr<fvScalarMatrix> kConv;
    
    if (phiNames_.valid())
    {
        const surfaceScalarField& phiOwn = epsilon.mesh().objectRegistry::template
            lookupObject<surfaceScalarField>(phiNames_().first());

        const surfaceScalarField& phiNei = epsilon.mesh().objectRegistry::template
            lookupObject<surfaceScalarField>(phiNames_().second());
        
        epsConv.set
        (
            new fvScalarMatrix
            (
                fvm::div(phiOwn, epsilon)
                +
                fvm::div(phiNei, epsilon)
            )
        );

        kConv.set
        (
            new fvScalarMatrix
            (
                fvm::div(phiOwn, k)
                +
                fvm::div(phiNei, k)
            )
        );
    }
    else
    {
        epsConv.set
        (
            new fvScalarMatrix
            (
                fvm::div(alphaRhoPhi, epsilon)
            )
        );

        kConv.set
        (
            new fvScalarMatrix
            (
                fvm::div(alphaRhoPhi, k)
            )
        );
    }
    
    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon)
        + epsConv()
        - fvm::laplacian(alpha*rho*this->DepsilonEff(), epsilon)
        ==
        this->C1_*alpha()*rho()*G*epsilon()/k()
        - fvm::SuSp(((2.0/3.0)*this->C1_ + this->C3_)*alpha()*rho()*divU, epsilon)
        - fvm::Sp(this->C2_*alpha()*rho()*epsilon()/k(), epsilon)
        + epsilonSource()
        + fvOptions(alpha, rho, epsilon)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon);
    bound(epsilon, this->epsilonMin_);
    
    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k)
        + kConv()
        - fvm::laplacian(alpha*rho*this->DkEff(), k)
        ==
        alpha()*rho()*G
        - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k)
        - fvm::Sp(alpha()*rho()*epsilon()/k(), k)
        + rho()*(-this->alpha1_*MtSqr())*G
        - rho()*( this->alpha2_*MtSqr())*epsilon()
        + kSource()
        + fvOptions(alpha, rho, k)
    );
    
    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k);
    bound(k, this->kMin_);
    
    correctNut();

}

template<class BasicTurbulenceModel>
void
kEpsilonSarkarZeman<BasicTurbulenceModel>::correctNut()
{

    volScalarField MtSqr = Mt_*Mt_;

    this->nut_ = this->Cmu_*(1.0/(1.0 + alpha3_*MtSqr))*sqr(this->k_)/this->epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
    
    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<fvVectorMatrix>
kEpsilonSarkarZeman<BasicTurbulenceModel>::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    
    volTensorField gradUT = T(fvc::grad(U));
    
    dimensionedScalar etab
    (
        "etab",
        dimMass/dimLength/dimTime,
        this->etab_.value()
    );
    
    return
    (
        //- fvc::grad(etab*tr(gradUT))
        - fvc::div((this->alpha_*rho*this->nuEff())*dev2(gradUT))
        - fvm::laplacian(this->alpha_*rho*this->nuEff(), U)
    );
                        
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
