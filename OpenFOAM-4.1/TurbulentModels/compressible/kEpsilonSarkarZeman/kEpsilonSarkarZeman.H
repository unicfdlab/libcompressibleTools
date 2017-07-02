/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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

Class
    Foam::RASModels::kEpsilonSarkarZeman

Group
    grpRASTurbulence

Description
    Sarkar-Zeman modification of k-epsilon model to account for compressible
    flows
    
    The default model coefficients are
    \verbatim
        kEpsilonSarkarZemanCoeffs
        {
            Cg              1.0;
        }
    \endverbatim

See also
    Foam::RASModels::kEpsilon

SourceFiles
    kEpsilonSarkarZeman.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilonSarkarZeman_H
#define kEpsilonSarkarZeman_H

#include "kEpsilon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kEpsilonSarkarZeman Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kEpsilonSarkarZeman
:
    public kEpsilon<BasicTurbulenceModel>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        kEpsilonSarkarZeman(const kEpsilonSarkarZeman&);
        void operator=(const kEpsilonSarkarZeman&);


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar alpha1_;
            dimensionedScalar alpha2_;
            dimensionedScalar alpha3_;
            dimensionedScalar lambdaT_;
            dimensionedScalar etab_;
        
        // Model fields
        
            volScalarField Mt_;
            
        // Characteristic fluxes
        autoPtr<Pair<word> > phiNames_;

    // Protected Member Functions

        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;
        


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kEpsilonSarkarZeman");


    // Constructors

        //- Construct from components
        kEpsilonSarkarZeman
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kEpsilonSarkarZeman()
    {}


    // Member Functions

        //- Solve the turbulence equations and correct the turbulence viscosity
        void correct();

        //- Re-read model coefficients if they have changed
        bool read();
        
        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kEpsilonSarkarZeman.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
