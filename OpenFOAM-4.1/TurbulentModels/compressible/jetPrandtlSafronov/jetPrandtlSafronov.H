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
    Foam::RASModels::jetPrandtlSafronov

Group
    grpRASTurbulence

Description
    Prandtl-Safronov turbulent viscosity model for supersonic turbulent jets
    
    The default model coefficients are
    \verbatim
        jetPrandtlSafronovCoeffs
        {
            Cg              1.0;
        }
    \endverbatim

See also
    Foam::RASModels::kEpsilon

SourceFiles
    jetPrandtlSafronov.C

\*---------------------------------------------------------------------------*/

#ifndef jetPrandtlSafronov_H
#define jetPrandtlSafronov_H

#include "kEpsilon.H"
#include "meshSearch.H"
#include "faceOnlySet.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class jetPrandtlSafronov Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class jetPrandtlSafronov
:
    public kEpsilon<BasicTurbulenceModel>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        jetPrandtlSafronov(const jetPrandtlSafronov&);
        void operator=(const jetPrandtlSafronov&);


protected:

    // Protected data

        // Model coefficients

        word GasName_;
        word pName_;
        word UName_;
        word TName_;
        word hName_;
        vector axialDir_;
        vector radialDir_;
        word inletPatchName_;
        scalar pext_;
        scalar Text_;
        scalar Cpext_;
        scalar DIn_;
        meshSearch searchEngine_;
        point axisStart_;
        point axisEnd_;
        autoPtr<faceOnlySet> axisSamplePtr_;

    // Protected Member Functions

        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;
        


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("jetPrandtlSafronov");


    // Constructors

        //- Construct from components
        jetPrandtlSafronov
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
    virtual ~jetPrandtlSafronov()
    {}


    // Member Functions

        //- Solve the turbulence equations and correct the turbulence viscosity
        void correct();

        //- Re-read model coefficients if they have changed
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "jetPrandtlSafronov.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
