/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::powerLawWithConstThermalConductivityTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    powerLawWithConstThermalConductivityTransportI.H
    powerLawWithConstThermalConductivityTransport.C

\*---------------------------------------------------------------------------*/

#ifndef powerLawWithConstThermalConductivityTransport_H
#define powerLawWithConstThermalConductivityTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class powerLawWithConstThermalConductivityTransport;

template<class Thermo>
inline powerLawWithConstThermalConductivityTransport<Thermo> operator+
(
    const powerLawWithConstThermalConductivityTransport<Thermo>&,
    const powerLawWithConstThermalConductivityTransport<Thermo>&
);

template<class Thermo>
inline powerLawWithConstThermalConductivityTransport<Thermo> operator-
(
    const powerLawWithConstThermalConductivityTransport<Thermo>&,
    const powerLawWithConstThermalConductivityTransport<Thermo>&
);

template<class Thermo>
inline powerLawWithConstThermalConductivityTransport<Thermo> operator*
(
    const scalar,
    const powerLawWithConstThermalConductivityTransport<Thermo>&
);

template<class Thermo>
inline powerLawWithConstThermalConductivityTransport<Thermo> operator==
(
    const powerLawWithConstThermalConductivityTransport<Thermo>&,
    const powerLawWithConstThermalConductivityTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const powerLawWithConstThermalConductivityTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class powerLawWithConstThermalConductivityTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class powerLawWithConstThermalConductivityTransport
:
    public Thermo
{
    // Private data

        //- Constant dynamic viscosity [Pa.s]
        scalar mu0_;

        //- [K]
        scalar T0_;
        
        //- power coeff.
        scalar k_;

        //- Reciprocal Prandtl Number []
        scalar rPr_;
        
        //-
        scalar kappa0_;


    // Private Member Functions

        //- Construct from components
        inline powerLawWithConstThermalConductivityTransport
        (
            const Thermo& t,
            const scalar mu0,
            const scalar T0,
            const scalar kk,
            const scalar Pr,
            const scalar kappa0
        );


public:

    // Constructors

        //- Construct as named copy
        inline powerLawWithConstThermalConductivityTransport(const word&, const powerLawWithConstThermalConductivityTransport&);

        //- Construct from Istream
        powerLawWithConstThermalConductivityTransport(Istream&);

        //- Construct from dictionary
        powerLawWithConstThermalConductivityTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<powerLawWithConstThermalConductivityTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<powerLawWithConstThermalConductivityTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<powerLawWithConstThermalConductivityTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "powerLawWithConstThermalConductivity<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline powerLawWithConstThermalConductivityTransport& operator=(const powerLawWithConstThermalConductivityTransport&);

        inline void operator+=(const powerLawWithConstThermalConductivityTransport&);

        inline void operator-=(const powerLawWithConstThermalConductivityTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend powerLawWithConstThermalConductivityTransport operator+ <Thermo>
        (
            const powerLawWithConstThermalConductivityTransport&,
            const powerLawWithConstThermalConductivityTransport&
        );

        friend powerLawWithConstThermalConductivityTransport operator- <Thermo>
        (
            const powerLawWithConstThermalConductivityTransport&,
            const powerLawWithConstThermalConductivityTransport&
        );

        friend powerLawWithConstThermalConductivityTransport operator* <Thermo>
        (
            const scalar,
            const powerLawWithConstThermalConductivityTransport&
        );

        friend powerLawWithConstThermalConductivityTransport operator== <Thermo>
        (
            const powerLawWithConstThermalConductivityTransport&,
            const powerLawWithConstThermalConductivityTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const powerLawWithConstThermalConductivityTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "powerLawWithConstThermalConductivityTransportI.H"

#ifdef NoRepository
#   include "powerLawWithConstThermalConductivityTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
