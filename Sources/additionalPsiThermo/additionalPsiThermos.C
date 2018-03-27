/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "psiThermo.H"
#include "makeThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "PengRobinsonGas.H"
#include "hConstThermo.H"
#include "eConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "hPolynomialThermo.H"
#include "polynomialTransport.H"

#include "hePsiThermo.H"
#include "pureMixture.H"

#include "specie/transport/powerLaw/powerLawTransport.H"
#include "specie/transport/powerLawWithConstThermalConductivity/powerLawWithConstThermalConductivityTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * * Enthalpy-based * * * * * * * * * * * * * */

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    powerLawTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    powerLawWithConstThermalConductivityTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    powerLawTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    perfectGas,
    specie
);

//makeThermo
//(
//    psiThermo,
//    hePsiThermo,
//    pureMixture,
//    polynomial,
//    sensibleEnthalpy,
//    hPolynomialThermo,
//    perfectGas,
//    specie
//);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
