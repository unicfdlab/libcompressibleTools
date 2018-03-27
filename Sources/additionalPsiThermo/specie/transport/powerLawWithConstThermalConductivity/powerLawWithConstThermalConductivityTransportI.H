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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::powerLawWithConstThermalConductivityTransport<Thermo>::powerLawWithConstThermalConductivityTransport
(
    const Thermo& t,
    const scalar mu0,
    const scalar T0,
    const scalar kk,
    const scalar Pr,
    const scalar kappa0
)
:
    Thermo(t),
    mu0_(mu0),
    T0_(T0),
    k_(kk),
    rPr_(1.0/Pr),
    kappa0_(kappa0)
{}


template<class Thermo>
inline Foam::powerLawWithConstThermalConductivityTransport<Thermo>::powerLawWithConstThermalConductivityTransport
(
    const word& name,
    const powerLawWithConstThermalConductivityTransport& ct
)
:
    Thermo(name, ct),
    mu0_(ct.mu0_),
    T0_(ct.T0_),
    k_(ct.k_),
    rPr_(ct.rPr_),
    kappa0_(ct.kappa0_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::powerLawWithConstThermalConductivityTransport<Thermo> >
Foam::powerLawWithConstThermalConductivityTransport<Thermo>::clone() const
{
    return autoPtr<powerLawWithConstThermalConductivityTransport<Thermo> >
    (
        new powerLawWithConstThermalConductivityTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::powerLawWithConstThermalConductivityTransport<Thermo> >
Foam::powerLawWithConstThermalConductivityTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<powerLawWithConstThermalConductivityTransport<Thermo> >
    (
        new powerLawWithConstThermalConductivityTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::powerLawWithConstThermalConductivityTransport<Thermo> >
Foam::powerLawWithConstThermalConductivityTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<powerLawWithConstThermalConductivityTransport<Thermo> >
    (
        new powerLawWithConstThermalConductivityTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::powerLawWithConstThermalConductivityTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return mu0_ * Foam::pow(T / T0_, k_);
}


template<class Thermo>
inline Foam::scalar Foam::powerLawWithConstThermalConductivityTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return kappa0_;
}


template<class Thermo>
inline Foam::scalar Foam::powerLawWithConstThermalConductivityTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
    return kappa0_ / this->Cp(p,T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::powerLawWithConstThermalConductivityTransport<Thermo>& Foam::powerLawWithConstThermalConductivityTransport<Thermo>::operator=
(
    const powerLawWithConstThermalConductivityTransport<Thermo>& ct
)
{
    Thermo::operator=(ct);

    mu0_ = ct.mu0_;
    T0_ = ct.T0_;
    k_ = ct.k_;
    rPr_ = ct.rPr_;
    kappa0_ = ct.kappa0_;

    return *this;
}


template<class Thermo>
inline void Foam::powerLawWithConstThermalConductivityTransport<Thermo>::operator+=
(
    const powerLawWithConstThermalConductivityTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(st);

    if (mag(molr1) + mag(st.nMoles()) > SMALL)
    {
        molr1 /= this->nMoles();
        scalar molr2 = st.nMoles()/this->nMoles();

        rPr_ = 1.0/(molr1/rPr_ + molr2/st.rPr_);
    }
    
    notImplemented ("inline Foam::powerLawWithConstThermalConductivityTransport<Thermo> Foam::operator+=");
}


template<class Thermo>
inline void Foam::powerLawWithConstThermalConductivityTransport<Thermo>::operator-=
(
    const powerLawWithConstThermalConductivityTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(st);

    if (mag(molr1) + mag(st.nMoles()) > SMALL)
    {
        molr1 /= this->nMoles();
        scalar molr2 = st.nMoles()/this->nMoles();

        rPr_ = 1.0/(molr1/rPr_ - molr2/st.rPr_);
    }
    
    notImplemented ("inline Foam::powerLawWithConstThermalConductivityTransport<Thermo> Foam::operator-=");
}


template<class Thermo>
inline void Foam::powerLawWithConstThermalConductivityTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::powerLawWithConstThermalConductivityTransport<Thermo> Foam::operator+
(
    const powerLawWithConstThermalConductivityTransport<Thermo>& ct1,
    const powerLawWithConstThermalConductivityTransport<Thermo>& ct2
)
{
    notImplemented ("inline Foam::powerLawWithConstThermalConductivityTransport<Thermo> Foam::operator+");
}


template<class Thermo>
inline Foam::powerLawWithConstThermalConductivityTransport<Thermo> Foam::operator-
(
    const powerLawWithConstThermalConductivityTransport<Thermo>& ct1,
    const powerLawWithConstThermalConductivityTransport<Thermo>& ct2
)
{
    notImplemented ("inline Foam::powerLawWithConstThermalConductivityTransport<Thermo> Foam::operator-");
}


template<class Thermo>
inline Foam::powerLawWithConstThermalConductivityTransport<Thermo> Foam::operator*
(
    const scalar s,
    const powerLawWithConstThermalConductivityTransport<Thermo>& ct
)
{
    notImplemented ("inline Foam::powerLawWithConstThermalConductivityTransport<Thermo> Foam::operator*");
}


template<class Thermo>
inline Foam::powerLawWithConstThermalConductivityTransport<Thermo> Foam::operator==
(
    const powerLawWithConstThermalConductivityTransport<Thermo>& ct1,
    const powerLawWithConstThermalConductivityTransport<Thermo>& ct2
)
{
    notImplemented ("inline Foam::powerLawWithConstThermalConductivityTransport<Thermo> Foam::operator==");
}


// ************************************************************************* //
