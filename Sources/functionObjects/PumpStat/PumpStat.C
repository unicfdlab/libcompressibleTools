/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "PumpStat.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "wordReList.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(PumpStat, 0);
    
    addToRunTimeSelectionTable(functionObject, PumpStat, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::PumpStat::PumpStat
(
    const word& name,
    const Time& time,
    const dictionary& dict
)
:
    forces(name, time, dict),
    phiName_("phi"),
    inflowPatches_(),
    outflowPatches_(),
    timeStart_(-1.0),
    timeEnd_(-1.0),
    omega_(vector::zero),
    values_(0),
    vnames_(0)
{
    this->read(dict);
    this->initialize();
}

Foam::functionObjects::PumpStat::PumpStat
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    forces(name, obr, dict),
    phiName_("phi"),
    inflowPatches_(),
    outflowPatches_(),
    timeStart_(-1.0),
    timeEnd_(-1.0),
    omega_(vector::zero),
    values_(0),
    vnames_(0)
{
    this->read(dict);
    this->initialize();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::PumpStat::~PumpStat()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::PumpStat::initialize()
{
    vnames_.resize(8);
    
    vnames_[0] = "Ntotal";
    vnames_[1] = "Mx";
    vnames_[2] = "My";
    vnames_[3] = "Mz";
    vnames_[4] = "eta";
    vnames_[5] = "pratio";
    vnames_[6] = "flux";
    vnames_[7] = "Tout";
    
    values_.resize(8, 0.0);
}

bool Foam::functionObjects::PumpStat::read(const dictionary& dict)
{
    if (!forces::read(dict))
    {
        return false;
    }
    
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        
    outflowPatches_ = pbm.patchSet(wordReList(dict.lookup("outletPatches")));
    inflowPatches_  = pbm.patchSet(wordReList(dict.lookup("inletPatches")));

    if (dict.found("phi"))
    {
        dict.lookup("phi") >> phiName_;
    }

    dict.lookup("timeStart") >> timeStart_;
    
    dict.lookup("timeEnd") >> timeEnd_;
    
    dict.lookup("omega") >> omega_;
    
    return true;
}

Foam::tmp<Foam::volSymmTensorField> Foam::functionObjects::PumpStat::totalStress()
{
    tmp<volSymmTensorField> tStress
    (
        devRhoReff()
    );
    
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);
    
    if (p.dimensions() != dimPressure)
    {
        if (rhoName_ == "rhoInf")
        {
            tStress.ref() += p*dimensionedScalar("rhoRef",dimDensity,rhoRef_)*symmTensor::I;
        }
        else
        {
            const volScalarField& rho = obr_.lookupObject<volScalarField>(rhoName_);
            tStress.ref() += p*rho*symmTensor::I;
        }
    }
    else
    {
        tStress.ref() += p*symmTensor::I;
    }
    
    return tStress;
}

Foam::scalar Foam::functionObjects::PumpStat::patchFlow(const word& phiName, const labelHashSet& patches)
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    const surfaceScalarField& phi  = mesh.lookupObject<surfaceScalarField>(phiName);
    
    scalar totFlux = 0.0;

    forAllConstIter(labelHashSet, patches, iter)
    {
        label patchi = iter.key();
        totFlux += gSum(phi.boundaryField()[patchi]);
    }
        
    return totFlux;
}

Foam::scalar Foam::functionObjects::PumpStat::patchH(const labelHashSet& patches)
{
    const volScalarField& he = obr_.lookupObject<volScalarField>("h");
    
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    
    scalar Havr = 0.0;
    scalar mSf  = 0.0;

    forAllConstIter(labelHashSet, patches, iter)
    {
        label patchi = iter.key();
        Havr += gSum(he.boundaryField()[patchi] * mesh.magSf().boundaryField()[patchi]);
        mSf  += gSum(mesh.magSf().boundaryField()[patchi]);
    }
    
    if (mSf > VSMALL)
    {
        Havr /= mSf;
    }
    
    return Havr;
}

Foam::scalar Foam::functionObjects::PumpStat::patchP(const labelHashSet& patches)
{
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);
    
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    
    scalar Pavr = 0.0;
    scalar mSf  = 0.0;
    
    forAllConstIter(labelHashSet, patches, iter)
    {
        label patchi = iter.key();
        Pavr += gSum(p.boundaryField()[patchi] * mesh.magSf().boundaryField()[patchi]);
        mSf  += gSum(mesh.magSf().boundaryField()[patchi]);
    }
    
    if (mSf > VSMALL)
    {
        Pavr /= mSf;
    }
    
    return Pavr;
}

void Foam::functionObjects::PumpStat::correct()
{
    volSymmTensorField Pi(totalStress());
    
    vector rotAxis (omega_ / mag(omega_));
    vector Mtot(vector::zero);
    scalar Ntot(0.0);
    scalar outFlow(0.0);
    scalar inFlow(0.0);
    scalar outH(0.0);
    scalar inH(0.0);
    scalar cWork(0.0);
    scalar outP(0.0);
    scalar inP(0.0);
    scalar pratio(0.0);
    scalar eta(0.0);
    
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        
        vectorField Sf = Pi.mesh().Sf().boundaryField()[patchi];
        vectorField F = Pi.boundaryField()[patchi] & Sf;
        vectorField pCf = Pi.mesh().Cf().boundaryField()[patchi];
        vectorField B = rotAxis*(rotAxis & pCf);
        vectorField arm = pCf - B;
        vectorField moments= arm ^ F;
        Mtot += gSum(moments);
    }
    
    Ntot = Mtot & omega_;
    
    outFlow = patchFlow(phiName_, outflowPatches_);
    inFlow  = patchFlow(phiName_, inflowPatches_);
    outH    = patchH(outflowPatches_);
    inH     = patchH(inflowPatches_);
    cWork   = outFlow*(outH - inH);
    if (mag(Ntot) > VSMALL)
    {
        eta     = cWork / Ntot;
    }
    
    outP    = patchP(outflowPatches_);
    inP     = patchP(inflowPatches_);
    if (inP > VSMALL)
    {
        pratio = outP / inP;
    }
    
    values_[0] = Ntot;
    values_[1] = Mtot.x();
    values_[2] = Mtot.y();
    values_[3] = Mtot.z();
    values_[4] = eta;
    values_[5] = pratio;
    values_[6] = outFlow;
//        values_[7].append (Tout);

//
//    //calculate moments
//    vector Mtotal = -gatherMoments();
//    reduce (Mtotal, sumOp<vector>());
//    
//    //calculate power and efficiency
//    scalar Ntotal = 0.0;
//    scalar inArea = patchesArea(inflowPatches_); reduce (inArea, sumOp<scalar>());
//    scalar outArea= patchesArea(outflowPatches_); reduce(outArea, sumOp<scalar>());
//    scalar Tin    = volPatchIntegrate(TName_, inflowPatches_); reduce(Tin, sumOp<scalar>());
//    scalar Tout   = volPatchIntegrate(TName_, outflowPatches_); reduce(Tout, sumOp<scalar>());
//    scalar pin    = volPatchIntegrate(pName_, inflowPatches_); reduce(pin, sumOp<scalar>());
//    scalar pout   = volPatchIntegrate(pName_, outflowPatches_); reduce(pout, sumOp<scalar>());
//    scalar rhoIn  = volPatchIntegrate(rhoName_, inflowPatches_); reduce(rhoIn, sumOp<scalar>());
//    scalar rhoOut = volPatchIntegrate(rhoName_, outflowPatches_); reduce(rhoOut, sumOp<scalar>());
//    
//    scalar tFlux  = patchFlow(phiName_, outflowPatches_); reduce(tFlux, sumOp<scalar>());
//    scalar tFluxIn  = patchFlow(phiName_, inflowPatches_); reduce(tFluxIn, sumOp<scalar>());
//    
//    scalar cWork  = 0.0;
//    scalar eta    = 0.0;
//    scalar pratio = 0.0;
//    
//    if ( (Pstream::parRun() && Pstream::master()) || !Pstream::parRun() )
//    {
//	
//	Ntotal = Mtotal & omega_;
//	Tin /= inArea;
//	Tout /= outArea;
//	rhoIn /= inArea;
//	rhoOut /= outArea;
//	scalar pDynIn = sqr(tFluxIn / rhoIn / inArea)*rhoIn*0.5;
//	scalar pDynOut = sqr(tFlux / rhoOut / outArea)*rhoOut*0.5;
//	pin /= inArea; pin += pDynIn;
//	pout /= outArea; //pout += pDynOut;
//	
//	cWork = tFlux*heatCapacity_*(Tout - Tin);
//	eta = cWork / Ntotal;
//	pratio = pout / pin;
//	
//	values_[0].append (Ntotal);
//	values_[1].append (Mtotal.x());
//	values_[2].append (Mtotal.y());
//	values_[3].append (Mtotal.z());
//	values_[4].append (eta);
//	values_[5].append (pratio);
//	values_[6].append (tFlux);
//	values_[7].append (Tout);
//    }
}

void Foam::functionObjects::PumpStat::makeFile()
{
    if (Pstream::master())
    {
        if (pumpStatPtr_.valid())
        {
            return;
        }
    }
    
    fileName PumpStatDir;

    if (Pstream::master() && Pstream::parRun())
    {
        PumpStatDir = obr_.time().rootPath() + "/" + obr_.time().caseName().path()  + "/pumpData";
        mkDir(PumpStatDir);
    }
    else if (!Pstream::parRun())
    {
        PumpStatDir = obr_.time().rootPath() + "/" + obr_.time().caseName() + "/pumpData";
        mkDir(PumpStatDir);
    }
    else
    {
    }

    // File update
    if (Pstream::master() || !Pstream::parRun())
    {
        // Create the PumpStat file if not already created
        pumpStatPtr_.reset
        (
            new OFstream
            (
                PumpStatDir + "/" + (name() + "-time.dat")
            )
        );
        
        pumpStatPtr_()<< "Time ";
        
        forAll (vnames_, iName)
        {
            pumpStatPtr_()<< vnames_[iName] << " ";
        }
        
        pumpStatPtr_()<< endl;
    }
}

bool Foam::functionObjects::PumpStat::execute()
{
    return true;
}

bool Foam::functionObjects::PumpStat::write()
{
    // Create the PumpStat file if not already created
    makeFile();

    scalar cTime = obr_.time().value();
    
    if ((cTime < timeStart_) || (cTime > timeEnd_))
    {
        return false;
    }
    
    correct();
        
    if (Pstream::master() || !Pstream::parRun())
    {
        // time history output
        pumpStatPtr_() << (cTime - timeStart_) << " ";
        
        forAll(values_, iValue)
        {
            pumpStatPtr_() << values_[iValue] << " ";
        }
        
        pumpStatPtr_() << endl;
    }
    
    return true;
}




// ************************************************************************* //
