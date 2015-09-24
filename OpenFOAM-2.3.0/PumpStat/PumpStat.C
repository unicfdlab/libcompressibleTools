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

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "incompressible/LES/LESModel/LESModel.H"

#include "basicThermo.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "compressible/LES/LESModel/LESModel.H"

#include "FoamFftwDriver.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PumpStat, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::tmp<Foam::scalarField> Foam::PumpStat::normalStress(const word& patchName) const
{

    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    const volScalarField& p = mesh.lookupObject<volScalarField>(pName_);
    
    label patchId = mesh.boundary().findPatchID(patchName);
    
    scalarField pPatch = p.boundaryField()[patchId];
    
    if (p.dimensions() == dimPressure)
    {
	//return tmp<scalarField>
	//(
	//    new scalarField(pPatch)
	//);
    }
    else
    {
	if (rhoRef_ < 0) //density in volScalarField
	{
	    scalarField pRho = mesh.lookupObject<volScalarField>(rhoName_).
				boundaryField()[patchId];
	    pPatch *= pRho;
	}
	else //density is constant
	{
	    pPatch *= rhoRef_;
	}
    }

    return tmp<scalarField>
    (
	new scalarField(pPatch)
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PumpStat::PumpStat
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    probeFreq_(1),
    fftProbeFreq_(1024),
    log_(false),
    momentPatchNames_(0, word::null),
    inflowPatches_(0, word::null),
    outflowPatches_(0, word::null),
    timeStart_(-1.0),
    timeEnd_(-1.0),
    pName_(word::null),
    TName_(word::null),
    phiName_(word::null),
    rhoName_(word::null),
    rhoRef_(1.0),
    heatCapacity_(0.0),
    origin_(vector::zero),
    omega_(vector::zero),
    PumpStatFilePtr_(NULL),
    probeI_(0),
    fftProbeI_(0),
    values_(8),
    vnames_(8, word::null)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "Foam::PumpStat::PumpStat"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating."
            << endl;
    }

    read(dict);
    
    vnames_[0] = "Ntotal";
    vnames_[1] = "Mx";
    vnames_[2] = "My";
    vnames_[3] = "Mz";
    vnames_[4] = "eta";
    vnames_[5] = "pratio";
    vnames_[6] = "flux";
    vnames_[7] = "Tout";
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PumpStat::~PumpStat()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PumpStat::read(const dictionary& dict)
{
    if (!active_)
    {
	return;
    }

    log_ = dict.lookupOrDefault<Switch>("log", false);
    
    if (!log_)
    {
	Info << "Direct logging to stdio disabled" << endl
	    << " to enable, please insert string:" << endl
	    << "log\t\t true;" << endl
	    << "in dictionary" << endl;
    }
    
    dict.lookup("probeFrequency") >> probeFreq_;
    
    dict.lookup("fftProbeFrequency") >> fftProbeFreq_;

    dict.lookup("torquePatchNames") >> momentPatchNames_;
    
    dict.lookup("inflowPatches") >> inflowPatches_;
    
    dict.lookup("outflowPatches") >> outflowPatches_;
    
    dict.lookup("timeStart") >> timeStart_;
    
    dict.lookup("timeEnd") >> timeEnd_;

    dict.lookup("pName") >> pName_;
    
    dict.lookup("TName") >> TName_;
    
    dict.lookup("phiName") >> phiName_;
    
    dict.lookup("rhoName") >> rhoName_;
    
    dict.lookup("rhoRef") >> rhoRef_;
    
    dict.lookup("heatCapacity") >> heatCapacity_;
    
    dict.lookup("origin") >> origin_;
    
    dict.lookup("omega") >> omega_;
}

Foam::scalar Foam::PumpStat::patchesArea(const List<word>& patches)
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    scalar area = 0.0;
    
    forAll(patches, iPatch)
    {
	word patchName = patches[iPatch];
	label patchId = mesh.boundary().findPatchID(patchName);
	if (patchId < 0)
        {
	    FatalError
            << "Unable to find patch " << patchName << " to calculate patch area" << nl
            << exit(FatalError);
        }
        
        area += sum(mesh.magSf().boundaryField()[patchId]);
    }
    
    return area;
}

Foam::scalar Foam::PumpStat::patchFlow(const word& phiName, const List<word>& patches)
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    const surfaceScalarField& phi  = mesh.lookupObject<surfaceScalarField>(phiName);
    
    scalar totFlux = 0.0;
    forAll (patches, iPatch)
    {
	word patchName = patches[iPatch];
	label patchId = mesh.boundary().findPatchID(patchName);
	if (patchId < 0)
        {
	    FatalError
            << "Unable to find patch " << patchName << " to calculate patch integral of the flux: " << phiName  << nl
            << exit(FatalError);
        }
        
        totFlux += sum(phi.boundaryField()[patchId]);
    }
    
    return totFlux;
}

Foam::scalar Foam::PumpStat::volPatchIntegrate(const word& fieldName, const List<word>& patches)
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    const volScalarField& field = mesh.lookupObject<volScalarField>(fieldName);
    
    scalar fieldI = 0.0;
    
    forAll(patches, iPatch)
    {
	word patchName = patches[iPatch];
	label patchId = mesh.boundary().findPatchID(patchName);
	if (patchId < 0)
        {
	    FatalError
            << "Unable to find patch " << patchName << " to calculate patch integral of volField: " << fieldName  << nl
            << exit(FatalError);
        }

	scalarField mSf = mesh.magSf().boundaryField()[patchId];
	scalarField ff  = field.boundaryField()[patchId];
	
	fieldI += sum(mSf*ff);
    }
    
    return fieldI;
}

Foam::vector Foam::PumpStat::gatherMoments()
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    vector Mtotal (0.0, 0.0, 0.0);
    vector rotAxis (omega_ / mag(omega_));
    
    //loop over all patches
    forAll (momentPatchNames_, iPatch)
    {
	word patchName = momentPatchNames_[iPatch];
	label patchId = mesh.boundary().findPatchID(patchName);
	if (patchId < 0)
        {
	    FatalError
            << "Unable to find patch " << patchName << " to calculate moment" << nl
            << exit(FatalError);
        }
        
        //calculate total moment of the patch on each processor
	scalarField pp = normalStress(patchName);
	vectorField Fn = mesh.Sf().boundaryField()[patchId] * pp;
	//vectorField Ft = mesh.Sf().boundaryField()[patchi] & Reff().boundaryField()[patchi];
	//vectorField F  = Fn; /*+  Ft*/;
	//vectorField pCf  = mesh.Cf().boundaryField()[patchId];
	//vectorField Mf = mesh.Cf().boundaryField()[patchId];
	//scalarField t      = ((Mf & rotAxis) - (rotAxis & origin_)) / (rotAxis & rotAxis);
	//vectorField A      = origin_ + t * rotAxis;
	//vectorField arm    = Mf - A;
	vectorField pCf = mesh.Cf().boundaryField()[patchId];
	vectorField B = rotAxis*(rotAxis & pCf);
	vectorField arm = pCf - B;
	vectorField moments= arm ^ Fn;
	
	Mtotal += sum(moments);
    }
    
    return Mtotal;
}

void Foam::PumpStat::correct()
{

    //calculate moments
    vector Mtotal = -gatherMoments();
    reduce (Mtotal, sumOp<vector>());
    
    //calculate power and efficiency
    scalar Ntotal = 0.0;
    scalar inArea = patchesArea(inflowPatches_); reduce (inArea, sumOp<scalar>());
    scalar outArea= patchesArea(outflowPatches_); reduce(outArea, sumOp<scalar>());
    scalar Tin    = volPatchIntegrate(TName_, inflowPatches_); reduce(Tin, sumOp<scalar>());
    scalar Tout   = volPatchIntegrate(TName_, outflowPatches_); reduce(Tout, sumOp<scalar>());
    scalar pin    = volPatchIntegrate(pName_, inflowPatches_); reduce(pin, sumOp<scalar>());
    scalar pout   = volPatchIntegrate(pName_, outflowPatches_); reduce(pout, sumOp<scalar>());
    scalar rhoIn  = volPatchIntegrate(rhoName_, inflowPatches_); reduce(rhoIn, sumOp<scalar>());
    scalar rhoOut = volPatchIntegrate(rhoName_, outflowPatches_); reduce(rhoOut, sumOp<scalar>());
    
    scalar tFlux  = patchFlow(phiName_, outflowPatches_); reduce(tFlux, sumOp<scalar>());
    scalar tFluxIn  = patchFlow(phiName_, inflowPatches_); reduce(tFluxIn, sumOp<scalar>());
    
    scalar cWork  = 0.0;
    scalar eta    = 0.0;
    scalar pratio = 0.0;
    
    if ( (Pstream::parRun() && Pstream::master()) || !Pstream::parRun() )
    {
	scalar 
	
	Ntotal = Mtotal & omega_;
	Tin /= inArea;
	Tout /= outArea;
	rhoIn /= inArea;
	rhoOut /= outArea;
	scalar pDynIn = sqr(tFluxIn / rhoIn / inArea)*rhoIn*0.5;
	scalar pDynOut = sqr(tFlux / rhoOut / outArea)*rhoOut*0.5;
	pin /= inArea; pin += pDynIn;
	pout /= outArea; //pout += pDynOut;
	
	cWork = tFlux*heatCapacity_*(Tout - Tin);
	eta = cWork / Ntotal;
	pratio = pout / pin;
	
	values_[0].append (Ntotal);
	values_[1].append (Mtotal.x());
	values_[2].append (Mtotal.y());
	values_[3].append (Mtotal.z());
	values_[4].append (eta);
	values_[5].append (pratio);
	values_[6].append (tFlux);
	values_[7].append (Tout);
    }
}

void Foam::PumpStat::makeFile()
{

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
	if (PumpStatFilePtr_.empty())
	{
	    // Open new file at start up
	    PumpStatFilePtr_.reset
	    (
		new OFstream
		(
		    PumpStatDir + "/" + (name_ + "-time.dat")
		)
	    );
	    
	    writeFileHeader();
	}
    }
}


void Foam::PumpStat::writeFileHeader()
{
    if (PumpStatFilePtr_.valid())
    {
        PumpStatFilePtr_()
            << "Time ";
        
        forAll (vnames_, iName)
        {
	    PumpStatFilePtr_() << vnames_[iName] << " ";
        }

        PumpStatFilePtr_()<< endl;
    }
}

void Foam::PumpStat::writeFft()
{

    fftProbeI_ = values_[0].size();

    if ( mag(fftProbeI_ % fftProbeFreq_) > VSMALL  )
    {
	return;
    }

    fileName PumpStatDir;

    if (Pstream::master() && Pstream::parRun())
    {
	PumpStatDir = obr_.time().rootPath() + "/" + obr_.time().caseName().path()  + "/pumpData";
    }
    else if (!Pstream::parRun())
    {
	PumpStatDir = obr_.time().rootPath() + "/" + obr_.time().caseName() + "/pumpData";
    }
    
    if (Pstream::master() || !Pstream::parRun())
    {
	const fvMesh& mesh = refCast<const fvMesh>(obr_);
	scalar tau = (mesh.time().value() - timeStart_);
	forAll(vnames_, iName)
	{
	    const DynamicList<scalar>& cvalues = values_[iName];
	    
	    Info << "Executing fft for: " << vnames_[iName] << endl;

	    FoamFftwDriver fftw (cvalues, tau);
	    
	    autoPtr<Pair<List<scalar> > > valFftPtr = fftw.simpleForwardTransform();
	    
	    if (valFftPtr().first().size() > 0)
	    {
		fileName fftFile = PumpStatDir + "/fft-" + vnames_[iName] + ".dat";
		
		OFstream fftStream (fftFile);
		fftStream << "Freq " << vnames_[iName] << endl;
		
		forAll(valFftPtr().first(), k)
		{
		    fftStream << valFftPtr().first()[k] << " " << valFftPtr().second()[k] << endl;
		}
		
		fftStream.flush();
	    }
	}
    }
}

void Foam::PumpStat::execute()
{
    if (!active_)
    {
	return;
    }

    // Create the PumpStat file if not already created
    makeFile();

    scalar cTime = obr_.time().value();
    
    probeI_++;
    
    if ( mag(probeI_ % probeFreq_) > VSMALL  )
    {
	return;
    }
    else
    {
	if (log_)
	{
	    Info << "Evaluating pump statistics" << endl;
	}
	probeI_ = 0;
    }
    
    if ( (cTime < timeStart_) || (cTime > timeEnd_))
    {
	return;
    }
    
    correct();
    
    if (Pstream::master() || !Pstream::parRun())
    {
	label lidx = values_[0].size() - 1;
	// time history output
	PumpStatFilePtr_() << (cTime - timeStart_) << " ";
	
	forAll(values_, iValue)
	{
	    PumpStatFilePtr_() << values_[iValue][lidx] << " ";
	}
	
	PumpStatFilePtr_() << endl;
	
	//fft output
	writeFft();
	
	//output to stdio
	if (log_)
	{
	    label lidx = values_[0].size() - 1;
	    if (lidx >= 0)
	    {
		forAll (vnames_, iName)
		{
		    Info << vnames_[iName] << " = " << values_[iName][lidx] << endl;
		}
	    }
	}
    }
}


void Foam::PumpStat::end()
{
}


void Foam::PumpStat::write()
{
}

void Foam::PumpStat::timeSet()
{
}


// ************************************************************************* //
