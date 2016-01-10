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

Class
    Foam::PumpStat

Description
    
    Calculatates integral statistics for pump:
    - instant moment
    - instant power
    - coefficient of efficiency
    - pressure ratio
    - mass flow rate
    
    Writes data to specified file as the time history and fft

    Input data:
    - rotating patches
    - axis of rotation
    - rotation speed (in rpm)
    - flow and pressure patches
    - heat capacity of media

SourceFiles
    PumpStat.C
    IOPumpStat.H

\*---------------------------------------------------------------------------*/

#ifndef PumpStat_H
#define PumpStat_H

#include "fvMesh.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "Switch.H"
#include "pointFieldFwd.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PumpStat Declaration
\*---------------------------------------------------------------------------*/

class PumpStat
{
public:

protected:

    // Private data

        //- Name of this set of PumpStat,
        //  Also used as the name of the probes directory.
        word name_;
	
	//-
        const objectRegistry& obr_;

        //- on/off switch
        bool active_;
        
        //- probe frequency
        label probeFreq_;
        
        //-
        label fftProbeFreq_;

        //- Switch to send output to Info as well as to file
        Switch log_;

        // Read from dictionary

	//- Patch to integrate PumpStat over
        List<word> momentPatchNames_;
        
        //-
        List<word> inflowPatches_;
        
        //-
        List<word> outflowPatches_;
        
        //-
        scalar timeStart_;
        
        //-
        scalar timeEnd_;

        //- Name of the pressure field
        word pName_;
        
        //- Name of the temperature field
        word TName_;
        
        //- Name of the mass flux field
        word phiName_;

        //- Name of density field 
        word rhoName_;

	//- Reference density needed for incompressible calculations
	scalar rhoRef_;

        //-
        scalar heatCapacity_;
        
        //- origin of rotation
        vector origin_;
        
        //- speed of rotation (in radians)
        vector omega_;

        //- PumpStat/moment file ptr
        autoPtr<OFstream> PumpStatFilePtr_;
        
        //-
        label probeI_;
        
        //-
        label fftProbeI_;
        
        //-
        List<DynamicList<scalar> > values_;
        
        //-
        List<word> vnames_;


    // Private Member Functions

        //- If the PumpStat file has not been created create it
        void makeFile();

        //- Output file header information
        virtual void writeFileHeader();
        
	//-
	scalar patchFlow (const word&, const List<word>& );

        //-
        scalar patchesArea(const List<word>&);
        
        //-
        scalar volPatchIntegrate(const word&, const List<word>& );
        
        //-
        vector gatherMoments ();
        
        //-
        void correct();
        
        //-
        void writeFft();

	//- Returns normal stresses using pressure (and optionally density) field
	tmp<scalarField> normalStress(const word& patchName) const;

        //- Disallow default bitwise copy construct
        PumpStat(const PumpStat&);

        //- Disallow default bitwise assignment
        void operator=(const PumpStat&);


public:

    //- Runtime type information
    TypeName("PumpStat");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        PumpStat
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );

	//- Destructor
	virtual ~PumpStat();


    // Member Functions

        //- Return name of the set of PumpStat
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the PumpStat data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Write the PumpStat
        virtual void write();
	
        //- Update for changes of mesh
	virtual void updateMesh(const mapPolyMesh&)
        {};
	
	//- Called when time was set at the end of the Time::operator++
	virtual void timeSet();
	
	//- Update for changes of mesh
	virtual void movePoints(const polyMesh&)
	{};
	
	//- Update for changes of mesh due to readUpdate
	virtual void readUpdate(const polyMesh::readUpdateState state)
	{};

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
