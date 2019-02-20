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

#include "Time.H"
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
#include "forces.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class PumpStat Declaration
\*---------------------------------------------------------------------------*/

class PumpStat
    : public forces
{
public:

protected:

    // Private data

        //-
        word phiName_;
        
        //-
        labelHashSet inflowPatches_;
        
        //-
        labelHashSet outflowPatches_;
        
        //-
        scalar timeStart_;
        
        //-
        scalar timeEnd_;
        
        //- speed of rotation (in radians)
        vector omega_;
        
        //-
        List<scalar> values_;
        
        //-
        List<word> vnames_;
        
        //-
        autoPtr<OFstream> pumpStatPtr_;


    // Private Member Functions
    
        //-
        void initialize();

        //- If the PumpStat file has not been created create it
        void makeFile();
        
        //-
        scalar patchFlow (const word&, const labelHashSet& );
        
        //-
        scalar patchH (const labelHashSet& );
        
        //-
        scalar patchP (const labelHashSet& );

        //-
        void correct();

        //-
        tmp<volSymmTensorField> totalStress();
        
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
            const Time&,
            const dictionary&
        );

        //-
        PumpStat
        (
            const word& name,
            const objectRegistry&,
            const dictionary&
        );

        //- Destructor
        virtual ~PumpStat();

    // Member Functions

        //- Read the PumpStat data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the PumpStat
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
