#include "CompressibleTurbulenceModel.H"
#include "compressibleTransportModel.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "TurbulenceModels/compressible/turbulentFluidThermoModels/makeTurbulenceModel.H"

#include "TurbulenceModels/compressible/ThermalDiffusivity/ThermalDiffusivity.H"
#include "TurbulenceModels/compressible/EddyDiffusivity/EddyDiffusivity.H"

#include "laminar.H"
#include "RASModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeBaseTurbulenceModel
(
    geometricOneField,
    volScalarField,
    compressibleTurbulenceModel,
    CompressibleTurbulenceModel,
    ThermalDiffusivity,
    fluidThermo
);

#define makeRASModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (fluidThermoCompressibleTurbulenceModel, RAS, Type)

#include "kEpsilonSarkarZeman.H"
makeRASModel(kEpsilonSarkarZeman);

#include "jetPrandtlSafronov.H"
makeRASModel(jetPrandtlSafronov)

//
//END-OF-FILE
//


