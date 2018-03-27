#include "standardPIDController.H"
#include "volFields.H"
#include "fvMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(standardPIDController, 0);
    addToRunTimeSelectionTable
    (
	controllerModel,
	standardPIDController,
	dictionary
    );
}
}

namespace Foam
{
namespace fv
{

standardPIDController::standardPIDController(const word& name, const word& type, const dictionary& parentDict, const Time& time)
:
    controllerModel (name, type, parentDict, time),
    cError_(0.0),
    cErrorInt_(0.0),
    cErrorDer_(0.0),
    errorIntegral_(0.0),
    errorDerivative_(0.0),
    oldError_(0.0),
    Kp_(0.0),
    Ti_(1.0),
    Td_(0.0)
{
    read(coeffs_);
}

standardPIDController::~standardPIDController()
{
}

scalar standardPIDController::newIncrement(const scalar& error, const scalar& cTime)
{
    if (oldTime_ < cTime)
    {
        oldError_  = cError_;
        cError_ = 0.0;
        oldTime_  = cTime;
    }

    scalar
    deltaT = runTime_.deltaTValue();
    
    cErrorInt_ += (-cError_*deltaT + error * deltaT);
    cErrorDer_ =  (error - oldError_) / deltaT;
    
    errorIntegral_ += cErrorInt_;
    errorDerivative_ = cErrorDer_;
    
    cError_ = error;
    
    return Kp_ * (error + (1.0 / Ti_) * errorIntegral_ + Td_ * errorDerivative_);
}

void standardPIDController::writeData (Ostream& os) const
{
    controllerModel::writeData(os);
}

bool standardPIDController::read(const dictionary& dict)
{
    if (controllerModel::read(dict))
    {
	Kp_ = readScalar(dict.lookup("Kp"));
	Ti_ = readScalar(dict.lookup("Ti"));
	Td_ = readScalar(dict.lookup("Td"));
	return true;
    }
    else
    {
	return false;
    }
    
    return true;
}

}; //namespace fv

}; //namespace Foam


//END-OF-FILE

