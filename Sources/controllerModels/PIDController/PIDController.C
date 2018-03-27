#include "PIDController.H"
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
    defineTypeNameAndDebug(PIDController, 0);
    addToRunTimeSelectionTable
    (
	controllerModel,
	PIDController,
	dictionary
    );
}
}

namespace Foam
{
namespace fv
{

PIDController::PIDController(const word& name, const word& type, const dictionary& parentDict, const Time& time)
:
    controllerModel (name, type, parentDict, time),
    errorIntegral_(0.0),
    errorDerivative_(0.0),
    oldError_(0.0),
    Kp_(0.0),
    Ki_(0.0),
    Kd_(0.0)
{
    read(coeffs_);
}

PIDController::~PIDController()
{
}

scalar PIDController::newIncrement(const scalar& error, const scalar& cTime)
{
    if (oldTime_ < cTime)
    {
	errorIntegral_ += error * (cTime - oldTime_);
	errorDerivative_ = (error - oldError_) / (cTime - oldTime_);
	oldError_ = error;
	oldTime_  = cTime;
	return Kp_ * error + Ki_ * errorIntegral_ + Kd_ * errorDerivative_;
    }
    else
    {
	FatalErrorIn
	(
	    "scalar PIDController::newIncrement(const scalar& error, const scalar& cTime)"
	)   << "current time is not larger then old time" << endl
	    << "That means that controller is executed again at the same time step " << endl
	    << "It is not allowed in current implementation" << endl
	    << abort(FatalError);
    }
    
    return 0.0;
}

void PIDController::writeData (Ostream& os) const
{
    controllerModel::writeData(os);
}

bool PIDController::read(const dictionary& dict)
{
    if (controllerModel::read(dict))
    {
	Kp_ = readScalar(dict.lookup("Kp"));
	Ki_ = readScalar(dict.lookup("Ki"));
	Kd_ = readScalar(dict.lookup("Kd"));
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

