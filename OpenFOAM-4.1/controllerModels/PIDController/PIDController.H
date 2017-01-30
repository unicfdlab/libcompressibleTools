#ifndef PIDController_H
#define PIDController_H

#include "controllerModel.H"

namespace Foam
{

namespace fv
{


class PIDController 
: public controllerModel
{

private:

    //- forbid default constructor
    PIDController();
    
    //- forbid copy constructor
    PIDController(const PIDController& );

protected:

    //-
    scalar errorIntegral_;
    
    //-
    scalar errorDerivative_;

    //-
    scalar oldError_;
    
    //-
    scalar Kp_;
    
    //-
    scalar Ki_;
    
    //-
    scalar Kd_;
    

public:

    //-
    TypeName("PIDController");

    //- Construct from components
    PIDController
    (
	const word& name,
	const word& type,
	const dictionary& parentDict,
	const Time& time
    );
    
    //-
    virtual ~PIDController();
    
    //-
    virtual autoPtr<controllerModel> clone() const
    {
	autoPtr<PIDController> nPIDController
	(
	    dynamic_cast<PIDController*>
	    (
		controllerModel::clone().ptr()
	    )
	);
	
	nPIDController->Kp_ = this->Kp_;
	nPIDController->Ki_ = this->Ki_;
	nPIDController->Kd_ = this->Kd_;
	
	nPIDController->errorIntegral_ = this->errorIntegral_;
	nPIDController->errorDerivative_ = this->errorDerivative_;
	nPIDController->oldError_ = this->oldError_;
	
	return autoPtr<controllerModel>(nPIDController.ptr());
    }

    //-
    virtual scalar newIncrement(const scalar& error, const scalar& cTime);
    
    // I/O
    
    //-
    virtual void writeData (Ostream& ) const;
    
    //-
    virtual bool read (const dictionary& dict);
};

};

};

#endif

//END-OF-FILE

