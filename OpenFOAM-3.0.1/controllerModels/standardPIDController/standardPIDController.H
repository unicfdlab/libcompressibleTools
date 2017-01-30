#ifndef standardPIDController_H
#define standardPIDController_H

#include "controllerModel.H"

namespace Foam
{

namespace fv
{


class standardPIDController 
: public controllerModel
{

private:

    //- forbid default constructor
    standardPIDController();
    
    //- forbid copy constructor
    standardPIDController(const standardPIDController& );

protected:

    //-
    scalar cError_;

    //-
    scalar cErrorInt_;
    
    //-
    scalar cErrorDer_;

    //-
    scalar errorIntegral_;
    
    //-
    scalar errorDerivative_;

    //-
    scalar oldError_;
    
    //-
    scalar Kp_;
    
    //-
    scalar Ti_;
    
    //-
    scalar Td_;

public:

    //-
    TypeName("standardPIDController");

    //- Construct from components
    standardPIDController
    (
	const word& name,
	const word& type,
	const dictionary& parentDict,
	const Time& time
    );
    
    //-
    virtual ~standardPIDController();

    //-
    virtual autoPtr<controllerModel> clone() const
    {
	autoPtr<standardPIDController> nPIDController
	(
	    dynamic_cast<standardPIDController*>
	    (
		controllerModel::clone().ptr()
	    )
	);
	
	nPIDController->Kp_ = this->Kp_;
	nPIDController->Ti_ = this->Ti_;
	nPIDController->Td_ = this->Td_;
	
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

