#ifndef controllerModel_H
#define controllerModel_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "HashTable.H"



namespace Foam
{

class fvMesh;
class Time;

namespace fv
{


class controllerModel
{

private:

    //- forbid default constructor
    controllerModel();
    
    //- forbid copy constructor
    controllerModel(const controllerModel& );

protected:
    
    //-
    const word name_;
    
    //-
    const word cmType_;
    
    //-
    const Time& runTime_;
    
    //-
    const dictionary parent_;

    //-
    const dictionary coeffs_;
    
    //-
    scalar oldTime_;

public:

    //-
    TypeName("controllerModel");
    
    declareRunTimeSelectionTable
    (
	autoPtr,
	controllerModel,
	dictionary,
	(
	    const word& name,
	    const word& type,
	    const dictionary& parentDict,
	    const Time& runTime
	),
	(name, type, parentDict, runTime)
    );

    //- Construct from components
    controllerModel
    (
	const word& name,
	const word& type,
	const dictionary& parentDict,
	const Time& runTime
    );

    //- Return clone
    virtual autoPtr<controllerModel> clone() const
    {
	autoPtr<controllerModel> nControllerModel = 
	    controllerModel::New
	    (
		this->name_,
		this->parent_,
		this->runTime_
	    );
	    
	nControllerModel->oldTime_ = this->oldTime_;
	
	return nControllerModel;
    }

    //- Return pointer to new controllerModel object created
    //  on the freestore from an Istream
    class iNew
    {
	//- Reference to the name
	const word& name_;
	
	//-
	const word& type_;
	
	//-
	const Time& time_;

    public:
	iNew
	(
	    const word& name,
	    const word& cmType,
	    const Time& time
	)
	:
	name_(name),
	type_(cmType),
	time_(time)
	{}

	autoPtr<controllerModel> operator()(Istream& is) const
	{
	    //const word name(is);
	    const dictionary parentDict(is);

	    return autoPtr<controllerModel>
	    (
		controllerModel::New(name_, parentDict, time_)
	    );
	}
    };

    //- Return a reference to the selected fvOption model
    static autoPtr<controllerModel> New
    (
	const word& name,
	const dictionary& parentDict,
	const Time& runTime
    );
    
    //-
    virtual ~controllerModel();
    
    //-
    virtual scalar error(const scalar& reference, const scalar& actualValue);
    
    //-
    virtual scalar newIncrement(const scalar& error, const scalar& cTime) = 0;
    
    //-
    const word& name() const;
    
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

