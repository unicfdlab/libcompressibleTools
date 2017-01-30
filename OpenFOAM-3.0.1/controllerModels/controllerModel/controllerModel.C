#include "controllerModel.H"
#include "volFields.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(controllerModel, 0);
    defineRunTimeSelectionTable(controllerModel, dictionary);
}
}

namespace Foam
{
namespace fv
{

autoPtr<controllerModel> controllerModel::New
(
    const word& name,
    const dictionary& parentDict,
    const Time& time
)
{

    word modelType(parentDict.lookup("type"));

    Info<< "Selecting controllerModel type " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
	dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
	FatalErrorIn
	(
	    "controllerModel::New(const word&, const dictionary&, const fvMesh&)"
	)   << "Unknown Model type " << modelType << nl << nl
	<< "Valid model types are:" << nl
	<< dictionaryConstructorTablePtr_->sortedToc()
	<< exit(FatalError);
    }

    return autoPtr<controllerModel>(cstrIter()(name, modelType, parentDict, time));
}

controllerModel::controllerModel(const word& name, const word& type, const dictionary& parentDict, const Time& time)
:
    name_(name),
    cmType_(type),
    runTime_(time),
    parent_(parentDict),
    coeffs_(parentDict.subDict(cmType_ + "Coeffs")),
    oldTime_(-1.0)
{
}

controllerModel::~controllerModel()
{
}

void controllerModel::writeData (Ostream& os) const
{
    os << nl;
    os << token::BEGIN_BLOCK << nl;
    os.incrIndent();
    os << "type" << token::TAB << cmType_ << token::END_STATEMENT << nl << nl;
    os << (cmType_ + "Coeffs") << nl;
    coeffs_.write(os);
    os.decrIndent();
    os << token::END_BLOCK << nl;
}

bool controllerModel::read(const dictionary& dict)
{
    return true;
}

scalar controllerModel::error(const scalar& reference, const scalar& actualValue)
{
    return reference - actualValue;
}


const word& controllerModel::name() const
{
    return name_;
}

}; //namespace fv

}; //namespace Foam


//END-OF-FILE

