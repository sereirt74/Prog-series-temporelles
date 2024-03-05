#include "StdAfxRegArchLib.h"
/*!
	\file cAbstResiduals.cpp
	\brief sources for abstract class cAbstResiduals methods.

	\author Jean-Baptiste DURAND, Ollivier TARAMASCO 
	\date dec-18-2006 - Last change feb-18-2011
*/
namespace RegArchLib {
	/*!
	 * \fn cAbstResiduals::cAbstResiduals(eDistrTypeEnum theDistr, cDVector* theModel, bool theSimulFlag)
	 * \param const eDistrTypeEnum theDistr: code for the conditional distribution.
	 * \param const cDVector* theModel: vector of parameters
	 * \param const bool theSimulFlag: true if created for simulation
	 */
	cAbstResiduals::cAbstResiduals(eDistrTypeEnum theDistr, cDVector* theDistrParam, bool theSimulFlag)
	{
		mvDistr = theDistr ;
		if (theSimulFlag)
		{	// random generator simulation
		gsl_rng_env_setup() ;
			mtR = gsl_rng_alloc(gsl_rng_default) ;
		#ifndef _DEBUG
			gsl_rng_set(mtR, (unsigned long int)time(NULL)) ;
		#else
			gsl_rng_set(mtR, 0) ; // Pour avoir toujours la m�me s�rie simul�e quand on teste
		#endif // _DEBUG
		}
		else
			mtR = NULL ;

		if (theDistrParam != NULL)
			mDistrParameter = new cDVector(*theDistrParam) ;
		else
			mDistrParameter = NULL;
		MESS_CREAT("cAbstResiduals") ;
	}



	/*!
	 * \fn cAbstResiduals::~cAbstResiduals()
	 * \details mDistrParameter is deleted here.
	 */
	cAbstResiduals::~cAbstResiduals()
	{
		if (mDistrParameter != NULL)
		{
			delete mDistrParameter;
			mDistrParameter = NULL;
		}
		if (mtR != NULL)
		{	// random generator destructor
			gsl_rng_free(mtR) ;
			mtR = NULL ;
		}
		MESS_DESTR("cAbstResiduals") ;
	}

	/*!
	 * \fn cAbstResiduals::Delete(void)
	 * \param void
	 * \details mDistrParameter is deleted here.
	 */
	void cAbstResiduals::Delete(void)
	{
		if (mDistrParameter != NULL)
		{
			delete mDistrParameter;
			mDistrParameter = NULL;
		}
		if (mtR != NULL)
		{	// random generator destructor
			gsl_rng_free(mtR) ;
			mtR = NULL ;
		}
	}

	/*!
	 * \fn cAbstResiduals::SetSimul(void)
	 * \param void
	 */
	void cAbstResiduals::SetSimul(void)
	{
		if (mtR == NULL)
		{
			// random generator initialisation
			gsl_rng_env_setup() ;
			mtR = gsl_rng_alloc(gsl_rng_default) ;
		#ifndef _DEBUG
			gsl_rng_set(mtR, (unsigned long int)time(NULL)) ;
		#else
			gsl_rng_set(mtR, 0) ; // Pour avoir toujours la m�me s�rie simul�e quand on teste
		#endif // _DEBUG
		}
	}

	/*!
	 * \fn inline eDistrTypeEnum cAbstResiduals::GetDistrType(void) const
	  * \param void
	 */
	eDistrTypeEnum cAbstResiduals::GetDistrType(void) const
	{	return mvDistr ;
	}

	double cAbstResiduals::Get(const uint theIndex)
	{
		if (mDistrParameter->GetSize() > theIndex)
			return (*mDistrParameter)[theIndex] ;
		else
			throw cError("Wrong size in cAbstResiduals::Get") ;
	}

	void cAbstResiduals::Set(double theValue, const uint theIndex)
	{
		if (mDistrParameter->GetSize() <= theIndex)
		{	mDistrParameter->ReAlloc(theIndex+1) ;
		}
		(*mDistrParameter)[theIndex] = theValue ;
	}
	
	void cAbstResiduals::ReAlloc(uint theSize)
	{
		mDistrParameter->ReAlloc(theSize);
	}

#ifndef _RDLL_
	/*!
	 * \fn ostream& cAbstResiduals::operator <<(ostream& theOut, const cAbstResiduals& theAbstResisuals)
	 * \param ostream& theOut: output (file or screen).
	 * \param const cAbstResiduals& theAbstResisuals: the residuals model.
	 * \details Uses cAbstResiduals::Print method.
	 */
	ostream& operator <<(ostream& theOut, const cAbstResiduals& theAbstResisuals)
	{
		theAbstResisuals.Print(theOut) ;
		return theOut ;
	}
#endif // _RDLL_

	void cAbstResiduals::NumericComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, uint theBegIndex, cNumericDerivative& theNumDeriv)
	{

	}


	template<class T>
	T* TemplateCreateRealCondResiduals(cAbstResiduals* theAbstCondResiduals)
	{
	T* mySrc = static_cast<T *>(theAbstCondResiduals);
		if (mySrc)
		{
//			static std::shared_ptr<T> myCondResid(new T(*mySrc));
//			return &*myCondResid;
			return new T(*mySrc);
		}
		else
			throw cError("Wrong contional mean in TemplateCreateOneRealCondMean");
	}
	/*!
	* \fn template<class T> static T* TemplateCreateRealCondResiduals(const cDVector* theDistrParam, const bool theSimulFlag)
	* \param const cDVector* theDistrParam: vector of parameters
	* \param const bool theSimulFlag: true if created for simulation
	*/
	template<class T>
	T* TemplateCreateRealCondResiduals(cDVector* theDistrParam, bool theSimulFlag)
	{
		return new T(theDistrParam, theSimulFlag);
	}
	

	/*!
	* \fn cAbstResiduals* CreateRealCondResiduals(eDistrTypeEnum theType, cDVector* theDistrParam, bool theSimulFlag)
	* \param eDistrTypeEnum theType: type of conditional residuals.
	* \param cDVector* theDistrParam: distribution parameters. Default NULL
	* \param bool theSimulFlag. True if created for simulation. Default true.
	* \details
	* This function has to be changed when adding a new conditional residuals type.
	*/
	cAbstResiduals* CreateRealCondResiduals(eDistrTypeEnum theType, cDVector* theDistrParam, bool theSimulFlag)
	{
		switch (theType)
		{
		case eNormal:
			return TemplateCreateRealCondResiduals<cNormResiduals>(NULL, theSimulFlag);
		break;
		case eStudent:
			if (theDistrParam == NULL)
			{
				cDVector myParam = cDVector(1, 10);
				return TemplateCreateRealCondResiduals<cStudentResiduals>(&myParam, theSimulFlag);
			}
			else
				return TemplateCreateRealCondResiduals<cStudentResiduals>(theDistrParam, theSimulFlag);

		break;

		case eGed:
			if (theDistrParam == NULL)
			{
			cDVector myParam = cDVector(1, 3);
				return TemplateCreateRealCondResiduals<cGedResiduals>(&myParam, theSimulFlag);

			}
			else
				return TemplateCreateRealCondResiduals<cGedResiduals>(theDistrParam, theSimulFlag);
		break;

		case eMixNorm:
			if (theDistrParam == NULL)
			{
			cDVector myParam = cDVector(3, 0.5);
				return TemplateCreateRealCondResiduals<cMixNormResiduals>(&myParam, theSimulFlag);

			}
			else
				return TemplateCreateRealCondResiduals<cMixNormResiduals>(theDistrParam, theSimulFlag);
		break;

		default:
			cError myError("CreateRealCondResiduals: unknown conditional distribution type");
		break;
		}
	}

	/*!
	* \fn cAbstResiduals* CreateRealCondResiduals(eDistrTypeEnum theType, cDVector* theDistrParam, bool theSimulFlag)
	* \param eDistrTypeEnum theType: type of conditional residuals.
	* \param cDVector* theDistrParam: distribution parameters. Default NULL
	* \param bool theSimulFlag. True if created for simulation. Default true.
	* \details
	* This function has to be changed when adding a new conditional residuals type.
	*/
	cAbstResiduals* CreateRealCondResiduals(cAbstResiduals& theAbstCondResiduals)
	{
		switch (theAbstCondResiduals.GetDistrType())
		{
		case eNormal:
			return TemplateCreateRealCondResiduals<cNormResiduals>(&theAbstCondResiduals);
			break;
		case eStudent:
			return TemplateCreateRealCondResiduals<cStudentResiduals>(&theAbstCondResiduals);
			break;

		case eGed:
			return TemplateCreateRealCondResiduals<cGedResiduals>(&theAbstCondResiduals);
			break;

		case eMixNorm:
			return TemplateCreateRealCondResiduals<cMixNormResiduals>(&theAbstCondResiduals);
			break;

		default:
			cError myError("CreateRealCondResiduals: unknown conditional distribution type");
			break;
		}
	}

	/*
	cAbstResiduals* cAbstResiduals::PtrCopy(void) 	
	{
		return CreateRealCondResiduals(*this);
	}
	*/

} // namespace

