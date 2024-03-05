#include "StdAfxRegArchLib.h"

/*!
	\file cAbstCondVar.cpp
	\brief sources for abstract class cAbstCondVar methods.

	\author Jean-Baptiste DURAND, Ollivier TARAMASCO 
	\date dec-18-2006 - Last change feb-18-2011
*/

namespace RegArchLib {

	/*!
	 * \fn cAbstCondVar::cAbstCondVar(eCondVarEnum theType)
	 * \param eCondVarEnum theType: Conditional variance type code. Default eNotKnown.
	 */
	cAbstCondVar::cAbstCondVar(eCondVarEnum theType)
	{
		mvCondVar = theType;
	}

	/*!
	 * \fn cAbstCondVar::~cAbstCondVar()
	 * \details Nothing to do here
	 */
	cAbstCondVar::~cAbstCondVar()
	{
	}

	/*!
	 * \fn inline eCondVarEnum cAbstCondVar::GetCondVarType(void) const
	 * \param void
	 * \brief return the real conditional variance type code.
	*/
	eCondVarEnum cAbstCondVar::GetCondVarType(void) const
	{
		return mvCondVar;
	}

	void cAbstCondVar::SetCondVarType(eCondVarEnum theType)
	{
		mvCondVar = theType;
	}

#ifndef _RDLL_
	/*!
	 * \fn ostream& operator <<(ostream& theOut, const cAbstCondVar& theAbstCondVar)
	 * \param ostream& theOut: output stream (file or screen). Default cout.
	 * \param const cAbstCondVar& theAbstCondVar: the variance model.
	 * \details Uses cAbstCondVar::Print
	 */
	ostream& operator <<(ostream& theOut, const cAbstCondVar& theAbstCondVar)
	{
		theAbstCondVar.Print(theOut);
		return theOut;
	}
#endif // _RDLL_

	/*
		void cAbstCondVar::NumericComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, uint theBegIndex, cAbstResiduals* theResiduals, cNumericDerivative& theNumDeriv)
		{
		uint myNParamTot = theGradData.GetNParam();
		uint myNParam = GetNParam();
		cDVector myParam(myNParamTot);
			RegArchParamToVector(myParam, theBegIndex);
		double myh0 = theNumDeriv.Geth();
		double myF0 = theData.mHt[theDate];
			for (uint i = 0; i < myNParam; i++)
			{
			double myh1 = fabs(myh0*myParam[i + theBegIndex]);
				if (myh1 < 1e-16)
					myh1 = myh0;
			double myF1 = theNumDeriv.mValueForGrad[i + theBegIndex].mHt[theDate];
				theGradData.mCurrentGradVar[theBegIndex + i] = (myF1 - myF0) / myh1;
			}
		}

		void cAbstCondVar::NumericComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex, cAbstResiduals* theResiduals, cNumericDerivative& theNumDeriv)
		{
		uint myNParamTot = theGradData.GetNParam();
		uint myNParam = GetNParam();
			cDVector myParam(myNParamTot);
			RegArchParamToVector(myParam, theBegIndex);
			double myh0 = theNumDeriv.Geth();
			double myF0 = theData.mHt[theDate];
			double* myF1 = new double[myNParam];
			double* myh1 = new double[myNParam];
			for (uint i = 0; i < myNParam; i++)
			{
				myh1[i] = fabs(myh0*myParam[i + theBegIndex]);
				if (myh1[i] < 1e-16)
					myh1[i] = myh0;
				myF1[i] = theNumDeriv.mValueForGrad[i + theBegIndex].mHt[theDate];
				theGradData.mCurrentGradVar[theBegIndex + i] = (myF1[i] - myF0) / myh1[i];
			}
			for (uint i = 0; i < myNParam; i++)
			{
				for (uint j = i; j < myNParam; j++)
				{
					double myh2;
					if (j > i)
						myh2 = fabs(myh0*myParam[j]);
					else
						myh2 = fabs(myh0*(myParam[j] - myh1[j]));
					if (myh2 < 1e-16)
						myh2 = myh0;
					double myF2 = theNumDeriv.mValueForHess[i + theBegIndex][j + theBegIndex].mHt[theDate];
					theHessData.mCurrentHessVar[theBegIndex + i][theBegIndex + j] = theHessData.mCurrentHessVar[theBegIndex + j][theBegIndex + i] = (myF2 - myF1[i] - myF1[j] + myF0) / (myh1[i] * myh2);
				}
			}
			delete[] myh1;
			delete[] myF1;
		}
	*/

	/*!
	* \fn template<class T> static T* TemplateCreateOneRealCondVar(cAbstCondVar& theAbstCondVar)
	* \param cAbstCondVar& theAbstCondVar
	*/
	template<class T>
	T* TemplateCreateOneRealCondVar(cAbstCondVar* theAbstCondVar)
	{
	T*	mySrc = static_cast<T *>(theAbstCondVar);
		if (mySrc)
		{
			return new T(*mySrc);
		}
		else
		{
			throw cError("Wrong Conditional Variance in TemplateCreateOneRealCondVar");
		}
	}

	/*!
	* \fn template<class T> static T* TemplateCreateOneRealCondVar(void)
	* \param void
	*/
	template<class T>
	T* TemplateCreateOneRealCondVar(void)
	{
		return new T();
	}

	/*!
	* \fn CreateOneRealCondVar* CreateOneRealCondVar(eCondVarEnum theType)
	* \param theType: type of conditional variance.
	* \par Details This function has to be changed when adding a new conditional variance type.
	*/
	cAbstCondVar* CreateOneRealCondVar(eCondVarEnum theType)
	{
		switch (theType)
		{
		case eCste:
			return TemplateCreateOneRealCondVar<cConstCondVar>();
			break;
		case eArch:
			return TemplateCreateOneRealCondVar<cArch>();
			break;
		case eGarch:
			return TemplateCreateOneRealCondVar<cGarch>();
			break;
		case eNgarch:
			return TemplateCreateOneRealCondVar<cNgarch>();
			break;
		case eEgarch:
			return TemplateCreateOneRealCondVar<cEgarch>();
			break;
		case eAparch:
			return TemplateCreateOneRealCondVar<cAparch>();
			break;
		case eTarch:
			return TemplateCreateOneRealCondVar<cTarch>();
			break;
		case eFigarch:
			return TemplateCreateOneRealCondVar<cFigarch>();
			break;
		case eUgarch:
			return TemplateCreateOneRealCondVar<cUgarch>();
			break;
		case eGtarch:
			return TemplateCreateOneRealCondVar<cGtarch>();
			break;
		case eTsgarch:
			return TemplateCreateOneRealCondVar<cTsgarch>();
			break;
		case eNagarch:
			return TemplateCreateOneRealCondVar<cNagarch>();
			break;
		case eSqrgarch:
			return TemplateCreateOneRealCondVar<cSqrgarch>();
			break;
		case eStgarch:
			return TemplateCreateOneRealCondVar<cStgarch>();
			break;
		case eLoggarch:
			return TemplateCreateOneRealCondVar<cLoggarch>();
			break;
		case eNotKnown:
		default:
			throw cError("unknown conditional variance type");
			break;
		}
	}

	cAbstCondVar* CreateOneRealCondVar(cAbstCondVar& theAbstCondVar) 
	{
		switch (theAbstCondVar.GetCondVarType())
		{
		case eCste:
			return TemplateCreateOneRealCondVar<cConstCondVar>(&theAbstCondVar);
			break;
		case eArch:
			return TemplateCreateOneRealCondVar<cArch>(&theAbstCondVar);
			break;
		case eGarch:
			return TemplateCreateOneRealCondVar<cGarch>(&theAbstCondVar);
			break;
		case eNgarch:
			return TemplateCreateOneRealCondVar<cNgarch>(&theAbstCondVar);
			break;
		case eEgarch:
			return TemplateCreateOneRealCondVar<cEgarch>(&theAbstCondVar);
			break;
		case eAparch:
			TemplateCreateOneRealCondVar<cAparch>(&theAbstCondVar);
			break;
		case eTarch:
			return TemplateCreateOneRealCondVar<cTarch>(&theAbstCondVar);
			break;
		case eFigarch:
			return TemplateCreateOneRealCondVar<cFigarch>(&theAbstCondVar);
			break;
		case eUgarch:
			return TemplateCreateOneRealCondVar<cUgarch>(&theAbstCondVar);
			break;		
		case eGtarch:
			return TemplateCreateOneRealCondVar<cGtarch>(&theAbstCondVar);
			break;
		case eTsgarch:
			return TemplateCreateOneRealCondVar<cTsgarch>(&theAbstCondVar);
			break;		
		case eNagarch:
			return TemplateCreateOneRealCondVar<cNagarch>(&theAbstCondVar);
			break;
		case eSqrgarch:
			return TemplateCreateOneRealCondVar<cSqrgarch>(&theAbstCondVar);
			break;
		case eStgarch:
			return TemplateCreateOneRealCondVar<cStgarch>(&theAbstCondVar);
			break;
		case eLoggarch:
			return TemplateCreateOneRealCondVar<cLoggarch>(&theAbstCondVar);
			break;
		case eNotKnown:
		default:
			throw cError("unknown conditional variance type");
			break;
		}
	}

/*
	cAbstCondVar* cAbstCondVar::PtrCopy(void) 
	{
		return CreateOneRealCondVar(*this);
	}
*/
}
