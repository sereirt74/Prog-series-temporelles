#include "StdAfxRegArchLib.h"
/*!
 \file cVarInMean.cpp
 \brief sources for class cVarInMean methods.

 \author Jean-Baptiste Durand, Ollivier TMAAMASCO
 \date dec-18-2006 - Last change feb-18-2011
*/

namespace RegArchLib {
	cVarInMean::cVarInMean(double theVarInMean):cAbstCondMean(eVarInMean)
	{
		mvVarInMean = theVarInMean ;
		MESS_CREAT("cVarInMean") ;
	}

	/*!
	 * \fn cVarInMean::cVarInMean(cAbstCondMean& theAbstCondMean)
	 * \param cVarInMean cAbstCondMean& theAbstCondMean: the cVarInMean source.
	 */
	cVarInMean::cVarInMean(const cVarInMean& theAbstCondMean) :cAbstCondMean(eVarInMean)
	{
		*this = theAbstCondMean ;
		MESS_CREAT("cVarInMean") ;
	}

	cVarInMean::~cVarInMean()
	{
		MESS_DESTR("cVarInMean") ;
	}

	/*!
	* \fn cAbstCondMean cVarInMean::PtrCopy(void)
	* \param void
	*/
/*	cAbstCondMean* cVarInMean::PtrCopy(void)
	{
		//	cArfima* myArfima = new cArfima(*this);

		//		return myArfima;
		return cAbstCondMeanPtrCopy<cVarInMean>();
	}
*/
	void cVarInMean::Delete(void)
	{
		MESS_DESTR("cVarInMean") ;
	}

#ifndef _RDLL_	
	void cVarInMean::Print(ostream& theOut) const
	{
		theOut << "Variance In Mean model with:" << endl ;
			theOut << "\tDelta = " << mvVarInMean << endl ;
	}
#else		
		void cVarInMean::Print(void)
		{
			Rprintf("Variance In Mean model with:\n");
			Rprintf("\tDelta = %f\n", mvVarInMean);
		}
#endif //_RDLL_
	

	void cVarInMean::SetDefaultInitPoint(double theMean, double theVar)
	{
		mvVarInMean = 0.0 ;
	}

	void cVarInMean::SetDefaultInitPoint(cRegArchValue& theValue)
	{
		mvVarInMean = 0.0;
	}

	void cVarInMean::Set(const cDVector& theVectParam, const uint theNumParam) 
	{
		mvVarInMean = theVectParam[0] ;
	}

	void cVarInMean::Set(const double theValue, const uint theIndex, const uint theNumParam)
	{
		mvVarInMean = theValue ;
	}

	double cVarInMean::Get(const uint theIndex, const uint theNumParam)
	{
		return mvVarInMean ;
	}

	cDVector& cVarInMean::Get(const uint theNumParam)
	{
	cDVector* myAux;
		myAux = new cDVector(1, mvVarInMean);
		return *myAux;
	}

	void cVarInMean::ReAlloc(const uint theSize, const uint theNumParam)
	{
	}

	void cVarInMean::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{
	}

	cVarInMean& cVarInMean::operator =(const cVarInMean& theSrc)
	{
		mvVarInMean = theSrc.mvVarInMean ;
		return *this ;
	} 

	double cVarInMean::ComputeMean(uint theDate, const cRegArchValue& theData) const
	{
		return mvVarInMean * theData.mHt[theDate] ; 
	}

	uint cVarInMean::GetNParam(void) const
	{
		return 1 ;
	}

	uint cVarInMean::GetNLags(void) const
	{
		return 0 ;
	}

	uint cVarInMean::GetNu(void) const
	{
		return 0;
	}

	uint cVarInMean::GetNh(void) const
	{
		return 1;
	}

	void cVarInMean::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{
		theDestVect[theIndex] = mvVarInMean ;
	}

	void  cVarInMean::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
		mvVarInMean = theSrcVect[theIndex] ;
	}

	void cVarInMean::ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, uint theBegIndex, cAbstResiduals* theResiduals)
	{
		theGradData.mCurrentGradMu[theBegIndex] += theData.mHt[theDate] ;
		theGradData.mCurrentGradMu += mvVarInMean * theGradData.mCurrentGradVar ;
	}

	void cVarInMean::GetNParamF(uint theNParam[3]) const
	{
		theNParam[0] = 1;
		theNParam[1] = 0;
		theNParam[2] = 1;
	}

	void cVarInMean::ComputeGradF(uint theDate, const cRegArchValue& theValue, cDVector& theGradF)
	{
		theGradF[0] = theValue.mHt[theDate];
		theGradF[1] = mvVarInMean * 2 *sqrt(theValue.mHt[theDate]);
	}

	void cVarInMean::ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cDVector& theGradTheta, cDVector& theGradU, double& theGradH, uint& theNu, uint& theNh)
	{
		theNh = 1;
		theGradTheta = 0.0;
		theNu = 0 ;
		theGradTheta[theBegIndex] = theValue.mHt[theDate];
		theGradH = mvVarInMean;
	}

	void cVarInMean::ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cFuncMeanAndVar& theDerivM)
	{
		theDerivM.mdFx = 0.0;
		theDerivM.mdFx[theBegIndex] = theValue.mHt[theDate];
		theDerivM.mdFh[0] = mvVarInMean;
	}

	void cVarInMean::ComputeHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData,cRegArchHessien& theHessData, uint theBegIndex, cAbstResiduals* theResiduals)
	{
	uint myNParam = theGradData.GetNParam();
	cDMatrix myMat = cDMatrix(myNParam, myNParam);
		myMat.SetRow(theBegIndex, theGradData.mCurrentGradVar);
		theHessData.mCurrentHessMu += myMat+Transpose(myMat) ;
		theHessData.mCurrentHessMu += mvVarInMean * theHessData.mCurrentHessVar;
	}

	void cVarInMean::ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, const cDVector& theGradTheta, const cDVector& theGradU, double theGradH, uint theNu, uint theNh, cDMatrix& theHessTheta2, cDVector* theHessThetaU, cDVector& theHessThetaH, cDVector* theHessU2, cDVector& theHessUH, double& theHessH2)
	{
//		ComputeGradForM(theDate, theValue, theBegIndex, theGradTheta, theGradU, theGradH, theNu, theNh);
		theHessThetaH[theBegIndex] = 1;
	}

	void cVarInMean::ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, cFuncMeanAndVar& theDerivM)
	{
		theDerivM.md2Fxh[0][theBegIndex] = 1;
	}


	void cVarInMean::ComputeHessF(uint theDate, const cRegArchValue& theData, cDMatrix& theHessF)
	{
		theHessF = 0.0;
		theHessF[0][1] = theHessF[1][0] = 2.0 * sqrt(theData.mHt[theDate]);
	}

	void cVarInMean::GetParamName(uint theIndex, char** theName)
	{
		uint myIndex = theIndex;
		sprintf(theName[myIndex++], "DELTA");
	}

	void cVarInMean::GetParamName(uint theIndex, string theName[])
	{
	uint myIndex = theIndex;
	char myChar[100];
		sprintf(myChar, "DELTA");
		theName[myIndex++] = myChar;
	}
}//namespace
