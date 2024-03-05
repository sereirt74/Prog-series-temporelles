#include "StdAfxRegArchLib.h"
/*!
 \file cStdDevInMean.cpp
 \brief sources for class cStdDevInMean methods.

 \author Jean-Baptiste Durand, Ollivier TMAAMASCO
 \date dec-18-2006 - Last change feb-18-2011
*/
namespace RegArchLib {
	cStdDevInMean::cStdDevInMean(double theStdDevInMean):cAbstCondMean(eStdDevInMean)
	{
		mvStdDevInMean = theStdDevInMean ;
		MESS_CREAT("cStdDevInMean") ;
	}

	/*!
	 * \fn cStdDevInMean::cStdDevInMean(cAbstCondMean& theAbstCondMean)
	 * \param const cAbstCondMean& theAbstCondMean: the cStdDevInMean source.
	 */
	cStdDevInMean::cStdDevInMean(const cStdDevInMean& theStdDevInMean) :cAbstCondMean(eStdDevInMean)
	{
		*this = theStdDevInMean;
		MESS_CREAT("cStdDevInMean") ;
	}

	cStdDevInMean::~cStdDevInMean()
	{
		MESS_DESTR("cStdDevInMean") ;
	}

	/*!
	 * \fn cAbstCondMean cStdDevInMean::PtrCopy(void)
	 * \param void
	 */
/*	cAbstCondMean* cStdDevInMean::PtrCopy(void)
	{
		//	cArfima* myArfima = new cArfima(*this);

		//		return myArfima;
		return cAbstCondMeanPtrCopy<cStdDevInMean>();
	}
*/
	void cStdDevInMean::Delete(void)
	{
		MESS_DESTR("cStdDevInMean") ;
	}

#ifndef _RDLL_
	void cStdDevInMean::Print(ostream& theOut) const
	{
		theOut << "Standard Deviation In Mean model with:" << endl ;
		theOut << "\tDelta = " << mvStdDevInMean << endl ;
	}
#else
	void cStdDevInMean::Print(void)
	{
		Rprintf("Standard Deviation In Mean model with:%d\n") ;
		Rprintf("\tDelta = %f\n", mvStdDevInMean) ;
	}
#endif // _RDLL_


	void cStdDevInMean::SetDefaultInitPoint(double theMean, double theVar)
	{
		mvStdDevInMean = 0.0 ;
	}

	void cStdDevInMean::SetDefaultInitPoint(cRegArchValue& theValue)
	{
		mvStdDevInMean = 0.0;
	}

	void cStdDevInMean::Set(const cDVector& theVectParam, const uint theNumParam)
	{
		mvStdDevInMean = theVectParam[0] ;
	}

	void cStdDevInMean::Set(const double theValue, const uint theIndex, const uint theNumParam) 
	{
		mvStdDevInMean = theValue ;
	}

	double cStdDevInMean::Get(const uint theIndex, const uint theNumParam)
	{
		return mvStdDevInMean ;
	}

	cDVector& cStdDevInMean::Get(const uint theNumParam)
	{
	cDVector* myAux;
		myAux = new cDVector(1, mvStdDevInMean);
		return *myAux;
	}

	void cStdDevInMean::ReAlloc(const uint theSize, const uint theNumParam)
	{
	}

	void cStdDevInMean::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{
	}

	cStdDevInMean& cStdDevInMean::operator =(const cStdDevInMean& theSrc)
	{
		mvStdDevInMean = theSrc.mvStdDevInMean ;
		return *this ;
	} 

	double cStdDevInMean::ComputeMean(uint theDate, const cRegArchValue& theData) const
	{
		return mvStdDevInMean * sqrt(theData.mHt[theDate]) ; 
	}

	uint cStdDevInMean::GetNParam(void) const
	{
		return 1 ;
	}
	
	uint cStdDevInMean::GetNLags(void) const
	{
		return 0 ;
	}

	uint cStdDevInMean::GetNu(void) const
	{
		return 0;
	}

	uint cStdDevInMean::GetNh(void) const
	{
		return 1;
	}

	void cStdDevInMean::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{
		theDestVect[theIndex] = mvStdDevInMean ;
	}

	void  cStdDevInMean::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
		mvStdDevInMean = theSrcVect[theIndex] ;
	}

	void cStdDevInMean::ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, uint theBegIndex, cAbstResiduals* theResiduals)
	{
		theGradData.mCurrentGradMu[theBegIndex] += sqrt(theData.mHt[theDate]) ;
		theGradData.mCurrentGradMu += mvStdDevInMean * theGradData.mCurrentGradSigma ;
	}

	void cStdDevInMean::GetNParamF(uint theNParam[3]) const
	{
		theNParam[0] = 1;
		theNParam[1] = 0;
		theNParam[2] = 1;
	}

	void cStdDevInMean::ComputeGradF(uint theDate, const cRegArchValue& theValue, cDVector& theGradF)
	{
		theGradF[0] = sqrt(theValue.mHt[theDate]);
		theGradF[1] = mvStdDevInMean;
	}

	void cStdDevInMean::ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cDVector& theGradTheta, cDVector& theGradU, double& theGradH, uint& theNu, uint& theNh)
	{
		theNh = 1;
		theNu = 0;
		double myAux = theGradTheta[theBegIndex] = sqrt(theValue.mHt[theDate]);
		theGradH = mvStdDevInMean/(2*myAux);
	}

	void cStdDevInMean::ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cFuncMeanAndVar& theDerivM)
	{
		theDerivM.mdFx = 0.0;
		double myAux = theDerivM.mdFx[theBegIndex] = sqrt(theValue.mHt[theDate]);
			theDerivM.mdFh[0] = mvStdDevInMean / (2 * myAux);
	}

	void cStdDevInMean::ComputeHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData,cRegArchHessien& theHessData, uint theBegIndex, cAbstResiduals* theResiduals)
	{
	uint myNParam = theGradData.GetNParam();
	cDMatrix myMat = cDMatrix(myNParam, myNParam);
		myMat.SetRow(theBegIndex, theGradData.mCurrentGradSigma);
		
		theHessData.mCurrentHessMu += myMat + Transpose(myMat) ;
		theHessData.mCurrentHessMu += mvStdDevInMean * theHessData.mCurrentHessSigma;
	}

	void cStdDevInMean::ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, const cDVector& theGradTheta, const cDVector& theGradU, double theGradH, uint theNu, uint theNh, cDMatrix& theHessTheta2,cDVector* theHessThetaU, cDVector& theHessThetaH, cDVector* theHessU2, cDVector& theHessUH, double& theHessH2)
	{
		theHessThetaH[theBegIndex] = 1/(2*sqrt(theValue.mHt[theDate]));
	}

	void cStdDevInMean::ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, cFuncMeanAndVar& theDerivM)
	{
		theDerivM.md2Fxh[0][theBegIndex] = 1 / (2 * sqrt(theValue.mHt[theDate]));
		theDerivM.md2Fhh[0][0] = -mvStdDevInMean / (4 * pow(theValue.mHt[theDate], 1.5));
	}

	void cStdDevInMean::ComputeHessF(uint theDate, const cRegArchValue& theData, cDMatrix& theHessF)
	{
		theHessF = 0.0;
		theHessF[0][1] = theHessF[1][0] = 1.0;
	}
	
	void cStdDevInMean::GetParamName(uint theIndex, char** theName)
	{
		uint myIndex = theIndex;
		sprintf(theName[myIndex++], "DELTA");
	}

	void cStdDevInMean::GetParamName(uint theIndex, string theName[])
	{
	uint myIndex = theIndex;
	char myChar[100];

		sprintf(myChar, "DELTA");
		theName[myIndex++] = myChar;
	}

}//namespace
