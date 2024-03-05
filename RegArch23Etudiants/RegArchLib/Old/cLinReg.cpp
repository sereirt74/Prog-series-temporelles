#include "StdAfxRegArchLib.h"
/*!
	\file cLinReg.cpp
	\brief sources for class cLinReg methods.

	\author Jean-Baptiste DURAND, Ollivier TARAMASCO 
	\date dec-18-2006 - Last change feb-18-2011
*/

namespace RegArchLib {
	cLinReg::cLinReg(int theNLinReg):cAbstCondMean(eLinReg)
	{
		mvBeta.ReAlloc(theNLinReg) ;
		MESS_CREAT("cLinReg")
	}

	cLinReg::cLinReg(const cDVector& theBeta):cAbstCondMean(eLinReg)
	{
		mvBeta = theBeta ;
		MESS_CREAT("cLinReg")
	}

	/*!
	 * \fn cLinReg::cLinReg(cAbstCondMean& theAbstCondMean)
	 * \param const cAbstCondMean& theAbstCondMean: the cLinReg source.
	 */
	cLinReg::cLinReg(const cLinReg& theLinReg) :cAbstCondMean(eLinReg)
	{
		*this = theLinReg;
		MESS_CREAT("cLinReg")
	}

	cLinReg::~cLinReg()
	{
		mvBeta.Delete() ;
		MESS_DESTR("cLinReg") ;
	}

	/*!
	 * \fn cAbstCondMean cLinReg::PtrCopy(void)
	 * \param void
	 */
/*
cAbstCondMean* cLinReg::PtrCopy(void)
	{
		//	cArfima* myArfima = new cArfima(*this);

		//		return myArfima;
		return cAbstCondMeanPtrCopy<cLinReg>();
	}
*/
	void cLinReg::Delete(void)
	{
		mvBeta.Delete() ;
	}

#ifndef _RDLL_
	void cLinReg::Print(ostream& theOut) const
	{	cout << "Linear Regression with:" << endl ;
		for (uint i = 0 ; i < mvBeta.GetSize() ; i++)
				theOut << "\tLinReg[" << i+1 << "] = " << mvBeta[i] << endl ;
	}
#else
	void cLinReg::Print(void)
	{
		Rprintf("Linear Regression with:\n");
		for (uint i = 0; i < mvBeta.GetSize(); i++)
			Rprintf("\tLinReg[%d]=%f", i + 1, mvBeta[i]);

	}
#endif // _RDLL_

	void cLinReg::SetDefaultInitPoint(double theMean, double theVar)
	{
		for (uint i = 0 ; i < mvBeta.GetSize() ; i++)
			mvBeta[i] = 0.0 ;
	}

	void cLinReg::SetDefaultInitPoint(cRegArchValue& theValue)
	{
		mvBeta = Inv(Transpose(theValue.mXt)*theValue.mXt)*Transpose(theValue.mXt)*theValue.mYt;
	}

	void cLinReg::Set(const cDVector& theVectParam, const uint theNumParam)
	{
		mvBeta = theVectParam ;
	}

	void cLinReg::Set(const double theValue, const uint theIndex, const uint theNumParam)
	{
		mvBeta[theNumParam] = theValue ;
	}

	double cLinReg::Get(const uint theIndex, const uint theNumParam)
	{
		return mvBeta[theIndex] ;
	}

	cDVector& cLinReg::Get(const uint theNumParam)
	{
		return mvBeta;
	}

	void cLinReg::ReAlloc(const cDVector& theModel, const uint theNumParam)
	{
		mvBeta.ReAlloc((int)theModel[0]) ;
	}

	void cLinReg::ReAlloc(const uint theModel, const uint theNumParam)
	{
		mvBeta.ReAlloc(theModel) ;
	}

	cLinReg& cLinReg::operator =(const cLinReg &theSrc)
	{
		mvBeta =theSrc.mvBeta ;
		return *this ;
	} 

	double cLinReg::ComputeMean(uint theDate, const cRegArchValue& theData) const
	{
		int myp = (int)mvBeta.GetSize() ;
		if (myp > 0)
		{	
			double myRes = 0.0 ;
				for (int i = 0 ; i < myp ; i++)
						myRes += mvBeta[i] * theData.mXt[theDate][i] ; 
				return myRes ;
		}
		else 
			return 0.0 ;
	}

	uint cLinReg::GetNParam(void) const
	{
		return (int)mvBeta.GetSize() ;
	}

	uint cLinReg::GetNu(void) const
	{
		return 0 ;
	}

	uint cLinReg::GetNLags(void) const
	{
		return 0;
	}

	uint cLinReg::GetNh(void) const
	{
		return 0;
	}

	void cLinReg::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{
	uint mySize = mvBeta.GetSize() ;
		if (theDestVect.GetSize() < mySize + theIndex)
			throw cError("Wrong size") ;
		mvBeta.SetSubVectorWithThis(theDestVect, theIndex) ;
	}

	void cLinReg::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	uint mySize = theSrcVect.GetSize() ;
		if (mvBeta.GetSize() + theIndex > mySize)
			throw cError("Wrong size") ;
		mvBeta.SetThisWithSubVector(theSrcVect, theIndex) ;
	}

	void cLinReg::ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, uint theBegIndex, cAbstResiduals* theResiduals)
	{	for (int i = 0 ; i < (int)mvBeta.GetSize()  ; i++)
			theGradData.mCurrentGradMu[theBegIndex+i] += theData.mXt[theDate][i];
	}

	void cLinReg::GetNParamF(uint theNParam[3]) const
	{
		theNParam[0] = mvBeta.GetSize();
		theNParam[1] = theNParam[2] = 0;
	}
	
	void cLinReg::ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cDVector& theGradTheta, cDVector& theGradU, double& theGradH, uint& theNu, uint& theNh)
	{
		theNh = theNu = 0;
		theGradTheta = 0.0;
		theGradH = 0.0;
		uint myNBeta =  mvBeta.GetSize();
		for (uint i = 0; i < myNBeta; i++)
			theGradTheta[theBegIndex + i] = theValue.mXt[theDate][i];
	}

	void cLinReg::ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cFuncMeanAndVar& theDerivM)
	{
		theDerivM.mdFx = 0.0;
		uint myNBeta = mvBeta.GetSize();
		for (uint i = 0; i < myNBeta; i++)
			theDerivM.mdFx[theBegIndex + i] = theValue.mXt[theDate][i];
	}

	void cLinReg::ComputeGradF(uint theDate, const cRegArchValue& theValue, cDVector& theGradF)
	{
	uint myq = mvBeta.GetSize();
		for (uint i = 0; i < myq ; i++)
			theGradF[i] = theValue.mXt[theDate][i];
	}
	
	void cLinReg::ComputeHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData,cRegArchHessien& theHessData, uint theBegIndex, cAbstResiduals* theResiduals)
	{
	// Hess = 0 
	}

	void cLinReg::ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, const cDVector& theGradTheta, const cDVector& theGradU, double theGradH, uint theNu, uint theNh, cDMatrix& theHessTheta2,cDVector* theHessThetaU, cDVector& theHessThetaH, cDVector* theHessU2, cDVector& theHessUH, double& theHessH2)
	{
//		ComputeGradForM(theDate, theValue, theBegIndex, theGradTheta, theGradU, theGradH, theNu, theNh);

	}

	void cLinReg::ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, cFuncMeanAndVar& theDerivM)
	{

	}

	void cLinReg::ComputeHessF(uint theDate, const cRegArchValue& theData, cDMatrix& theHessF)
	{
		theHessF = 0.0;
	}

	void cLinReg::GetParamName(uint theIndex, char** theName)
	{
		uint myIndex = theIndex;
		for (uint i = 0; i < mvBeta.GetSize(); i++)
		{
			sprintf(theName[myIndex++], "LINREG[%d]", i + 1);

		}
	}

	void cLinReg::GetParamName(uint theIndex, string theName[])
	{
	uint myIndex = theIndex;
	char myChar[100] ;
		for (uint i = 0; i < mvBeta.GetSize(); i++)
		{
			sprintf(myChar, "LINREG[%d]", i + 1);
			theName[myIndex++] = myChar;

		}
	}

}//namespace
