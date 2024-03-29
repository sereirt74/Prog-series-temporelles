#include "StdAfxRegArchLib.h"
/*!
 \file cMa.cpp
 \brief sources for class cMa methods.

 \author Jean-Baptiste Durand, Ollivier TMAAMASCO
 \date dec-18-2006 - Last change feb-18-2011
*/
namespace RegArchLib {
	/*!
	 * \fn cMa::cMa(uint theNMa):cAbstCondMean(eMa)
	 * \param int theNMa: number of lags.
	 */

	cMa::cMa(uint theNMa):cAbstCondMean(eMa)
	{
		mvMa.ReAlloc(theNMa) ;
		MESS_CREAT("cMa")
	}

	/*!
	 * \fn cMa::cMa(const cDVector& theMa):cAbstCondMean(eMa)
	 * \param const cDVector& theAr: vector of AR coefficients.
	 */
	cMa::cMa(const cDVector& theMa):cAbstCondMean(eMa)
	{
		mvMa = theMa ;
		MESS_CREAT("cMa")
	}

	/*!
	 * \fn cMa::cMa(cAbstCondMean& theAbstCondMean)
	 * \param const cAbstCondMean& theAbstCondMean: the cMa source.
	 */
	cMa::cMa(const cMa& theMa):cAbstCondMean(eMa)
	{
		*this = theMa;
		MESS_CREAT("cMa")
	}

	/*!
	 * \fn cMa::~cMa()
	 */
	cMa::~cMa()
	{
		mvMa.Delete() ;
		MESS_DESTR("cMa")
	}

	/*!
	 * \fn cAbstCondMean cMa::PtrCopy(void)
	 * \param void
	 */
/*
cAbstCondMean* cMa::PtrCopy(void)
	{
		//	cArfima* myArfima = new cArfima(*this);

		//		return myArfima;
		return cAbstCondMeanPtrCopy<cMa>();
	}
*/
	/*!
	 * \fn void cMa::Delete(void)
	 * \param void
	 */
	void cMa::Delete(void)
	{
		mvMa.Delete() ;
	}

	/*!
	 * \fn void cMa::Print(ostream& theOut) const
	 * \param ostream& theOut: the output stream, default cout.
	 */
#ifndef _RDLL_
	void cMa::Print(ostream& theOut) const
	{
	uint myNMa = mvMa.GetSize();
		theOut << "MA(" << myNMa << ") mean with:" << endl ;
		for (uint i = 0 ; i <  myNMa ; i++)
			theOut << "\tMA[" << i+1 << "]=" << mvMa[i] << endl ;
	}
#else
	void cMa::Print(void)
	{
	uint myNMa = mvMa.GetSize();
		Rprintf("MA(%d) model with:", myNMa);
		for (uint i = 0; i < myNMa; i++)
			Rprintf("\tMA[%d]=%f", i + 1, mvMa[i]);
	}
#endif //_RDLL_

	void cMa::SetDefaultInitPoint(double theMean, double theVar)
	{
	uint myNMa = mvMa.GetSize() ;
		for (uint i = 0 ; i < myNMa ; i++)
			mvMa[i] = 0.0 ;
	}

	void cMa::SetDefaultInitPoint(cRegArchValue& theValue)
	{
		uint myNMa = mvMa.GetSize();
		for (uint i = 0; i < myNMa; i++)
			mvMa[i] = 0.0;
	}

	/*!
	 * \fn void cMa::Set(const double theValue, const uint theIndex=0, const uint theNumParam)
	 * \param const double theValue: the theIndex th value
	 * \param const uint theIndex: the index
	 * \param const uint theNumParam: not used here
	 * \details mvMa[theIndex] = theValue
	 */
	void cMa::Set(const double theValue, const uint theIndex, const uint theNumParam)
	{
		if (theIndex >= mvMa.GetSize())
			throw cError("Bad index") ;
		else
			mvMa[theIndex]=theValue ;
	}

	/*!
	 * \fn void cMa::Set(const cDVector& theVectParam, const uint theNumParam)
	 * \param const cDVector& theVectParam: the vector of MA coefficients
	 * \param const uint theNumParam: not used here
	 * \details mvMa = theVectParam
	 */
	void cMa::Set(const cDVector& theVectParam, const uint theNumParam)
	{
		mvMa=theVectParam ;
	}

	double cMa::Get(const uint theIndex, const uint theNumParam)
	{
		return mvMa[theIndex] ;
	}

	cDVector& cMa::Get(const uint theNumParam)
	{
		return mvMa;
	}

	/*!
	 * \fn void cMa::ReAlloc(const uint theSize, const uint theNumParam)
	 * \param const uint theSize: new size of mvMA
	 * \param const uint theNumParam; not used here.
	 * \details new allocation of mvMa 
	 */
	void cMa::ReAlloc(const uint theSize, const uint theNumParam)
	{
		mvMa.ReAlloc(theSize) ;
	}

	/*!
	 * \fn void cAr::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	 * \param const cDVector& theVectParam: the vector of AR coefficients
	 * \param const uint theNumParam: not used here.
	 * \details new allocation of mvAr
	 */
	void cMa::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{
		mvMa = theVectParam ;
	}

	/*!
	 * \fn cAbstCondMean& cAr::operator =(cAbstCondMean& theSrc)
	 * \param cAbstCondMean& theSrc: source to be recopied
	 * \details An error occurs if theSrc is not an cAr class parameter
	 */
	cMa& cMa::operator =(const cMa& theSrc)
	{
		mvMa = theSrc.mvMa ;
		return *this ;
	}

	/*!
	 * \fn cMa::ComputeMean(uint theDate, const cRegArchValue& theData) const
	 * \param int theDate: date of the computation
	 * \param cRegArchValue& theData: past datas.
	 * \details theData is not updated here.
	 */
	double cMa::ComputeMean(uint theDate, const cRegArchValue& theData) const
	{
	uint myq = mvMa.GetSize() ;

	double myRes = 0.0 ;
		for (uint i = 1 ; i <= MIN(myq, theDate) ; i++)
			myRes += mvMa[i-1] * theData.mUt[theDate-i] ;
		return myRes ;
	}

	uint cMa::GetNParam(void) const
	{
		return mvMa.GetSize() ;
	}
	
	uint cMa::GetNLags(void) const
	{
		return mvMa.GetSize() ;
	}

	uint cMa::GetNu(void) const
	{
		return mvMa.GetSize();
	}

	uint cMa::GetNh(void) const
	{
		return 0;
	}

	void cMa::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, uint theBegIndex, cAbstResiduals* theResids)
	{
	uint myq = mvMa.GetSize() ;
	uint i ;
		for (i = 1 ; i <= MIN(myq, theDate) ; i++)
			theGradData.mCurrentGradMu[theBegIndex+i-1] += theValue.mUt[theDate - i] ;
		for (i = 0 ; i < MIN(myq, theDate) ; i++)
			theGradData.mCurrentGradMu -=  mvMa[i] * theGradData.mGradMt[i] ;
	}

	void cMa::ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cDVector& theGradTheta, cDVector& theGradU, double& theGradH, uint& theNu, uint& theNh)
	{
		theNh = 0;
		theGradTheta = 0.0;
		theGradU = 0.0;
		theGradH = 0.0;
		uint myNMa = theNu = mvMa.GetSize();
		for (uint i = 1; i <= MIN(myNMa, theDate); i++)
		{
			theGradTheta[theBegIndex + i - 1] = theValue.mUt[theDate - i];
			theGradU[i - 1] = mvMa[i - 1];
		}
	}

	void cMa::ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cFuncMeanAndVar& theDerivM)
	{
		theDerivM.mdFx = 0.0;
		theDerivM.mdFu = 0.0;
		for (uint i = 1; i <= MIN(theDerivM.mNu, theDate); i++)
		{
			theDerivM.mdFx[theBegIndex + i - 1] = theValue.mUt[theDate - i];
			theDerivM.mdFu[i - 1] = mvMa[i - 1];
		}
	}

	void cMa::GetNParamF(uint theNParam[3]) const
	{
		theNParam[0] = theNParam[1] = mvMa.GetSize();
		theNParam[2] = 0;
	}

	void cMa::ComputeGradF(uint theDate, const cRegArchValue& theValue, cDVector& theGradF)
	{
		uint myq = mvMa.GetSize();
		uint i;
		for (i = 1; i <= MIN(myq, theDate); i++)
		{
			theGradF[i - 1] = theValue.mUt[theDate - i];
			theGradF[myq + i - 1] = mvMa[i-1];
		}
	}

	void cMa::ComputeHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData,cRegArchHessien& theHessData, uint theBegIndex, cAbstResiduals* theResiduals)
	{
	uint myq = mvMa.GetSize();
	uint myNParam = theGradData.GetNParam();
	cDMatrix myMat1(myNParam, myNParam);
	uint i ;
		for (i = 0; i < MIN(myq, theDate); i++)
		{	
			myMat1.SetColumn(theBegIndex + i, theGradData.mGradMt[i]);
		}
		theHessData.mCurrentHessMu -= myMat1 + Transpose(myMat1);
		for (i = 0; i < MIN(myq, theDate); i++)
			theHessData.mCurrentHessMu -= mvMa[i] * theHessData.mHessMt[i];
	}

	void cMa::ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, const cDVector& theGradTheta, const cDVector& theGradU, double theGradH, uint theNu, uint theNh, cDMatrix& theHessTheta2,cDVector* theHessThetaU, cDVector& theHessThetaH, cDVector* theHessU2, cDVector& theHessUH, double& theHessH2)
	{
		for (uint k = 0; k < theNu; k++)
			theHessThetaU[k][theBegIndex + k] = 1 ;
	}

	void cMa::ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, cFuncMeanAndVar& theDerivM)
	{
		theDerivM.md2Fxx = 0.0;
		theDerivM.md2Fxu = 0.0;
		theDerivM.md2Fxh = 0.0;
		theDerivM.md2Fuu = 0.0;
		theDerivM.md2Fuh = 0.0;
		theDerivM.md2Fhh = 0.0;

		for (uint k = 0; k < theDerivM.mNu; k++)
			theDerivM.md2Fxu[theBegIndex + k][k] = 1 ;

	}

	void cMa::ComputeHessF(uint theDate, const cRegArchValue& theData, cDMatrix& theHessF)
	{
	uint myq = mvMa.GetSize();
		theHessF = 0;
		for (uint i = 0; i < myq; i++)
			theHessF[i][myq + i] = theHessF[myq + i][i] = 1.0;
	}

	void cMa::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{
	uint mySize = mvMa.GetSize() ;
		if (theDestVect.GetSize() < mySize + theIndex)
			throw cError("Wrong size") ;
		mvMa.SetSubVectorWithThis(theDestVect, theIndex) ;
	}

	void cMa::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	uint mySize = theSrcVect.GetSize() ;
		if (mvMa.GetSize() + theIndex > mySize)
			throw cError("Wrong size") ;
		mvMa.SetThisWithSubVector(theSrcVect, theIndex) ;
	}

	void cMa::GetParamName(uint theIndex, char** theName)
	{
		uint myIndex = theIndex;
		for (uint i = 0; i < mvMa.GetSize(); i++)
		{
			sprintf(theName[myIndex++], "MA[%d]", i + 1);

		}
	}

	void cMa::GetParamName(uint theIndex, string theName[])
	{
	uint myIndex = theIndex;
	char myChar[100];
		for (uint i = 0; i < mvMa.GetSize(); i++)
		{
			sprintf(myChar, "MA[%d]", i + 1);
			theName[myIndex++] = myChar;

		}
	}

}//namespace
