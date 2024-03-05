#include "StdAfxRegArchLib.h"
/*!
	\file cUgarch.cpp
	\brief sources for class cUgarch methods.

	\author Jean-Baptiste DURAND, Ollivier TARAMASCO 
	\date  sept-26-2016 - Last change feb-26-2016
*/
namespace RegArchLib {
	/*!
	 * \fn cUgarch::cUgarch(uint theNArch, uint theNGarch):cAbstCondVar(eUgarch)
     * \param bool theExistConst
	 * \param uint theNArch: number of ARCH lags
	 * \param uint theNGarch: number of GARCH lags
	*/
	cUgarch::cUgarch(bool theExistConst, uint theNBeta, uint theNArch, uint theNGarch)
	:cAbstCondVar(eUgarch)  // call constructor of cAbstCondVar with type eUgarch
	{
	}

	/*!
	 * \fn cUgarch::cUgarch(double theConst, cDVector& theArch, cDVector& theGarch):cAbstCondVar(eUgarch)
	 * \param bool theExistConst: true if const exists
	 * \param double theConst: constant part of the UGARCH(p, q) model if exists
	 * \param cDVector& theBetaS : Beta parameters
	 * \param cDVector& theGarch theArch: ARCH parameters
	 * \param cDVector& theGarch theGarch: GARCH parameters
	*/
	cUgarch::cUgarch(bool theExistConst, double theConst, cDVector& theBetaS, cDVector& theArch, cDVector& theGarch)
	:cAbstCondVar(eUgarch)
	{
	}

	/*!
	* \fn cUgarch::cUgarch(const cAbsCondVar& theGarch):cAbstCondVar(eUgarch)
	* \param cAbsCondVar& theEgarch: theGarch class
	*/
	cUgarch::cUgarch(const cUgarch& theUgarch)
	:cAbstCondVar(eUgarch)
	{
	}

	/*!
	 * \fn cUgarch::~cUgarch()
	*/
	cUgarch::~cUgarch()
	{
	}

	/*!
	 * \fn void cUgarch::Delete(void)
	 * \param void
	 * \details Free memory
	 */
	void cUgarch::Delete(void)
	{
	}

	/*!
	 * \fn void cUgarch::Print(ostream& theOut) const
	 * \param ostream& theOut: the output stream, default cout.
	 */
#ifndef _RDLL_
	void cUgarch::Print(ostream& theOut) const
	{
	}
#else
	void cUgarch::Print(void)
	{
	uint myNArch = mvArch.GetSize();
	uint myNUgarch = mvGarch.GetSize();
    uint myNBeta = mvBetaS.GetSize() ;
    	Rprintf("UGARCH(%d, %d) model with:\n", myNArch, myNUgarch);
		if (mvExistConst)
        {   Rprintf("\tCste=%f\n", mvConst);
		}
        for (uint i = 0; i < myNBeta; i++)
			Rprintf("\tBeta[%d]=%f\n", i + 1, mvBetaS[i]);
		for (uint i = 0; i < myNArch; i++)
			Rprintf("\tARCH[%d]=%f\n", i + 1, mvArch[i]);
		for (uint j = 0; j < myNUgarch; j++)
			Rprintf("\tGARCH[%d]=%f\n", j + 1, mvGarch[j]);
	}
#endif //_RDLL_

	void cUgarch::SetDefaultInitPoint(double theMean, double theVar)
	{
	}

	void cUgarch::SetDefaultInitPoint(cRegArchValue& theValue)
	{
	}

	/*!
	 * \fn void cUgarch::ReAlloc(const uint theSize, const uint theNumParam)
	 * \param const uint theSize: new size of mvArch or mvGarch
	 * \param const uint theNumParam: 0 for mvArch, 1 for mvGarch.
	 * \details new allocation of mvArch or mvGarch
	 */
	void cUgarch::ReAlloc(const uint theSize, const uint theNumParam)
	{
	}

	/*!
	 * \fn void cUgarch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	 * \param const cDVector& theVectParam: the vector of Const, ARCH or GARCH coefficients
	 * \param const uint theNumParam: =0, the constant part; =1 the ARCH coefficients; =2 theGARCH Coefficients
	 * \details new allocation of mvArch or mvConst
	 */
	void cUgarch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{
	}

	/*!
	 * \fn void cUgarch::Set(const double theValue, const uint theIndex, const uint theNumParam)
	 * \brief fill the parameters vector
	 * \param const double theValue: the value of the "theIndex" th lag. Default 0.
	 * \param const uint theIndex: the index.
	 * \param const uint theNumParam: =0, mvConst, =1, ARCH parameters; =2, GARCH parameters
	 * \details mvArch[theIndex] = theValue or mvGarch[theIndex]= theValue or mvConst = theValue
	 */
	void cUgarch::Set(const double theValue, const uint theIndex, const uint theNumParam)
	{			
	}


	/*!
	 * \fn void cUgarch::Set(const cDVector& theVectParam, const uint theNumParam)
	 * \brief fill the parameters vector
	 * \param const cDVector& theVectParam: the vector of values
	 * \param const uint theNumParam: =0, mvConst, =1, ARCH parameters; =2, GARCH parameters
	 * \details mvAr = theValue
	 */
	void cUgarch::Set(const cDVector& theVectParam, const uint theNumParam)
	{
	}

	double  cUgarch::Get(const uint theIndex, const uint theNumParam)
	{
		return 0;

	}

	cDVector& cUgarch::Get(const uint theNumParam)
	{
		return cDVector(0);
	
	}

	cUgarch& cUgarch::operator =(const cUgarch& theSrc)
	{
		return cUgarch(theSrc);
	}

	/*!
	 * \fn double cUgarch::ComputeVar(uint theDate, const cRegArchValue& theData) const
	 * \param int theDate: date of computation
	 * \param const cRegArchValue& theData: past datas
	 * \details theData is not updated here.
	*/
	double cUgarch::ComputeVar(uint theDate, const cRegArchValue& theData) const 
	{
		return 0;

	}

	uint cUgarch::GetNParam(void) const
	{
		return 0;
	}

	uint cUgarch::GetNLags(void) const
	{
		return 0;
	}

	void cUgarch::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, cAbstResiduals* theResiduals)
	{
	}

	void cUgarch::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{
	}

	void cUgarch::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	}

	void cUgarch::ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	}
	
	void cUgarch::ComputeGradAndHess(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	}

	void cUgarch::GetParamName(uint theIndex, char** theName)
	{
	}

	void cUgarch::GetParamName(uint theIndex, string theName[])
	{
	}


}//namespace
