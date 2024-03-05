#include "StdAfxRegArchLib.h"
/*!
	\file cNgarch.cpp
	\brief sources for class cNgarch methods.

	\author Jean-Baptiste DURAND, Ollivier TARAMASCO 
	\date dec-18-2006 - Last change feb-18-2011
*/
namespace RegArchLib {
	/*!
	 * \fn cNgarch::cNgarch(uint theNArch, uint theNGarch):cAbstCondVar(eNgarch)
	 * \param uint theNArch: number of ARCH lags
	 * \param uint theNGarch: number of GARCH lags
	*/
	cNgarch::cNgarch(uint theNArch, uint theNGarch)
		:cAbstCondVar(eNgarch)  // call constructor of cAbstCondVar with type eGarch
	{
	}

	/*!
	 * \fn cNgarch::cNgarch(double theConst, cDVector& theArch, cDVector& theGarch):cAbstCondVar(eGarch)
	 * \param double theConst: constant part of the GARCH(p, q) model
	 * \param cDVector& theGarch theArch: ARCH parameters
	 * \param cDVector& theGarch theGarch: GARCH parameters
	*/
	cNgarch::cNgarch(double theConst, double theTheta, cDVector& theArch, cDVector& theGarch)
		:cAbstCondVar(eNgarch)
	{
	}

	/*!
	* \fn cNgarch::cNgarch(const cAbsCondVar& theGarch):cAbstCondVar(eGarch)
	* \param cAbsCondVar& theEgarch: theGarch class
	*/
	cNgarch::cNgarch(const cNgarch& theNgarch)
		:cAbstCondVar(eNgarch)
	{
	}

	/*!
	 * \fn cNgarch::~cNgarch()
	*/
	cNgarch::~cNgarch()
	{
	}

	/*!
	 * \fn cAbstCondVar* cNgarch::PtrCopy()
	 */
/*	cAbstCondVar* cNgarch::PtrCopy() const
	{
		//		 cConstCondVar *myConstCondVar = new cConstCondVar(*this);
		//		 return myConstCondVar;
		return cAbstCondVarPtrCopy<cNgarch>();
	}
*/
	/*!
	 * \fn void cNgarch::Delete(void)
	 * \param void
	 * \details Free memory
	 */
	void cNgarch::Delete(void)
	{
	}

	/*!
	 * \fn void cNgarch::Print(ostream& theOut) const
	 * \param ostream& theOut: the output stream, default cout.
	 */
#ifndef _RDLL_
	void cNgarch::Print(ostream& theOut) const
	{
	}
#else
	void cNgarch::Print(void)
	{
	uint myNArch = mvArch.GetSize();
	uint myNGarch = mvGarch.GetSize();
		Rprintf("NGARCH(%d, %d) model with:\n", myNArch, myNGarch);
		Rprintf("\tCste=%f\n", mvConst);
		Rprintf("\tTheta=%f\n", mvTheta);
		for (uint i = 0; i < myNArch; i++)
			Rprintf("\tARCH[%d]=%f\n", i + 1, mvArch[i]);
		for (uint j = 0; j < myNGarch; j++)
			Rprintf("\tGARCH[%d]=%f\n", j + 1, mvGarch[j]);
	}
#endif //_RDLL_

	void cNgarch::SetDefaultInitPoint(double theMean, double theVar)
	{
	}

	void cNgarch::SetDefaultInitPoint(cRegArchValue& theValue)
	{
	}

	/*!
	 * \fn void cNgarch::ReAlloc(const uint theSize, const uint theNumParam)
	 * \param const uint theSize: new size of mvArch or mvGarch
	 * \param const uint theNumParam: 0 for mvArch, 1 for mvGarch.
	 * \details new allocation of mvArch or mvGarch
	 */
	void cNgarch::ReAlloc(const uint theSize, const uint theNumParam)
	{
	}

	/*!
	 * \fn void cNgarch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	 * \param const cDVector& theVectParam: the vector of Const, ARCH or GARCH coefficients
	 * \param const uint theNumParam: =0, the constant part; =1 the ARCH coefficients; =2 theGARCH Coefficients
	 * \details new allocation of mvArch or mvConst
	 */
	void cNgarch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{
	}

	/*!
	 * \fn void cNgarch::Set(const double theValue, const uint theIndex, const uint theNumParam)
	 * \brief fill the parameters vector
	 * \param const double theValue: the value of the "theIndex" th lag. Default 0.
	 * \param const uint theIndex: the index.
	 * \param const uint theNumParam: =0, mvConst, =1, ARCH parameters; =2, GARCH parameters
	 * \details mvArch[theIndex] = theValue or mvGarch[theIndex]= theValue or mvConst = theValue
	 */
	void cNgarch::Set(const double theValue, const uint theIndex, const uint theNumParam)
	{
	}

	/*!
	 * \fn void cNgarch::Set(const cDVector& theVectParam, const uint theNumParam)
	 * \brief fill the parameters vector
	 * \param const cDVector& theVectParam: the vector of values
	 * \param const uint theNumParam: =0, mvConst, =1, ARCH parameters; =2, GARCH parameters
	 * \details mvAr = theValue
	 */
	void cNgarch::Set(const cDVector& theVectParam, const uint theNumParam)
	{	
	}

	double  cNgarch::Get(const uint theIndex, const uint theNumParam)
	{
		return 0;
	}

	cDVector& cNgarch::Get(const uint theNumParam)
	{
		return cDVector(0);
	}

	cNgarch& cNgarch::operator =(const cNgarch& theSrc)
	{
		return cNgarch(theSrc);
	}

	/*!
	 * \fn double cNgarch::ComputeVar(uint theDate, const cRegArchValue& theData) const
	 * \param int theDate: date of computation
	 * \param const cRegArchValue& theData: past datas
	 * \details theData is not updated here.
	*/
	double cNgarch::ComputeVar(uint theDate, const cRegArchValue& theData) const 
	{
		return 0;
	}

	uint cNgarch::GetNParam(void) const
	{
		return 0;
	}

	uint cNgarch::GetNLags(void) const
	{
		return 0;
	}

	void cNgarch::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, cAbstResiduals* theResiduals)
	{
	}

	void cNgarch::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{
	}

	void cNgarch::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	}

	void cNgarch::ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	}

	void cNgarch::ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	}

	void cNgarch::GetParamName(uint theIndex, char** theName)
	{
	}

	void cNgarch::GetParamName(uint theIndex, string theName[])
	{
	}

}//namespace
