#include "StdAfxRegArchLib.h"
/*!
	\file cAbstCondMean.cpp
	\brief sources for abstract class cConst methods.

	\author Jean-Baptiste Durand, Ollivier TARAMASCO
	\date dec-18-2006 - Last change feb-18-2011
*/
namespace RegArchLib {
	/*!
	 * \fn cConst::cConst(double theVal):cAbstCondMean(eConst)
	 * \param double theVal: constant value, default 0.0L.
	 */
	cConst::cConst(double theVal):cAbstCondMean(eConst)
	{	mvConst = theVal ;
		MESS_CREAT("cConst") ;
	}

	/*!
	 * \fn cConst::cConst(cAbstCondMean& theAbstCondMean)
	 * \param const cAbstCondMean& theAbstCondMean: the cConst source.
	 */
	cConst::cConst(const cConst& theConst):cAbstCondMean(eConst)
	{
		*this = theConst;
		MESS_CREAT("cConst") ;
	}

	/*!
	 * \fn cConst::~cConst()
	 * \details Nothing to do.
	 */
	cConst::~cConst()
	{	MESS_DESTR("cConst") ;
	}

	/*!
	 * \fn cAbstCondMean cConst::PtrCopy(void)
	 * \param void
	 */
/*
cAbstCondMean* cConst::PtrCopy(void)
 	{
//	cConst* myConst = new cConst(*this) ;

//		return myConst ;
	
		return cAbstCondMeanPtrCopy<cConst>();
	}
*/

	/*!
	 * \fn void cConst::Delete(void)
	 * \\details Delete. Nothing to do.
	 */
	void cConst::Delete(void)
	{	MESS_DESTR("cConst") ;
	}

	/*!
	 * \fn void cConst::Print(ostream& theOut) const
	 * \param ostream& theOut: output stream (screen or file). Default cout.
	 */

#ifndef _RDLL_
	void cConst::Print(ostream& theOut) const
	{
		theOut << "Constant mean with:" << endl;
		theOut << "\tConstant=" << mvConst << endl ;
	}

#else
	void cConst::Print(void)
	{
		Rprintf("Constant mean with:\n");
		Rprintf("\tConstant=%f\n", mvConst);
	}
#endif // _RDLL_


	void cConst::SetDefaultInitPoint(double theMean, double theVar)
	{
		mvConst = theMean ;
	}

	void cConst::SetDefaultInitPoint(cRegArchValue& theValue)
	{
	double myMean, myVar;
		theValue.ComputeMeanAndVar(myMean, myVar);
		mvConst = myMean;
	}

	/*!
	 * \fn void cConst::Set(const double theValue, const uint theIndex, const uint theNumParam)
	 * \brief fill the parameters vector
	 * \param const double theValue: the constant value.
	 * \param const uint theIndex: not used here. Default 0.
	 * \param const uint theNumParam: not used for cConst model. Default 0.
	 * \details mvConst = theValue
	 */
	void cConst::Set(const double theValue, const uint theIndex, const uint theNumParam)
	{	mvConst = theValue ;
	}

	/*!
	 * \fn void cConst::Set(const cDVector& theVectParam, const uint theNumParam)
	 * \brief fill the parameters vector
	 * \param const cDVector& theVectParam: the constant value is in theVectParam[0].
	 * \param const uint theIndex: not used here. Default 0.
	 * \param const uint theNumParam: not used for cConst model. Default 0.
	 * \details mvConst = theVectParam[0]
	 */

	void cConst::Set(const cDVector& theVectParam, const uint theNumParam)
	{	if (theVectParam.GetSize() > 0)
			mvConst = theVectParam[0] ;
		else
			throw cError("the size of theVectParam must be > 0") ;
	}

	double  cConst::Get(const uint theIndex, const uint theNumParam)
	{
		return mvConst ;
	}

	cDVector& cConst::Get(const uint theNumParam)
	{
	cDVector* myAux = new cDVector(1, mvConst);
		return *myAux;
	}

	/*!
	 * \fn void cConst::ReAlloc(const uint theSize, const uint theNumParam=0)
	 * \param const uint theSize: not used. Not used for cConstClass
	 * \param const uint theNumParam: not used for cConst class
	 * \details Nothing to do for cConst Class.
	 */
	void cConst::ReAlloc(const uint theSize, const uint theNumParam)
	{
	}

	/*!
	 * \fn void cConst::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	 * \param const cDVector& theVectParam: the constant value is in theVectParam[0]
	 * \param const uint theNumParam: not used for cConst class
	 * \details Here, mvConst = theVectParam[0]
	 */
	void cConst::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{
		if (theVectParam.GetSize() > 0)
			mvConst = theVectParam[0] ;
		else
			throw cError("Size of 'theVectParam' must be > 0") ;
	}

	/*!
	 * \fn cAbstCondMean& cConst::operator =(cAbstCondMean &theSrc)
	 * \param cAbstCondMean &theSrc
	 * \details theSrc must be a cConst parameter
	 */
	cConst& cConst::operator =(const cConst& theSrc)
	{
		if (GetCondMeanType() == eUnknown)
		{
			SetCondMeanType(eConst);
		}
		mvConst = theSrc.mvConst;
		return *this ;
	}

	/*!
	 * \param int theDate: date of the computation
	 * \param cRegArchValue& theData: past datas.
	 * \details theData must not be updated here.
	 */
	double cConst::ComputeMean(uint theDate, const cRegArchValue& theData) const
	{
		return mvConst ;
	}

	uint cConst::GetNParam(void) const
	{	return 1 ;
	}

	uint cConst::GetNu(void) const
	{
		return 0;
	}

	uint cConst::GetNh(void) const
	{
		return 0;
	}

	uint cConst::GetNLags(void) const
	{	return 0 ;
	}

	void cConst::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, uint theBegIndex, cAbstResiduals* theResids)
	{
		theGradData.mCurrentGradMu[theBegIndex] = 1.0 ;
	}

	void cConst::GetNParamF(uint theNParam[3]) const
	{
		theNParam[0] =  1;
		theNParam[1] = theNParam[2] = 0;
	}

	void cConst::ComputeGradF(uint theDate, const cRegArchValue& theValue, cDVector& theGradF)
	{
		theGradF[0] = 1;
	}

	void cConst::ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cDVector& theGradTheta, cDVector& theGradU, double& theGradH, uint& theNu, uint& theNh)
	{
		theNu = theNh = 0;
		theGradTheta = 0.0;
		theGradH = 0.0;
		theGradTheta[theBegIndex] = 1;
	}

	void cConst::ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cFuncMeanAndVar& theDerivM)
	{
		theDerivM.mdFx = 0.0;
		theDerivM.mdFx[theBegIndex] = 1;
	}

	void cConst::ComputeHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex, cAbstResiduals* theResiduals)
	{
	}

	void cConst::ComputeHessF(uint theDate, const cRegArchValue& theData, cDMatrix& theHessF)
	{
		theHessF[0][0] = 0;
	}

	void cConst::ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, const cDVector& theGradTheta, const cDVector& theGradU, double theGradH, uint theNu, uint theNh, cDMatrix& theHessTheta2,cDVector* theHessThetaU, cDVector& theHessThetaH, cDVector* theHessU2, cDVector& theHessUH, double& theHessH2)
	{
	}

	void cConst::ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, cFuncMeanAndVar& theDerivM)
	{
		theDerivM.md2Fxx = 0.0;

	}
	
	void cConst::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{
		if ((int)theDestVect.GetSize() + 1 < (int)theIndex)
			throw cError("Wrong size") ;
		theDestVect[theIndex] = mvConst ;
	}
	
	void cConst::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
		if (1 + theIndex > theSrcVect.GetSize())
			throw cError("Wrong size") ;
		mvConst = theSrcVect[theIndex] ;
	}

	void cConst::GetParamName(uint theIndex, char** theName)
	{
		uint myIndex = theIndex;
		sprintf(theName[myIndex++], "CONST");
	}

	void cConst::GetParamName(uint theIndex, string theName[])
	{
	uint myIndex = theIndex;
	char myChar[100];
		sprintf(myChar, "CONST");
		theName[myIndex++] = myChar;
	}

}//namespace
