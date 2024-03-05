#include "StdAfxRegArchLib.h"
/*!
 \file cLoggarch.cpp
 \brief sources for class cLoggarch methods.

 \author Jean-Baptiste Durand, Ollivier TMAAMASCO
 \date feb-18-2024 - Last change feb-18-2024
*/
namespace RegArchLib {
	cLoggarch::cLoggarch(uint theNArch, uint theNGarch)
	:cAbstCondVar(eLoggarch)  // call constructor of cAbstCondVar with type eLogarch
	{
	}

	cLoggarch::cLoggarch(double theConst, cDVector& theArch, cDVector& theGarch)
	:cAbstCondVar(eLoggarch)  // call constructor of cAbstCondVar with type eLogarch
	{
	}

	cLoggarch::cLoggarch(const cLoggarch& theLoggarch)
	: cAbstCondVar(eLoggarch)  // call constructor of cAbstCondVar with type eLogarch
	{
	}

	cLoggarch::~cLoggarch()
	{
	}

	/*!
	 * \fn cAbstCondVar* cLoggarch::PtrCopy()
	 */
/*	cAbstCondVar* cLoggarch::PtrCopy() const
	{
		//		 cConstCondVar *myConstCondVar = new cConstCondVar(*this);
		//		 return myConstCondVar;
		return cAbstCondVarPtrCopy<cLoggarch>();
	}
*/
	void cLoggarch::Delete(void)
	{
	}

#ifndef _RDLL_
	void cLoggarch::Print(ostream& theOut) const
	{
	}
#else
	void cLoggarch::Print(void)
	{
		Rprintf("TARCH(%d) model with:\n", mvArchPos.GetSize());
		Rprintf("Const=%f\n", mvCste);
		for (uint i = 0; i < mvArchPos.GetSize(); i++)
			Rprintf("ARCH+[%d]=%f\n", i+1,  mvArchPos[i]);
		for (uint i = 0; i < mvArchNeg.GetSize(); i++)
			Rprintf("ARCH-[%d]=%f\n", i+1,  mvArchNeg[i]);
	}


#endif // -RDLL-

	void cLoggarch::SetDefaultInitPoint(double theMean, double theVar)
	{
	}

	void cLoggarch::SetDefaultInitPoint(cRegArchValue& theValue)
	{
	}

	void cLoggarch::ReAlloc(const uint theSize, const uint theNumParam)
	{
	}

	void cLoggarch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{

	}

	void cLoggarch::Set(const cDVector& thecDVector, const uint thePlace)
	{
	}

	void cLoggarch::Set(const double theValue, const uint theIndex, const uint theNumParam)
	{
	}

	double cLoggarch::Get(const uint theIndex, const uint theNumParam)
	{
		return 0;
	}

	cDVector& cLoggarch::Get(const uint theNumParam)
	{
		return cDVector(0);
	}

	cLoggarch& cLoggarch::operator =(const cLoggarch& theSrc)
	{
		return cLoggarch(theSrc);
	}

	double cLoggarch::ComputeVar(uint theDate, const cRegArchValue& theValue) const
	{
		return 0;
	}

	uint cLoggarch::GetNParam(void) const
	{
		return 0;
	}

	uint cLoggarch::GetNLags(void) const
	{
		return 0;
	}

	void cLoggarch::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{

	}

	void cLoggarch::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	}

	void cLoggarch::ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cAbstResiduals* theResiduals)
	{

	}
	
	void cLoggarch::ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{

	}

	void cLoggarch::ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	}

	void cLoggarch::GetParamName(uint theIndex, char** theName)
	{
	}

	void cLoggarch::GetParamName(uint theIndex, string theName[])
	{
	}


}//namespace
