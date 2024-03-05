#include "StdAfxRegArchLib.h"
/*!
 \file cGtarch.cpp
 \brief sources for class cGtarch methods.

 \author Jean-Baptiste Durand, Ollivier TMAAMASCO
 \date feb-18-2024 - Last change feb-18-2024
*/
namespace RegArchLib {
	cGtarch::cGtarch(uint theNArch, uint theNGarch)
	:cAbstCondVar(eGtarch)  // call constructor of cAbstCondVar with type eGtarch
	{
	}

	cGtarch::cGtarch(double theConst, cDVector& theArchPos, cDVector& theArchNeg, cDVector& theGarch)
	:cAbstCondVar(eGtarch)  // call constructor of cAbstCondVar with type eGtarch
	{
	}

	cGtarch::cGtarch(const cGtarch& theTarch)
	:cAbstCondVar(eGtarch) // call constructor of cAbstCondVar with type eGtarch
	{
	}

	cGtarch::~cGtarch()
	{
	}

	/*!
	 * \fn cAbstCondVar* cGtarch::PtrCopy()
	 */
/*	cAbstCondVar* cGtarch::PtrCopy() const
	{
		//		 cConstCondVar *myConstCondVar = new cConstCondVar(*this);
		//		 return myConstCondVar;
		return cAbstCondVarPtrCopy<cGtarch>();
	}
*/
	void cGtarch::Delete(void)
	{
	}

#ifndef _RDLL_
	void cGtarch::Print(ostream& theOut) const
	{
	}
#else
	void cGtarch::Print(void)
	{
		Rprintf("TARCH(%d) model with:\n", mvArchPos.GetSize());
		Rprintf("Const=%f\n", mvCste);
		for (uint i = 0; i < mvArchPos.GetSize(); i++)
			Rprintf("ARCH+[%d]=%f\n", i+1,  mvArchPos[i]);
		for (uint i = 0; i < mvArchNeg.GetSize(); i++)
			Rprintf("ARCH-[%d]=%f\n", i+1,  mvArchNeg[i]);
	}


#endif // -RDLL-

	void cGtarch::SetDefaultInitPoint(double theMean, double theVar)
	{
	}

	void cGtarch::SetDefaultInitPoint(cRegArchValue& theValue)
	{
	}

	void cGtarch::ReAlloc(const uint theSize, const uint theNumParam)
	{
	}

	void cGtarch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{

	}

	void cGtarch::Set(const cDVector& thecDVector, const uint thePlace)
	{
	}

	void cGtarch::Set(const double theValue, const uint theIndex, const uint theNumParam)
	{
	}

	double cGtarch::Get(const uint theIndex, const uint theNumParam)
	{
		return 0;
	}

	cDVector& cGtarch::Get(const uint theNumParam)
	{
		return cDVector(0);
	}

	cGtarch& cGtarch::operator =(const cGtarch& theSrc)
	{
		return cGtarch(theSrc);
	}

	double cGtarch::ComputeVar(uint theDate, const cRegArchValue& theValue) const
	{
		return 0;
	}

	uint cGtarch::GetNParam(void) const
	{
		return 0;
	}

	uint cGtarch::GetNLags(void) const
	{
		return 0;
	}

	void cGtarch::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{

	}

	void cGtarch::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	}

	void cGtarch::ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cAbstResiduals* theResiduals)
	{

	}
	
	void cGtarch::ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{

	}

	void cGtarch::ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	}

	void cGtarch::GetParamName(uint theIndex, char** theName)
	{
	}

	void cGtarch::GetParamName(uint theIndex, string theName[])
	{
	}


}//namespace
