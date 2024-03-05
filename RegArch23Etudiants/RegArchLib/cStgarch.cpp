#include "StdAfxRegArchLib.h"
/*!
 \file cStgarch.cpp
 \brief sources for class cStdDevInMean methods.

 \author Jean-Baptiste Durand, Ollivier TMAAMASCO
 \date feb-18-2024 - Last change feb-18-2024
*/
namespace RegArchLib {
	cStgarch::cStgarch(uint theNArch, uint theNGarch)
	:cAbstCondVar(eStgarch)  // call constructor of cAbstCondVar with type eStgarch
	{
	}

	cStgarch::cStgarch(double theConst, cDVector& theArch, cDVector& theGarch, double theGamma)
	:cAbstCondVar(eStgarch)  // call constructor of cAbstCondVar with type eStgarch
	{
	}

	cStgarch::cStgarch(const cStgarch& theTarch)
	:cAbstCondVar(eStgarch)  // call constructor of cAbstCondVar with type eStgarch
	{
	}

	cStgarch::~cStgarch()
	{
	}

	/*!
	 * \fn cAbstCondVar* cStgarch::PtrCopy()
	 */
/*	cAbstCondVar* cStgarch::PtrCopy() const
	{
		//		 cConstCondVar *myConstCondVar = new cConstCondVar(*this);
		//		 return myConstCondVar;
		return cAbstCondVarPtrCopy<cStgarch>();
	}
*/
	void cStgarch::Delete(void)
	{
	}

#ifndef _RDLL_
	void cStgarch::Print(ostream& theOut) const
	{
	}
#else
	void cStgarch::Print(void)
	{
		Rprintf("TARCH(%d) model with:\n", mvArchPos.GetSize());
		Rprintf("Const=%f\n", mvCste);
		for (uint i = 0; i < mvArchPos.GetSize(); i++)
			Rprintf("ARCH+[%d]=%f\n", i+1,  mvArchPos[i]);
		for (uint i = 0; i < mvArchNeg.GetSize(); i++)
			Rprintf("ARCH-[%d]=%f\n", i+1,  mvArchNeg[i]);
	}


#endif // -RDLL-

	void cStgarch::SetDefaultInitPoint(double theMean, double theVar)
	{
	}

	void cStgarch::SetDefaultInitPoint(cRegArchValue& theValue)
	{
	}

	void cStgarch::ReAlloc(const uint theSize, const uint theNumParam)
	{
	}

	void cStgarch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{

	}

	void cStgarch::Set(const cDVector& thecDVector, const uint thePlace)
	{
	}

	void cStgarch::Set(const double theValue, const uint theIndex, const uint theNumParam)
	{
	}

	double cStgarch::Get(const uint theIndex, const uint theNumParam)
	{
		return 0;
	}

	cDVector& cStgarch::Get(const uint theNumParam)
	{
		return cDVector(0);
	}

	cStgarch& cStgarch::operator =(const cStgarch& theSrc)
	{
		return cStgarch(theSrc);
	}

	double cStgarch::ComputeVar(uint theDate, const cRegArchValue& theValue) const
	{
		return 0;
	}

	uint cStgarch::GetNParam(void) const
	{
		return 0;
	}

	uint cStgarch::GetNLags(void) const
	{
		return 0;
	}

	void cStgarch::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{

	}

	void cStgarch::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	}

	void cStgarch::ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cAbstResiduals* theResiduals)
	{

	}
	
	void cStgarch::ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{

	}

	void cStgarch::ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	}

	void cStgarch::GetParamName(uint theIndex, char** theName)
	{
	}

	void cStgarch::GetParamName(uint theIndex, string theName[])
	{
	}


}//namespace
