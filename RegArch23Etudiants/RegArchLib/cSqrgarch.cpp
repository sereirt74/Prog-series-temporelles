#include "StdAfxRegArchLib.h"
/*!
 \file cSqrgarch.cpp
 \brief sources for class cStdDevInMean methods.

 \author Jean-Baptiste Durand, Ollivier TMAAMASCO
 \date feb-18-2024 - Last change feb-18-2024
*/
namespace RegArchLib {
	cSqrgarch::cSqrgarch(uint theNArch, uint theNGarch)
	:cAbstCondVar(eSqrgarch)  // call constructor of cAbstCondVar with type eTarch
	{
	}

	cSqrgarch::cSqrgarch(double theConst, cDVector& theArch, cDVector& theGarch, cDVector& theGamma)
	:cAbstCondVar(eSqrgarch)  // call constructor of cAbstCondVar with type eTarch
	{
	}

	cSqrgarch::cSqrgarch(const cSqrgarch& theSqrarch)
	:cAbstCondVar(eSqrgarch)  // call constructor of cAbstCondVar with type eTarch
	{
	}

	cSqrgarch::~cSqrgarch()
	{
	}

	/*!
	 * \fn cAbstCondVar* cSqrgarch::PtrCopy()
	 */
/*	cAbstCondVar* cSqrgarch::PtrCopy() const
	{
		//		 cConstCondVar *myConstCondVar = new cConstCondVar(*this);
		//		 return myConstCondVar;
		return cAbstCondVarPtrCopy<cSqrgarch>();
	}
*/
	void cSqrgarch::Delete(void)
	{
	}

#ifndef _RDLL_
	void cSqrgarch::Print(ostream& theOut) const
	{
	}
#else
	void cSqrgarch::Print(void)
	{
		Rprintf("TARCH(%d) model with:\n", mvArchPos.GetSize());
		Rprintf("Const=%f\n", mvCste);
		for (uint i = 0; i < mvArchPos.GetSize(); i++)
			Rprintf("ARCH+[%d]=%f\n", i+1,  mvArchPos[i]);
		for (uint i = 0; i < mvArchNeg.GetSize(); i++)
			Rprintf("ARCH-[%d]=%f\n", i+1,  mvArchNeg[i]);
	}


#endif // -RDLL-

	void cSqrgarch::SetDefaultInitPoint(double theMean, double theVar)
	{
	}

	void cSqrgarch::SetDefaultInitPoint(cRegArchValue& theValue)
	{
	}

	void cSqrgarch::ReAlloc(const uint theSize, const uint theNumParam)
	{
	}

	void cSqrgarch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{

	}

	void cSqrgarch::Set(const cDVector& thecDVector, const uint thePlace)
	{
	}

	void cSqrgarch::Set(const double theValue, const uint theIndex, const uint theNumParam)
	{
	}

	double cSqrgarch::Get(const uint theIndex, const uint theNumParam)
	{
		return 0;
	}

	cDVector& cSqrgarch::Get(const uint theNumParam)
	{
		return cDVector(0);
	}

	cSqrgarch& cSqrgarch::operator =(const cSqrgarch& theSrc)
	{
		return cSqrgarch(theSrc);
	}

	double cSqrgarch::ComputeVar(uint theDate, const cRegArchValue& theValue) const
	{
		return 0;
	}

	uint cSqrgarch::GetNParam(void) const
	{
		return 0;
	}

	uint cSqrgarch::GetNLags(void) const
	{
		return 0;
	}

	void cSqrgarch::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{

	}

	void cSqrgarch::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	}

	void cSqrgarch::ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cAbstResiduals* theResiduals)
	{

	}
	
	void cSqrgarch::ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{

	}

	void cSqrgarch::ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	}

	void cSqrgarch::GetParamName(uint theIndex, char** theName)
	{
	}

	void cSqrgarch::GetParamName(uint theIndex, string theName[])
	{
	}


}//namespace
