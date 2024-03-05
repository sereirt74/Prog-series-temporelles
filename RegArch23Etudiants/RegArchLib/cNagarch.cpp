#include "StdAfxRegArchLib.h"
/*!
 \file cNagarch.cpp
 \brief sources for class cStdDevInMean methods.

 \author Jean-Baptiste Durand, Ollivier TMAAMASCO
 \date feb-18-2024 - Last change feb-18-2024
*/
namespace RegArchLib {
	cNagarch::cNagarch(uint theNArch, uint theNGarch)
	:cAbstCondVar(eNagarch)  // call constructor of cAbstCondVar with type eTarch
	{
	}

	cNagarch::cNagarch(double theConst, cDVector& theArch, cDVector& theGarch, cDVector& theGamma)
	:cAbstCondVar(eNagarch)  // call constructor of cAbstCondVar with type eTarch
	{
	}

	cNagarch::cNagarch(const cNagarch& theNagarch)
	:cAbstCondVar(eNagarch)  // call constructor of cAbstCondVar with type eTarch
	{
	}

	cNagarch::~cNagarch()
	{
	}

	/*!
	 * \fn cAbstCondVar* cNagarch::PtrCopy()
	 */
/*	cAbstCondVar* cNagarch::PtrCopy() const
	{
		//		 cConstCondVar *myConstCondVar = new cConstCondVar(*this);
		//		 return myConstCondVar;
		return cAbstCondVarPtrCopy<cNagarch>();
	}
*/
	void cNagarch::Delete(void)
	{
	}

#ifndef _RDLL_
	void cNagarch::Print(ostream& theOut) const
	{
	}
#else
	void cNagarch::Print(void)
	{
		Rprintf("TARCH(%d) model with:\n", mvArchPos.GetSize());
		Rprintf("Const=%f\n", mvCste);
		for (uint i = 0; i < mvArchPos.GetSize(); i++)
			Rprintf("ARCH+[%d]=%f\n", i+1,  mvArchPos[i]);
		for (uint i = 0; i < mvArchNeg.GetSize(); i++)
			Rprintf("ARCH-[%d]=%f\n", i+1,  mvArchNeg[i]);
	}


#endif // -RDLL-

	void cNagarch::SetDefaultInitPoint(double theMean, double theVar)
	{
	}

	void cNagarch::SetDefaultInitPoint(cRegArchValue& theValue)
	{
	}

	void cNagarch::ReAlloc(const uint theSize, const uint theNumParam)
	{
	}

	void cNagarch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{

	}

	void cNagarch::Set(const cDVector& thecDVector, const uint thePlace)
	{
	}

	void cNagarch::Set(const double theValue, const uint theIndex, const uint theNumParam)
	{
	}

	double cNagarch::Get(const uint theIndex, const uint theNumParam)
	{
		return 0;
	}

	cDVector& cNagarch::Get(const uint theNumParam)
	{
		return cDVector(0);
	}

	cNagarch& cNagarch::operator =(const cNagarch& theSrc)
	{
		return cNagarch(theSrc);
	}

	double cNagarch::ComputeVar(uint theDate, const cRegArchValue& theValue) const
	{
		return 0;
	}

	uint cNagarch::GetNParam(void) const
	{
		return 0;
	}

	uint cNagarch::GetNLags(void) const
	{
		return 0;
	}

	void cNagarch::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{

	}

	void cNagarch::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	}

	void cNagarch::ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cAbstResiduals* theResiduals)
	{

	}
	
	void cNagarch::ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{

	}

	void cNagarch::ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	}

	void cNagarch::GetParamName(uint theIndex, char** theName)
	{
	}

	void cNagarch::GetParamName(uint theIndex, string theName[])
	{
	}


}//namespace
