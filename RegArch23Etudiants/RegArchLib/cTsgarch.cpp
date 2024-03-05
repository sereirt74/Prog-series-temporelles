#include "StdAfxRegArchLib.h"
/*!
 \file cTsgarch.cpp
 \brief sources for class cTsgarch methods.

 \author Jean-Baptiste Durand, Ollivier TMAAMASCO
 \date feb-18-2024 - Last change feb-18-2024
*/

namespace RegArchLib {
	cTsgarch::cTsgarch(uint theNArch, uint theNGarch)
	:cAbstCondVar(eTsgarch)  // call constructor of cAbstCondVar with type eTsgarch
	{
	}

	cTsgarch::cTsgarch(double theConst, cDVector& theArch, cDVector& theGarch)
	:cAbstCondVar(eTsgarch)  // call constructor of cAbstCondVar with type eTsgarch
	{
	}

	cTsgarch::cTsgarch(const cTsgarch& theTsgarch)
	: cAbstCondVar(eTsgarch)  // call constructor of cAbstCondVar with type eTsgarch
	{
	}

	cTsgarch::~cTsgarch()
	{
	}

	/*!
	 * \fn cAbstCondVar* cTsgarch::PtrCopy()
	 */
/*	cAbstCondVar* cTsgarch::PtrCopy() const
	{
		//		 cConstCondVar *myConstCondVar = new cConstCondVar(*this);
		//		 return myConstCondVar;
		return cAbstCondVarPtrCopy<cTsgarch>();
	}
*/
	void cTsgarch::Delete(void)
	{
	}

#ifndef _RDLL_
	void cTsgarch::Print(ostream& theOut) const
	{
	}
#else
	void cTsgarch::Print(void)
	{
		Rprintf("TARCH(%d) model with:\n", mvArchPos.GetSize());
		Rprintf("Const=%f\n", mvCste);
		for (uint i = 0; i < mvArchPos.GetSize(); i++)
			Rprintf("ARCH+[%d]=%f\n", i+1,  mvArchPos[i]);
		for (uint i = 0; i < mvArchNeg.GetSize(); i++)
			Rprintf("ARCH-[%d]=%f\n", i+1,  mvArchNeg[i]);
	}


#endif // -RDLL-

	void cTsgarch::SetDefaultInitPoint(double theMean, double theVar)
	{
	}

	void cTsgarch::SetDefaultInitPoint(cRegArchValue& theValue)
	{
	}

	void cTsgarch::ReAlloc(const uint theSize, const uint theNumParam)
	{
	}

	void cTsgarch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{

	}

	void cTsgarch::Set(const cDVector& thecDVector, const uint thePlace)
	{
	}

	void cTsgarch::Set(const double theValue, const uint theIndex, const uint theNumParam)
	{
	}

	double cTsgarch::Get(const uint theIndex, const uint theNumParam)
	{
		return 0;
	}

	cDVector& cTsgarch::Get(const uint theNumParam)
	{
		return cDVector(0);
	}

	cTsgarch& cTsgarch::operator =(const cTsgarch& theSrc)
	{
		return cTsgarch(theSrc);
	}

	double cTsgarch::ComputeVar(uint theDate, const cRegArchValue& theValue) const
	{
		return 0;
	}

	uint cTsgarch::GetNParam(void) const
	{
		return 0;
	}

	uint cTsgarch::GetNLags(void) const
	{
		return 0;
	}

	void cTsgarch::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{

	}

	void cTsgarch::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	}

	void cTsgarch::ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cAbstResiduals* theResiduals)
	{

	}
	
	void cTsgarch::ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{

	}

	void cTsgarch::ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	}

	void cTsgarch::GetParamName(uint theIndex, char** theName)
	{
	}

	void cTsgarch::GetParamName(uint theIndex, string theName[])
	{
	}


}//namespace
