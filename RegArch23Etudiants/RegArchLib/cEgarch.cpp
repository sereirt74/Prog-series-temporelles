#include "StdAfxRegArchLib.h"
/*!
	\file cEgarch.cpp
	\brief sources for class cEgarch methods.

	\author Jean-Baptiste DURAND, Ollivier TARAMASCO 
	\date dec-18-2006 - Last change feb-18-2011
*/
namespace RegArchLib {

	cEgarch::cEgarch(int theNArch, int theNGarch)
	:cAbstCondVar(eEgarch)  // call constructor of cAbstCondVar with type eEgarch
	{
	}

	cEgarch::cEgarch(cAbstResiduals* theResiduals, int theNArch, int theNGarch)
	:cAbstCondVar(eEgarch)
	{
	}

	cEgarch::cEgarch(const cEgarch& theEgarch)
	:cAbstCondVar(eEgarch) 
	{
	}

	cEgarch::~cEgarch()
	{
	}

	void cEgarch::Delete(void)
	{
	}

	void cEgarch::ReAlloc(const uint theSize, const uint theNumParam)
	{
	}

	void cEgarch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{
	}

#ifndef _RDLL_
	void cEgarch::Print(ostream& theOut) const
	{
	}
#else
	void cEgarch::Print(void)
	{
		Rprintf("EGARCH(%d, %d) model with:\n", mvArch.GetSize(), mvGarch.GetSize());
		Rprintf("Const=%f\n", mvCste);
	uint i;
		for (i = 0; i < mvArch.GetSize(); i++)
			Rprintf("ARCH[%d]=%f\n", i + 1, mvArch[i]);
		for (i = 0; i < mvGarch.GetSize(); i++)
			Rprintf("GARCH[%d]=%f\n", i + 1, mvGarch[i]);
		Rprintf("Teta=%f\n", mvTeta);
		Rprintf("Gamma=%f\n", mvGamma);
	}
#endif // _RDLL_

	void cEgarch::SetDefaultInitPoint(double theMean, double theVar)
	{
	}

	void cEgarch::SetDefaultInitPoint(cRegArchValue& theValue)
	{
	}

	void cEgarch::Set(const cDVector& theVectParam, const uint theNumParam)
	{
	}

	void cEgarch::Set(const double theValue, const uint theIndex, const uint theNumParam) 
	{
	}

	double cEgarch::Get(const uint theIndex, const uint theNumParam)
	{
		return 0;
	}

	cDVector& cEgarch::Get(const uint theNumParam)
	{
		return cDVector(0);
	}

	cEgarch& cEgarch::operator =(const cEgarch& theSrc)
	{
		return cEgarch(theSrc);
	}

	double cEgarch::ComputeVar(uint theDate, const cRegArchValue& theValue) const
	{
		return 0;
	}

	uint cEgarch::GetNParam(void) const 
	{
		return 0;
	}

	uint cEgarch::GetNLags(void) const
	{
		return 0;
	}

	void cEgarch::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{
	}

	void cEgarch::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	}

	void cEgarch::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, cAbstResiduals* theResiduals)
	{
	}

	void cEgarch::SetEspAbsEps(double theEspAbsEps)
	{
	}

	void cEgarch::ComputeHess(uint theDate, const cRegArchValue& theValue, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	}

	
	void cEgarch::ComputeGradAndHess(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	}

	void cEgarch::GetParamName(uint theIndex, char** theName)
	{
	}
	
	void cEgarch::GetParamName(uint theIndex, string theName[])
	{
	}

}//namespace
