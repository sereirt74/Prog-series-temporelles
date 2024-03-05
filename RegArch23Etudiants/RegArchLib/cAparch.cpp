#include "StdAfxRegArchLib.h"

/*!
	\file cAparch.cpp
	\brief sources for abstract class cAparch methods.

	\author Jean-Baptiste DURAND, Ollivier TARAMASCO 
	\date dec-18-2006 - Last change feb-18-2011
*/

namespace RegArchLib {
	cAparch::cAparch(int theNArch, int theNGarch)
	:cAbstCondVar(eAparch)  // call constructor of cAbstCondVar with type eAparch
	{
	}

	cAparch::cAparch(const cAparch& theAparch)
	:cAbstCondVar(eAparch)  // call constructor of cAbstCondVar with type eAparch
	{
	}

	cAparch::~cAparch()
	{
	}

	/*!
	 * \fn cAbstCondVar* cAparch::PtrCopy()
	 */
/*
cAbstCondVar* cAparch::PtrCopy() const
	{

cAparch* myAparch = new cAparch(*this);

		 return myAparch;

		return cAbstCondVarPtrCopy<cAparch>();
	}
*/
	void cAparch::Delete(void)
	{
		mvArch.Delete() ;
		mvGamma.Delete() ;
		mvGarch.Delete() ;
	}

#ifdef _RDLL_
	void cAparch::Print(void)
	{
		Rprintf("APARCH(%d, %d) model with:\n", mvArch.GetSize(), mvGarch.GetSize());
		Rprintf("Const=%f\n", mvCste);
		Rprintf("Delta=%f\n", mvDelta);
		uint i;
		for (i = 0; i < mvArch.GetSize(); i++)
			Rprintf("Arch[%d]=%f\n", i + 1, mvArch[i]);
		for (i = 0; i < mvGamma.GetSize(); i++)
			Rprintf("Gamma[%d]=%f\n", i + 1, mvGamma[i]);
		for (i = 0; i < mvGarch.GetSize(); i++)
			Rprintf("Garch[d]=%f\n", i + 1, mvGarch[i]);
	}
#else
	void cAparch::Print(ostream& theOut) const
		{
		}
#endif // _RDLL_	

	void cAparch::SetDefaultInitPoint(double theMean, double theVar)
	{
	}

	void cAparch::SetDefaultInitPoint(cRegArchValue& theValue)
	{
	}

	void cAparch::ReAlloc(const uint theSize, const uint theNumParam)
	{
	}

	void cAparch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{
	}

	void cAparch::Set(const cDVector& theDVector, const uint thePlace) 
	{
	}

	void cAparch::Set(const double theValue, const uint theIndex, const uint theNumParam) 
	{
	}

	double cAparch::Get(const uint theIndex, const uint theNumParam)
	{
		return 0;
	}

	cDVector& cAparch::Get(const uint theNumParam)
	{
		return cDVector(0);
	}

	cAparch& cAparch::operator =(const cAparch& theSrc)
	{
		return cAparch(theSrc);
	}

	double cAparch::ComputeVar(uint theDate, const cRegArchValue& theValue) const
	{
		return 0;
	}

	uint cAparch::GetNParam(void) const
	{
		return 0;
	}

	uint cAparch::GetNLags(void) const
	{
		return 0;
	}

	void cAparch::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{
	}

	void cAparch::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	}

	void cAparch::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, cAbstResiduals* theResiduals)
	{
	}

	void cAparch::ComputeHess(uint theDate, const cRegArchValue& theValue, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	}

	void cAparch::ComputeGradAndHess(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{

	}

	void cAparch::GetParamName(uint theIndex, char** theName)
	{

	}


	void cAparch::GetParamName(uint theIndex, string theName[])
	{

	}

}//namespace
