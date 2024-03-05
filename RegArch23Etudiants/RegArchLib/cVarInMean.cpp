#include "StdAfxRegArchLib.h"
/*!
 \file cVarInMean.cpp
 \brief sources for class cVarInMean methods.

 \author Jean-Baptiste Durand, Ollivier TMAAMASCO
 \date dec-18-2006 - Last change feb-18-2011
*/

namespace RegArchLib {
	cVarInMean::cVarInMean(double theVarInMean)
		:cAbstCondMean(eVarInMean)
	{
	}

	/*!
	 * \fn cVarInMean::cVarInMean(cAbstCondMean& theAbstCondMean)
	 * \param const cAbstCondMean& theAbstCondMean: the cVarInMean source.
	 */
	cVarInMean::cVarInMean(const cVarInMean& theVarInMean)
		:cAbstCondMean(eVarInMean)
	{
	}

	cVarInMean::~cVarInMean()
	{
	}

	void cVarInMean::Delete(void)
	{
	}

#ifndef _RDLL_
	void cVarInMean::Print(ostream& theOut) const
	{
	}
#else
	void cVarInMean::Print(void)
	{
		Rprintf("Var In Mean model with:%d\n");
		Rprintf("\tDelta = %f\n", mvVarInMean);
	}
#endif // _RDLL_


	void cVarInMean::SetDefaultInitPoint(double theMean, double theVar)
	{
	}

	void cVarInMean::SetDefaultInitPoint(cRegArchValue& theValue)
	{
	}

	void cVarInMean::Set(const cDVector& theVectParam, const uint theNumParam)
	{
	}

	void cVarInMean::Set(const double theValue, const uint theIndex, const uint theNumParam)
	{
		mvVarInMean = theValue;
}

	double cVarInMean::Get(const uint theIndex, const uint theNumParam)
	{
		return 0;
	}

	cDVector& cVarInMean::Get(const uint theNumParam)
	{
		return cDVector(0);
	}

	void cVarInMean::ReAlloc(const uint theSize, const uint theNumParam)
	{
	}

	void cVarInMean::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{
	}

	cVarInMean& cVarInMean::operator =(const cVarInMean& theSrc)
	{
		return cVarInMean(theSrc);
	}

	double cVarInMean::ComputeMean(uint theDate, const cRegArchValue& theData) const
	{
		return 0;
	}

	uint cVarInMean::GetNParam(void) const
	{
		return 1;
	}

	uint cVarInMean::GetNLags(void) const
	{
		return 0;
	}

	void cVarInMean::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{
	}

	void  cVarInMean::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	}

	void cVarInMean::ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, uint theBegIndex)
	{
	}

	void cVarInMean::ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex)
	{
	}

	void cVarInMean::ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex)
	{
	}

	void cVarInMean::GetParamName(uint theIndex, char** theName)
	{
	}

	void cVarInMean::GetParamName(uint theIndex, string theName[])
	{
	}

}//namespace
