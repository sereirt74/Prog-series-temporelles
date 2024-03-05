#include "StdAfxRegArchLib.h"
/*!
 \file cStdDevInMean.cpp
 \brief sources for class cStdDevInMean methods.

 \author Jean-Baptiste Durand, Ollivier TMAAMASCO
 \date dec-18-2006 - Last change feb-18-2011
*/
namespace RegArchLib {
	cStdDevInMean::cStdDevInMean(double theStdDevInMean)
	:cAbstCondMean(eStdDevInMean)
	{
	}

	/*!
	 * \fn cStdDevInMean::cStdDevInMean(cAbstCondMean& theAbstCondMean)
	 * \param const cAbstCondMean& theAbstCondMean: the cStdDevInMean source.
	 */
	cStdDevInMean::cStdDevInMean(const cStdDevInMean& theStdDevInMean)
	:cAbstCondMean(eStdDevInMean)
	{
	}

	cStdDevInMean::~cStdDevInMean()
	{
	}

	void cStdDevInMean::Delete(void)
	{
	}

#ifndef _RDLL_
	void cStdDevInMean::Print(ostream& theOut) const
	{
	}
#else
	void cStdDevInMean::Print(void)
	{
		Rprintf("Standard Deviation In Mean model with:%d\n") ;
		Rprintf("\tDelta = %f\n", mvStdDevInMean) ;
	}
#endif // _RDLL_


	void cStdDevInMean::SetDefaultInitPoint(double theMean, double theVar)
	{
	}

	void cStdDevInMean::SetDefaultInitPoint(cRegArchValue& theValue)
	{
	}

	void cStdDevInMean::Set(const cDVector& theVectParam, const uint theNumParam)
	{
	}

	void cStdDevInMean::Set(const double theValue, const uint theIndex, const uint theNumParam) 
	{
		mvStdDevInMean = theValue ;
	}

	double cStdDevInMean::Get(const uint theIndex, const uint theNumParam)
	{
		return 0;
	}

	cDVector& cStdDevInMean::Get(const uint theNumParam)
	{
		return cDVector(0);
	}

	void cStdDevInMean::ReAlloc(const uint theSize, const uint theNumParam)
	{
	}

	void cStdDevInMean::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{
	}

	cStdDevInMean& cStdDevInMean::operator =(const cStdDevInMean& theSrc)
	{
		return cStdDevInMean(theSrc);
	} 

	double cStdDevInMean::ComputeMean(uint theDate, const cRegArchValue& theData) const
	{
		return 0;
	}

	uint cStdDevInMean::GetNParam(void) const
	{
		return 1 ;
	}
	
	uint cStdDevInMean::GetNLags(void) const
	{
		return 0 ;
	}

	void cStdDevInMean::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{
	}

	void  cStdDevInMean::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	}

	void cStdDevInMean::ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, uint theBegIndex)
	{
	}

	void cStdDevInMean::ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData,cRegArchHessien& theHessData, uint theBegIndex)
	{
	}

	void cStdDevInMean::ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex)
	{
	}

	void cStdDevInMean::GetParamName(uint theIndex, char** theName)
	{
	}

	void cStdDevInMean::GetParamName(uint theIndex, string theName[])
	{
	}

}//namespace
