#pragma once 
#ifndef _CVARINMEAN_H_
#define _CVARINMEAN_H_

#include "cAbstCondMean.h"
#include "cRegArchValue.h"
#include "cRegArchGradient.h"

/*!
 \file cVarInMean.h
 \brief header for class cVarInMean.

 \author Ollivier TARAMASCO
 \date dec-18-2006 - Last change feb-18-2011 void G
*/


namespace RegArchLib {
	class cVarInMean : public cAbstCondMean
	{
	private :
		double mvVarInMean ;
	public :
		cVarInMean(double theVarInMean = 0.0) ;
		cVarInMean(const cVarInMean& theVarInMean) ;
		~cVarInMean() ;
		void Delete(void) ;
		void Print(ostream& theOut) const ;
	#ifdef _RDLL_
		void Print(void);
	#endif //_RDLL_
		void SetDefaultInitPoint(double theMean, double theVar) ;
		void SetDefaultInitPoint(cRegArchValue& theValue);
		void Set(const double theValue, const uint theIndex = 0, const uint theNumParam = 0); ///< Set model parameters.
		void Set(const cDVector& theVectParam, const uint theNumParam=0) ; ///< Set model parameters.
		double Get(const uint theIndex, const uint theNumParam) ;
		cDVector& Get(const uint theNumParam);
		void ReAlloc(const uint theSize, const uint theNumParam=0) ;
		void ReAlloc(const cDVector& theVectParam, const uint theNumParam=0) ;
		cVarInMean& operator=(const cVarInMean& theSrc) ;
		void ReAllocProxyMeanParameters(uint theOldNParam=0) {};
		void UpdateProxyMeanParameters(void) {};
		double ComputeMean(uint theDate, const cRegArchValue& theData) const ;
		uint GetNParam(void) const ; ///< Return number of parameters
		uint GetNLags(void) const ;
		void RegArchParamToVector(cDVector& theDestVect, uint theIndex = 0) ;
		void VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex = 0) ;
		void ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, uint theBegIndex) ;
		void ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData,cRegArchHessien& theHessData, uint theBegIndex) ;
		void ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex);
		void GetParamName(uint theIndex, char** theName);
		void GetParamName(uint theIndex, string theName[]);
	} ;
}
#endif // _CVARINMEAN_H_
