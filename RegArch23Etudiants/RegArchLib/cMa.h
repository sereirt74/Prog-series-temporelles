#pragma once 
#ifndef _CMA_H_
#define _CMA_H_

#include "cAbstCondMean.h"
#include "cRegArchValue.h"
#include "cRegArchGradient.h"
#include "cRegArchHessien.h"

/*!
	\file cMa.h
	\brief header for class cMa.
	\author Jean-Baptiste DURAND, Ollivier TARAMASCO
	\date dec-18-2006 - Last change feb-18-2011
*/

/*
	mm(t, Teta) = F(Teta, ui, hj) = somme(i=1 à q, Teta(i) ui) où Teta(i) = MA(i)
	gradtetaF = (u1, ..., uq) 
	hesstetaF = 0 (qxq)
*/
namespace RegArchLib {
	/*! 
	 * \class cMa
	 * \brief  Class to implement an MA model
	 */
	class cMa: public cAbstCondMean
	{
	private :
		cDVector mvMa ;
	public :
		cMa(uint theNMa = 0) ; ///< A simple constructor
		cMa(const cDVector& theMa) ; ///< Another constructor
		cMa(const cMa& theMa) ; /// Recopy constructor
		virtual ~cMa() ;  ///< A simple destructor
		void Delete(void) ; ///< Free memory
	#ifndef _RDLL_
		void Print(ostream& theOut = cout) const; ///< print the parameters
	#else
		void Print(void);
	#endif //_RDLL_
		void SetDefaultInitPoint(double theMean, double theVar) ;
		void SetDefaultInitPoint(cRegArchValue& theValue);
		void ReAlloc(const uint theSize, const uint theNumParam = 0); ///< Allocation of the model parameters
		void ReAlloc(const cDVector& theVectParam, const uint theNumParam=0) ; ///< Allocation of the model parameters
		void Set(const double theValue, const uint theIndex=0, const uint theNumParam=0) ; ///< Set model parameters.
		void Set(const cDVector& theVectParam, const uint theNumParam=0) ; ///< Set model parameters.
		double Get(const uint theIndex, const uint theNumParam) ;
		cDVector& Get(const uint theNumParam);
		cMa& operator=(const cMa& theSrc) ; ///< Standard affectation
		void ReAllocProxyMeanParameters(uint theOldNParam=0) {};
		void UpdateProxyMeanParameters(void) {};
		double ComputeMean(uint theDate, const cRegArchValue& theData) const ; /// Compute the conditional mean value
		uint GetNParam(void) const ;
		uint GetNLags(void) const ;
		void ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, uint theBegIndex);
		void ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex);
		void ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex);
		void RegArchParamToVector(cDVector& theDestVect, uint theIndex) ;
		void VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex = 0) ;
		void GetParamName(uint theIndex, char** theName);
		void GetParamName(uint theIndex, string theName[]);
	} ;

}

#endif // _CMA_H_
