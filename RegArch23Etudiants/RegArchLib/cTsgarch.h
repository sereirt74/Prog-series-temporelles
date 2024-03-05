#pragma once 
#ifndef _CTSARCH_H_
#define _CTSARCH_H_

#include "cAbstCondVar.h"
#include "cRegArchValue.h"
#include "cRegArchGradient.h"
/*!
	\file cTsgarch.h
	\brief Definition of the TSGARCH(p, q) class
	
	\date feb-18-2024 - Last change feb-18-2024

*/

/**
 * Modèle TS-GARCH(p, q) Taylor-Schwert Garch model
 *
 * h(t) = Cste + Somme_{i=1 à p} {Arch+(i)* |U(t-i)|} + Somme[j=1 à q}{Garch(j) h(t-j)}
 *
 */


namespace RegArchLib {

	/*! 
	 * \class cTsgarch
	 * \brief Class to implement a 'pure' TS-GARCH(p, q) model
	 */
	class cTsgarch: public cAbstCondVar
	{
	private :
		double mvConst ; ///< Constant part of GARCH(p, q) variance model.
		cDVector mvArch; ///< Vector of ARCH+ coefficients. 
		cDVector mvGarch ; ///< Vector of GARCH coefficients.
	public :
		cTsgarch(uint theNArch = 0, uint theNGarch=0) ; ///< A simple constructor
		cTsgarch(double theConst, cDVector& theArch, cDVector& theGarch) ; ///< Another constructor
		cTsgarch(const cTsgarch& theTsgarch);
		virtual ~cTsgarch() ; ///< A simple destructor
		void Delete(void) ; /// Delete
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
		cTsgarch& operator=(const cTsgarch& theSrc); ///< Standard affectation
		void ReAllocProxyVarParameters(uint theOldNParam=0) {};
		void UpdateProxyVarParameters(void) {};
		double ComputeVar(uint theDate, const cRegArchValue& theData) const;	///< Return conditional variance.
		uint GetNParam(void) const ; ///< Number of parameters in that model part
		uint GetNLags(void) const ; ///< Number of past gradients required to compute gradient at current time t.
		void ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cAbstResiduals* theResiduals) ;
		void RegArchParamToVector(cDVector& theDestVect, uint theIndex);
		void VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex = 0) ;
		void ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals);
		void ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals);
		void GetParamName(uint theIndex, char** theName);
		void GetParamName(uint theIndex, string theName[]);
	} ;

}
#endif // _CTSARCH_H_
