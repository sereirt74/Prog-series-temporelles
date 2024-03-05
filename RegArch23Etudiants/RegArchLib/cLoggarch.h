#pragma once 
#ifndef _CLOGGARCH_H_
#define _CLOGGARCH_H_

#include "cAbstCondVar.h"
#include "cRegArchValue.h"
#include "cRegArchGradient.h"
/*!
	\file cLoggarch.h
	\brief Definition of the LOG-GARCH(p, q) class
	
	\date feb-18-2024 - Last change feb-18-2024

*/

/**
 * Modèle LOG-GARCH(p, q)  Smooth Transition GARCH
 *
 * Log(h(t)) = Const + somme_{i=1,p}{ARCH[i]*Log(U(t-i)^2)} + somme_{j=1,q}{GARCH[j]*log(h(t-j))}
 *
 */


namespace RegArchLib {

	/*! 
	 * \class cLogGarch
	 * \brief Class to implement a 'pure' LOG-ARCH(p, q) model
	 */
	class cLoggarch: public cAbstCondVar
	{
	private :
		double mvConst ; ///< Constant part of Loggarch(p, q) variance model.
		cDVector mvArch; ///< Vector of ARCH coefficients. 
		cDVector mvGarch ; ///< Vector of GARCH coefficients.
	public:  
		cLoggarch(uint theNArch = 0, uint theNGarch=0) ; ///< A simple constructor
		cLoggarch(double theConst, cDVector& theArch, cDVector& theGarch) ; ///< Another constructor
		cLoggarch(const cLoggarch& theLoggarch);
		virtual ~cLoggarch() ; ///< A simple destructor
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
		cLoggarch& operator=(const cLoggarch& theSrc); ///< Standard affectation
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
#endif // _CLOGGARCH_H_