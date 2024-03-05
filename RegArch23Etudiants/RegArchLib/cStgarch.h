#pragma once 
#ifndef _CSTGARCH_H_
#define _CSTGARCH_H_

#include "cAbstCondVar.h"
#include "cRegArchValue.h"
#include "cRegArchGradient.h"
/*!
	\file cStgarch.h
	\brief Definition of the ST-GARCH(p, q) class
	
	\date feb-18-2024 - Last change feb-18-2024

*/

/**
 * Modèle ST-GARCH(p, q)  Smooth Transition GARCH
 *
 * h(t) = Const + somme_{i=1,p}{(ARCH[i] +F(U(t-i); Gamma))*U(t-i)^2} + somme_{j=1,q}{GARCH[j]*h[t-j]}
 * et F(x; Gamma) = 1/(1+exp(Gamma*x))
 *
 */


namespace RegArchLib {

	/*! 
	 * \class cStgarch
	 * \brief Class to implement a 'pure' ST-GARCH(p, q) model
	 */
	class cStgarch: public cAbstCondVar
	{
	private :
		double mvConst ; ///< Constant part of Stgarch(p, q) variance model.
		cDVector mvArch; ///< Vector of ARCH coefficients. 
		cDVector mvGarch ; ///< Vector of GARCH coefficients.
		double mvGamma; ///< Gamma value
	public:  
		cStgarch(uint theNArch = 0, uint theNGarch=0) ; ///< A simple constructor
		cStgarch(double theConst, cDVector& theArch, cDVector& theGarch, double theGamma) ; ///< Another constructor
		cStgarch(const cStgarch& theGarch);
		virtual ~cStgarch() ; ///< A simple destructor
		//cAbstCondVar* PtrCopy() const ; /// < Return a copy of *this				
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
		cStgarch& operator=(const cStgarch& theSrc); ///< Standard affectation
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
	private :
		double StF(double thex, double thegamma) { return 1.0 / (1 + exp(thex * thegamma)); }
	} ;

}
#endif // _CSTGARCH_H_