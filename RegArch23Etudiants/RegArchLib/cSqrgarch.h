#pragma once 
#ifndef _CSQRGARCH_H_
#define _CSQRGARCH_H_

#include "cAbstCondVar.h"
#include "cRegArchValue.h"
#include "cRegArchGradient.h"
/*!
	\file cSqrgarch.h
	\brief Definition of the SQR-GARCH(p, q) class
	
	\date feb-18-2024 - Last change feb-18-2024

*/

/**
 * Modèle SQR-GARCH(p, q)
 *
 * h(t) = Const + somme_{i=1,p}{ARCH[i]*(U[t-i] - Gamma[i]*sqrt(h(t-i))^2} + somme_{j=1,q}{GARCH[j]*h[t-j]}
 *
 */


namespace RegArchLib {

	/*! 
	 * \class cSqrarch
	 * \brief Class to implement a 'pure' SQR-ARCH(p, q) model
	 */
	class cSqrgarch: public cAbstCondVar
	{
	private :
		double mvConst ; ///< Constant part of Sqrgarch(p, q) variance model.
		cDVector mvArch; ///< Vector of ARCH coefficients. 
		cDVector mvGarch ; ///< Vector of GARCH coefficients.
		cDVector mvGamma; ///< Vector of Gamma coefficients.
	public :
		cSqrgarch(uint theNArch = 0, uint theNGarch=0) ; ///< A simple constructor
		cSqrgarch(double theConst, cDVector& theArch, cDVector& theGarch, cDVector& theGamma) ; ///< Another constructor
		cSqrgarch(const cSqrgarch& theSqrgarch);
		virtual ~cSqrgarch() ; ///< A simple destructor
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
		cSqrgarch& operator=(const cSqrgarch& theSrc); ///< Standard affectation
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
#endif // _CSQRGARCH_H_