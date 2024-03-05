#pragma once 
#ifndef _CNGARCH_H_
#define _CNGARCH_H_

#include "cAbstCondVar.h"
#include "cRegArchValue.h"
#include "cRegArchGradient.h"
/*!
	\file cNgarch.h
	\brief Definition of the NGARCH(p, q) class
	
	\author Jean-Baptiste DURAND, Ollivier TARAMASCO 
	\date dec-18-2006 - Last change feb-18-2011
*/
// The variance model is given by:
// h(t) = sigma(t)^2 = mvConst + sum_{i=1}{p}{mvArch[i]*(U(t-i) - mvTheta*sigma(t-i))^2} + sum_{j=1}^{q}{mvGarch[j]*h(t-j)}
//

namespace RegArchLib {

	/*! 
	 * \class cNgarch
	 * \brief Class to implement a 'pure' NGARCH(p, q) model
	 */
	class cNgarch: public cAbstCondVar
	{
	private :
		double mvConst ; ///< Constant part of NGARCH(p, q) variance model.
		double mvTheta ; ///< Theta coefficient of NGARCH(p,q) variance model.
        cDVector mvArch ; ///< Vector of ARCH coefficients. 
		cDVector mvGarch ; ///< Vector of GARCH coefficients.
	public :
		cNgarch(uint theNArch = 0, uint theNNgarch=0) ; ///< A simple constructor
		cNgarch(double theConst, double theTheta, cDVector& theArch, cDVector& theGarch) ; ///< Another constructor
		cNgarch(const cNgarch& theNgarch);
		virtual ~cNgarch() ; ///< A simple destructor
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
		cNgarch& operator=(const cNgarch& theSrc); ///< Standard affectation
		void ReAllocProxyVarParameters(uint theOldNParam=0) {};
		void UpdateProxyVarParameters(void) {};
		double ComputeVar(uint theDate, const cRegArchValue& theData) const;	///< Return conditional variance.
		uint GetNParam(void) const ; ///< Number of parameters in that model part
		uint GetNLags(void) const ; ///< Number of past gradients required to compute gradient at current time t.
		void ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cAbstResiduals* theResiduals) ;
		void RegArchParamToVector(cDVector& theDestVect, uint theIndex) ;
		void VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex = 0) ;
		void ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals);
		void ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals);
		void GetParamName(uint theIndex, char** theName);
		void GetParamName(uint theIndex, string theName[]);
	} ;

}
#endif // _CNGARCH_H_
