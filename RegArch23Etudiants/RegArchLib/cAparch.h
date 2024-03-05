#pragma once 
#ifndef _CAPARCH_H_
#define _CAPARCH_H_

#include "cAbstCondVar.h"
#include "cRegArchValue.h"
#include "cRegArchGradient.h"
#include "cRegArchHessien.h"
/*!
	\file cAparch.h
	\brief header for class cAparch

	\author Jean-Baptiste DURAND, Ollivier TARAMASCO
	\date dec-18-2006 - Last change feb-18-2011
*/

namespace RegArchLib {
	class cAparch : public cAbstCondVar
	{
	private :
		double	mvCste ;
		double	mvDelta ;
		cDVector mvArch ;
		cDVector mvGamma ;
		cDVector mvGarch ;
	public :
		cAparch(int theNArch=0, int theNGarch=0) ; ///< a simple constructor
		cAparch(const cAparch& theAparch);
		virtual ~cAparch(); ///< a simple destructor
		void Delete(void);
#ifndef _RDLL_
		void Print(ostream& theOut = cout) const; ///< print the parameters
#else
		void Print(void);
#endif //_RDLL_
		void SetDefaultInitPoint(double theMean, double theVar);
		void SetDefaultInitPoint(cRegArchValue& theValue);
		void ReAlloc(const uint theSize, const uint theNumParam = 0); ///< Allocation of the model parameters
		void ReAlloc(const cDVector& theVectParam, const uint theNumParam = 0); ///< Allocation of the model parameters
		void Set(const double theValue, const uint theIndex = 0, const uint theNumParam = 0); ///< Set model parameters.
		void Set(const cDVector& theVectParam, const uint theNumParam = 0); ///< Set model parameters.
		double Get(const uint theIndex, const uint theNumParam);
		cDVector& Get(const uint theNumParam);
		cAparch& operator=(const cAparch& theSrc);
		void ReAllocProxyVarParameters(uint theOldNParam = 0) {};
		void UpdateProxyVarParameters(void) {};
		double ComputeVar(uint theDate, const cRegArchValue& theValue) const;
		uint GetNParam(void) const; ///< Number of parameters in that model part
		uint GetNLags(void) const; ///< Number of past gradients required to compute gradient at current time t.
		void RegArchParamToVector(cDVector& theDestVect, uint theIndex = 0);
		void VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex = 0);
		void ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cAbstResiduals* theResiduals);
		void ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals);
		void ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals);
		void GetParamName(uint theIndex, char** theName);
		void GetParamName(uint theIndex, string theName[]);
	} ;
}
#endif // _CAPARCH_H_
