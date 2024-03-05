#pragma once 
#ifndef _CSTDDEVINMEAN_H_
#define _CSTDDEVINMEAN_H_

#include "cAbstCondMean.h"
#include "cRegArchValue.h"
#include "cRegArchGradient.h"
/*!
	\file cStdDevInMean.h
	\brief header for class cStdDevInMean.
	\author Jean-Baptiste DURAND, Ollivier TARAMASCO
	\date dec-18-2006 - Last change feb-18-2011
*/

/*
	mm(t, Teta) = Teta x sqrt(h0) avec Teta = delta
	gradtetaF = sqrt(h0)
	hesstetaF = 0
*/
namespace RegArchLib {

	class cStdDevInMean : public cAbstCondMean
	{
	private :
		double mvStdDevInMean ;
	public :
		cStdDevInMean(double theStdDevInMean = 0.0) ;
		cStdDevInMean(const cStdDevInMean& theStdDevInMean) ; /// Recopy constructor
		~cStdDevInMean() ;
		//cAbstCondMean* PtrCopy(void) ; /// < Return a copy of *this
		void Delete(void) ;
	#ifndef _RDLL_
		void Print(ostream& theOut=cout) const ;
	#else
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
		cStdDevInMean& operator=(const cStdDevInMean& theSrc) ;
		void ReAllocProxyMeanParameters(uint theOldNParam=0) {};
		void UpdateProxyMeanParameters(void) {};
		double ComputeMean(uint theDate, const cRegArchValue& theData) const ;
		uint GetNParam(void) const ; ///< Return number of parameters
		uint GetNLags(void) const ;
		uint GetNu(void) const;
		uint GetNh(void) const;
		void RegArchParamToVector(cDVector& theDestVect, uint theIndex = 0) ;
		void VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex = 0) ;
		void ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, uint theBegIndex, cAbstResiduals* theResiduals) ;
		void GetNParamF(uint theNParam[3]) const;
		void ComputeGradF(uint theDate, const cRegArchValue& theData, cDVector& theGradF);
		void ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cDVector& theGradTheta, cDVector& theGradU, double& theGradH, uint& theNu, uint& theNh);
		void ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cFuncMeanAndVar& theDerivM);
		void ComputeHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData,cRegArchHessien& theHessData, uint theBegIndex, cAbstResiduals* theResiduals);
		void ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, const cDVector& theGradTheta, const cDVector& theGradU, double theGradH, uint theNu, uint theNh, cDMatrix& theHessTheta2, cDVector* theHessThetaU, cDVector& theHessThetaH, cDVector* theHessU2, cDVector& theHessUH, double& theHessH2);
		void ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, cFuncMeanAndVar& theDerivM);
		void ComputeHessF(uint theDate, const cRegArchValue& theData, cDMatrix& theHessF);
		void GetParamName(uint theIndex, char** theName);
		void GetParamName(uint theIndex, string theName[]);
	} ;
}

#endif // _CStdDevInMean_H_
