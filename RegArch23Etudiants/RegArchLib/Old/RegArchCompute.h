#pragma once 
#ifndef _REGARCHCOMPUTE_H_
#define _REGARCHCOMPUTE_H_

#include "StdAfxRegArchLib.h"

/*!
 \file RegArchCompute.h
 \brief Header for simulation / estimation of general RegArch models.

 \author Jean-Baptiste DURAND, Ollivier TARAMASCO
 \date dec-18-2006 - Last change feb-18-2011
*/
namespace RegArchLib {
	using namespace WrapperGslCppNameSpace;
	using namespace WrapperNloptCppNamespace;

extern void RegArchSimul(const uint theNSample, const cRegArchModel& theModel, cDVector& theYt, cDMatrix* theXt=NULL, cDMatrix* theXvt=NULL) ; ///< Simulation of a general RegArch Model
extern void RegArchSimul(const uint theNSample, const cRegArchModel& theModel, cRegArchValue& theData) ; ///< Simulation of a general RegArch Model
extern void FillValue(uint theDate, const cRegArchModel& theModel, cRegArchValue& theValue);
extern void FillValueForNumericGrad(uint theDate, const cRegArchModel& theModel, cRegArchValue& theValue, cNumericDerivative& theNumDeriv);
extern void FillValueForNumericGradAndHess(uint theDate, const cRegArchModel& theModel, cRegArchValue& theValue, cNumericDerivative& theNumDeriv);
extern double RegArchLLH(const cRegArchModel& theModel, cDVector* theYt, cDMatrix* theXt=NULL) ; ///< Log-Likelihood of a general RegArch Model
extern double RegArchLLH(const cRegArchModel& theModel, cRegArchValue& theData) ; ///< Log-Likelihood of a general RegArch Model
extern void RegArchGradLt(int theDate, cRegArchModel& theModel, cRegArchValue& theData, cRegArchGradient& theGradData, cDVector& theGradlt) ;
extern void RegArchLtAndGradLt(int theDate, cRegArchModel& theModel, cRegArchValue& theValue, cRegArchGradient& theGradData, double& theLt, cDVector& theGradlt);
extern void NumericRegArchGradLt(uint theDate, cRegArchModel& theModel, cRegArchValue* theValue, cDVector& theGradlt, double theh = 1e-6);
extern void RegArchGradLLH(cRegArchModel& theModel, cRegArchValue& theData, cDVector& theGradLLH);
extern void RegArchLLHAndGradLLH(cRegArchModel& theModel, cRegArchValue& theValue, double& theLLH, cDVector& theGradLLH) ;
extern void NumericRegArchHessLt(int theDate, cRegArchModel& theModel, cRegArchValue* theValue, cRegArchGradient* theGradData, cDMatrix& theHesslt, double theh=1e-6);
extern void RegArchHessLt(int theDate, cRegArchModel& theModel, cRegArchValue& theValue, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cDMatrix& theHesslt);
extern void RegArchGradAndHessLt(int theDate, cRegArchModel& theModel, cRegArchValue& theValue, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cDVector& theGradlt, cDMatrix& theHesslt); 
extern void RegArchLtGradAndHessLt(int theDate, cRegArchModel& theModel, cRegArchValue& theValue, double& thelt, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cDVector& theGradlt, cDMatrix& theHesslt);
extern void RegArchHessLLH(cRegArchModel& theModel, cRegArchValue& theValue, cDMatrix& theHessLLH);
extern void NumericRegArchGradLLH(cRegArchModel& theModel, cRegArchValue& theValue, cDVector& theGradLLH, double theh = 1e-3);
extern void NumericRegArchHessLLH(cRegArchModel& theModel, cRegArchValue& theValue, cDMatrix& theHessLLH, double theh = 1e-3) ;

#//ifndef _RDLL_
	extern void RegArchEstim(cRegArchModel& theModel, cRegArchValue& theValue, sGSLMultiMinResult& theResStruct, cRegArchModel& theResModel, cDVector* theInitPoint = NULL, eGSLMultiMinAlgoEnum theAlgo = eConjugateFR, double theStopValue = 1e-6, int theMaxIter = 100, bool theVerbose = false);
	extern void RegArchEstim(cRegArchModel& theModel, cRegArchValue& theValue, sGSLMultiMinResult& theResStruct, cRegArchModel& theResModel, cDVector* theInitPoint = NULL, sGSLMultiMinAlgoParam* theAlgoParam = NULL);
	extern void RegArchEstim(cRegArchModel& theModel, cRegArchValue& theValue, cNLOPTResult& theResStruct, cRegArchModel& theResModel, cDVector* theInitPoint = NULL, nlopt_algorithm theAlgo = NLOPT_GN_DIRECT, double theStopValue = 1e-6, double thefTol = 1e-6, double thexTol = 1e-6, double theNMaxSec = 60, int theMaxFuncEval = 9999999, bool theMinimisation = false, bool theVerbose = false);
	extern void RegArchEstim(cRegArchModel& theModel, cRegArchValue& theValue, cNLOPTResult& theResStruct, cRegArchModel& theResModel, cDVector* theInitPoint, cNLOPTAlgoParam* theAlgoParam);


//#endif // _RDLL_

extern void NumericComputeCov(cRegArchModel &theModel, cRegArchValue &theData, cDMatrix &theCov);
extern void RegArchComputeCov(cRegArchModel& theModel, cRegArchValue& theValue, cDMatrix& theCov);
extern void RegArchComputeCov(cRegArchModel& theModel, cRegArchValue& theValue, cDMatrix& theCov, int& theError);
extern void RegArchComputeI(cRegArchModel &theModel, cRegArchValue &theData, cDMatrix &theI);
extern void NumericRegArchHessLLHold(cRegArchModel& theModel, cRegArchValue& theValue, cDMatrix& theHessLLH, double theh = 1e-3);
extern void RegArchGradAndHessLt(int theDate, cRegArchModel& theModel, cRegArchValue& theValue, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cDVector& theGradlt, cDMatrix& theHesslt);
extern void RegArchComputeIAndJ(cRegArchModel &theModel, cRegArchValue &theData, cDMatrix &theI, cDMatrix &theJ);
extern void RegArchComputeCov(cRegArchModel& theModel, cRegArchValue& theValue, cDMatrix& theCov);
extern void RegArchStatTable(cRegArchModel &theModel, cRegArchValue& theValue, cDMatrix& theTable);
}

#endif //_REGARCHCOMPUTE_H_
