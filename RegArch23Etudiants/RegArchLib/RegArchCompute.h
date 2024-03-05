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

}

#endif //_REGARCHCOMPUTE_H_
