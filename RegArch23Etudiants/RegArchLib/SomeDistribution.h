#pragma once 
#ifndef _SOMEDISTRIBUTION_H_
#define _SOMEDISTRIBUTION_H_

#include "StdAfxError.h"
#include "StdAfxVectorAndMatrix.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include "cStudentResiduals.h"
/*!
\file StudentDistribution.h
\brief Somme functions for Student conditional distribution.

\author Jean-Baptiste DURAND, Ollivier TARAMASCO
\date mar-15-2018 - Last change mar-21-2019
*/
namespace RegArchLib {
	double StudentLogDensity(double theX, double theDof);
	void StudentGradLogDensity(double theX, double theDof, cDVector& theGrad);
	double SkewtLogDensity(double theX, double theDof, double theGamma);
	void SkewtGradLogDensity(double theX, double theDof, double theGamma, cDVector& theGrad);
	double SkewtDiffLogDensity(double theX, double theDof, double theGamma);
	double SkewtExpect(double theDof, double theGamma);
	double SkewtVar(double theDof, double theGamma);
	void GradSkewtExpect(double theDof, double theGamma, cDVector& theGrad);
	void GradSkewtVar(double theDof, double theGamma, cDVector& theGrad, cDVector* theGradExpVal=NULL);
	double MixNormSigma(double thep, double theSigma1, double theSigma2);
	double MixNormSigma(const cDVector& theDistrParam);
	void MixNormGradSigma(double thep, double theSigma1, double theSigma2, cDVector& theGrad);
	void MixNormGradSigma(const cDVector& theDistrParam, cDVector& theGrad);
	void MixNormHessSigma(double thep, double theSigma1, double theSigma2, cDMatrix& theHess);
	void MixNormHessSigma(const cDVector& theDistrParam, cDMatrix& theHess);
	double MixNormLogDensity(double theX, double thep, double theSigma1, double theSigma2);
	double MixNormLogDensity(double theX, const cDVector &theDistrParam);
	double MixNormDiffLogDensity(double theX, double thep, double theSigma1, double theSigma2);
	double MixNormDiffLogDensity(double theX, const cDVector& theDistrParameter);
	double MixNormDiff2LogDensity(double theX, double thep, double theSigma1, double theSigma2);
	double MixNormDiff2LogDensity(double theX, const cDVector& theDistrParameter);
	void MixNormGradDiffLogDens(double theX, double thep, double theSigma1, double theSigma2, cDVector& theGrad);
	void MixNormGradDiffLogDens(double theX, const cDVector& theDistrParam, cDVector& theGrad);
	void MixNormGradLogDensity(double theX, double thep, double theSigma1, double theSigma2, cDVector& theGrad);
	void MixNormGradLogDensity(double theX, const cDVector &theDistrParam, cDVector& theGrad);
	void MixNormHessLogDensity(double theX, double thep, double theSigma1, double theSigma2, cDMatrix& theHess);
	void MixNormHessLogDensity(double theX, const cDVector &theDistrParam, cDMatrix& theHess);
} // namespace

#endif //_SOMEDISTRIBUTION_H_
