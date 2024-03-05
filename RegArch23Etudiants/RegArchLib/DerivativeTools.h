#pragma once
#ifndef _DERIVATIVETOOLS_H_
#define _DERIVATIVETOOLS_H_
/*!
\file cNumericDerivative.h
\brief Definition of the cNumericDerivative class
This class is used to compute the gradient and the Hessian
\author Jean-Baptiste DURAND, Ollivier TARAMASCO
\date Jan-31-2018 - Last change Jan-31-2018
*/
namespace RegArchLib {

	typedef double(*StdevFunction)(const cDVector& theTeta);
	typedef double(*LogDensFunction)(double theX, const cDVector& theTeta);
	typedef void(*GradStdevFunction)(const cDVector& theTeta, cDVector& theGrad);
	typedef void(*GradLogDensFunction)(double theX, const cDVector& theTeta, cDVector& theGrad);
	typedef void(*HessStdevFunction)(const cDVector& theTeta, cDMatrix& theHess);
	typedef void(*HessLogDensFunction)(double theX, const cDVector& theTeta, cDMatrix& theHess);


	double TrueLogDensity(double theX, const cDVector theTeta, LogDensFunction theLogDens, StdevFunction theStdev);
	double DiffTrueDensity(double theX, const cDVector& theTeta, LogDensFunction theDiffLogDens, StdevFunction theStdev);
	double Diff2TrueDensity(double theX, const cDVector& theTeta, LogDensFunction theDiff2LogDens, StdevFunction theStdev);
	void GradTrueLogDensity(double theX, const cDVector& theTeta, LogDensFunction theDiffLogDens, StdevFunction theStdev, GradStdevFunction theGradStdev, GradLogDensFunction theGradLogDensity, cDVector& theGrad);
	void GradTrueDiffLogDensity(double theX, const cDVector& theTeta, LogDensFunction theDiffLogDens, LogDensFunction theDiff2LogDens, StdevFunction theStdev, GradStdevFunction theGradStdev, GradLogDensFunction theGradDiffLogDensity, cDVector& theGrad);
	void HessTrueLogDensity(double theX, const cDVector& theTeta, LogDensFunction theDiffLogDens, LogDensFunction theDiff2LogDens, StdevFunction theStdev, GradStdevFunction theGradStdev, GradLogDensFunction theGradDiffLogDensity, HessStdevFunction theHessStdev, HessLogDensFunction theHessLogDens, cDMatrix& theHess);

}

#endif // _DERIVATIVETOOLS_H_
