#pragma once 
#ifndef _CSKEWTRESIDUALS_H_
#define _CSKEWTRESIDUALS_H_

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_cdf.h>
#include "cAbstResiduals.h"
#include "cRegArchValue.h"
#include "cRegArchGradient.h"

/*!
\file cSkewtResiduals.h
\brief Definition of the class for Skewt conditional distribution.

\author Jean-Baptiste DURAND, Ollivier TARAMASCO
\date dec-18-2006 - Last change feb-18-2011
*/


namespace RegArchLib {

	/*!
	* \class cSkewtResiduals
	* \brief  Class to implement a Skewt distribution with unitary variance for residuals
	*/
	class cSkewtResiduals : public cAbstResiduals
	{
	public:
		cSkewtResiduals(double theDof, double theGamma, bool theSimulFlag = true); ///< A constructor
		cSkewtResiduals(scDVector* theDistrParameter = NULL, bool theSimulFlag = true);	///< A constructor
		virtual ~cSkewtResiduals(); ///< A simple destructor
		virtual cAbstResiduals* PtrCopy() const; /// < Return a copy of *this		
	#ifdef _RDLL_
		void Print(void);
	#else
		void Print(ostream& theOut) const; ///< Print the distribution type
	#endif // _RDLL_
		void SetDefaultInitPoint(void);
		void Generate(const uint theNSample, cDVector& theYt) const; ///< Draw a sample from residual distribution 
		double LogDensity(double theX) const;
		/** Return the number of parameters in distribution */
		uint GetNParam(void) const;
		/** Compute the derivative of log density with respect to the random variable (theGradData[0]) \e and the gradient
		of log density with respect to the model parameters (other components in theGradData) */
		double DiffLogDensity(double theX) const;
		void ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData) const;
		void RegArchParamToVector(cDVector& theDestVect, uint theIndex) const;
		void VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex = 0);
		double ComputeEspAbsEps(void);
		void ComputeGradBetaEspAbsEps(cDVector &theGrad);
		void ComputeHessBetaEspAbsEps(cDMatrix &theHess);
		double Diff2LogDensity(double theX) const;
		void GradDiffLogDensity(double theX, const cDVector& theDistrParam, cDVector& theGrad);
		void ComputeHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals);
		void GetParamName(uint theIndex, char** theName);
		void GetParamName(uint theIndex, string theName[]);
	};
}
#endif //_CSKEWTRESIDUALS_H_
