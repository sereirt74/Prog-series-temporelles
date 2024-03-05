#pragma once
#ifndef _CNUMERICDERIVATIVE_H_
#define _CNUMERICDERIVATIVE_H_
/*!
\file cNumericDerivative.h
\brief Definition of the cNumericDerivative class
This class is used to compute the gradient and the Hessian 
\author Jean-Baptiste DURAND, Ollivier TARAMASCO
\date Jan-31-2018 - Last change Jan-31-2018
*/
namespace RegArchLib {
	/*!
	* \class cNumericDerivative
	* \brief Data structure to store histories of gradients for empirical residuals,
	* conditional means and variances as stacks
	*/

	class cNumericDerivative
	{
	private:
		uint		mvNParam; ///< Numbert of parameters
		uint		mvNDistrParam; ///< Numbert of parameters of the conditional distribution
		double		mvh;
	public:
		double*			mh1;
		double**		mh2;
		cRegArchValue*	mValueForGrad;
		cRegArchValue**	mValueForHess;
		cDVector		mLogDensForGrad;
		cDMatrix		mLogDensForHess;
		cDVector		mGradDiffForGrad;

	public:
		cNumericDerivative(uint theNParam = 0, uint theNDistrParam=0, double theh = 1e-2, uint theNObs = 0);
		cNumericDerivative(uint theNParam, uint theNDistrParam, double theh, cDVector& theYt, cDMatrix* theXt, cDMatrix* theXvt);
		virtual ~cNumericDerivative();
		uint GetNParam(void);
		uint GetNDistrParam(void);
		double Geth(void);
	};

} // namespace

#endif // _CNUMERICDERIVATIVE_H_
