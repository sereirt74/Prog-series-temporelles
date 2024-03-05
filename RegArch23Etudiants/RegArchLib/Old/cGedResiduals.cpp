#include "StdAfxRegArchLib.h"
/*!
\file cGedResiduals.cpp
\brief implementation of the class for Ged conditional distribution.

\author Jean-Baptiste DURAND, Ollivier TARAMASCO
\date dec-18-2006 - Last change apr-7-2019
*/
namespace RegArchLib {
	/*!
	* \fn cGedResiduals::cGedResiduals(double theDof, bool theSimulFlag)
	* \param double theBeta: beta parameter
	* \param bool theSimulFlag: true if created for simulation
	* \details: mvBool is initialised by ce cAbstResiduals constructor
	*/
	cGedResiduals::cGedResiduals(double theBeta, bool theSimulFlag) : cAbstResiduals(eGed, NULL, theSimulFlag)
	{
		mDistrParameter = new cDVector(1, theBeta);
		(*mDistrParameter)[0] = theBeta;
		MESS_CREAT("cGedResiduals")
	}

	/*!
	* \fn cGedResiduals::cGedResiduals(const cDVector* theDistrParameter, bool theSimulFlag): cAbstResiduals(eGed, theDistrParameter, theSimulFlag)
	* \param const cDVector* theDistrParameter: theDistrParameter[0] = d.o.f.
	* \param bool theSimulFlag: true if created for simulation
	* \details: mvBool is initialised by ce cAbstResiduals constructor
	*/
	cGedResiduals::cGedResiduals(cDVector* theDistrParameter, bool theSimulFlag) : cAbstResiduals(eGed, theDistrParameter, theSimulFlag)
	{
		MESS_CREAT("cGedResiduals")
	}

	cGedResiduals::cGedResiduals(const cGedResiduals& theSrc) : cAbstResiduals(eGed, theSrc.mDistrParameter, (theSrc.mtR != NULL))
	{
		MESS_CREAT("cGedResiduals")
	}

	/*!
	* \fn cGedResiduals::~cGedResiduals()
	*/
	cGedResiduals::~cGedResiduals()
	{
		MESS_DESTR("cGedResiduals")
	}

	/*!
	* \fn cAbstCondVar* cNormResiduals::::PtrCopy()
	*/


	/*
	cAbstResiduals* cGedResiduals::PtrCopy() const
	{
		return cAbstResidualsPtrCopy<cGedResiduals>();
	}
	*/

	/*!
	* \fn void cGedResiduals::Generate(const uint theNSample, cDVector& theYt) const
	* \param const uint theNSample: the sample size
	* \param cDVector& theYt: the output vector
	*/
	void cGedResiduals::Generate(const uint theNSample, cDVector& theYt) const
	{
		theYt.ReAlloc(theNSample);

	double myAlpha = sqrt(gsl_sf_gamma(1 /(*mDistrParameter)[0]) / gsl_sf_gamma(3 /(*mDistrParameter)[0]));
		for (uint t = 0; t < theNSample; t++)
			theYt[t] = gsl_ran_exppow(mtR, myAlpha, (*mDistrParameter)[0]) ;
	}

	/*!
	* \fn void cGedResiduals::Print(ostream& theOut) const
	* \param ostream& theOut: the output stream, default cout.
	*/
#ifndef _RDLL_
	void cGedResiduals::Print(ostream& theOut) const
	{
		theOut << "Conditional Ged Distribution with: " << endl;
		theOut << " Beta=" << (*mDistrParameter)[0] << endl;
	}
#else
	void cGedResiduals::Print(void)
	{
		Rprintf("Conditional Ged Distribution with:\n");
		Rprintf(" Beta=%d\n", (*mDistrParameter)[0]);

	}
#endif // _RDLL_

	void cGedResiduals::SetDefaultInitPoint(void)
	{
		(*mDistrParameter)[0] = 2.5;
	}

	static double GedAlpha(double theBeta)
	{
		return sqrt(gsl_sf_gamma(1.0 / theBeta) / gsl_sf_gamma(3.0 / theBeta)); 
	}
	
	static double GedDensity(double theX, double theBeta)
	{
	double myAlpha = GedAlpha(theBeta);
		return theBeta*exp(-(pow(fabs(theX) / myAlpha, theBeta))) / (2 * myAlpha*gsl_sf_gamma(1.0 / theBeta));

	}

	double cGedResiduals::LogDensity(double theX) const
	{
		return log(GedDensity(theX, (*mDistrParameter)[0]));
	}

	/*!
	* \fn double cGedResiduals::GetNParam(void) const
	* \param void.
	* \brief return 2 : two parameters for Ged residuals.
	*/
	uint cGedResiduals::GetNParam(void) const
	{
		return 1;
	}

	double cGedResiduals::DiffLogDensity(double theX) const
	{
	double myBeta = (*mDistrParameter)[0];
	double myAlpha = GedAlpha(myBeta);

		return -myBeta*pow(fabs(theX), myBeta) / (theX*pow(myAlpha, myBeta));
	}

	/*!
	* \fn static void GradLogDensity(double theX, cDVector& theGrad, cDVector& theDistrParam)
	* \brief Compute the derivative of log density with respect to the random variable (theGrad[0]) \e and the gradient
	* of log density with respect to the model parameters (other components in theGrad)
	* \param theX double: value of the random variable
	* \param theGrad cDVector&: concatenation of derivatives with respect to the random variable and the model parameters
	* \param theDistrParam cDVector&: value of the distribution parameters
	*/
	static void GradLogDensity(double theX, cDVector& theGrad, const cDVector& theDistrParam, uint theBegIndex)
	{
	double myBeta = theDistrParam[0];
	double myGamma1 = gsl_sf_gamma(1.0 / myBeta);
	double myGamma3 = gsl_sf_gamma(3.0 / myBeta);
	double myPsi01 = gsl_sf_psi(1.0 / myBeta);
	double myPsi03 = gsl_sf_psi(3.0 / myBeta);
	double myAlpha = sqrt(myGamma1 / myGamma3);
	double myAlphaPBeta = pow(myAlpha, myBeta);
	double myXPBeta = pow(fabs(theX), myBeta);
	double myAux = myXPBeta * 2 * (log(myAlpha) - log(abs(theX))) * myBeta * myBeta;
		myAux += (2 * myAlphaPBeta + myXPBeta * (3 * myPsi03 - myPsi01)) * myBeta;
		myAux += (3 * myPsi01 - 3 * myPsi03) * myAlphaPBeta;
		myAux /= 2 * myAlphaPBeta * myBeta * myBeta;
		theGrad[theBegIndex] = myAux;
	}

	/*!
	* \fn void cGedResiduals::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData)
	* \brief Compute the derivative of log density with respect to the random variable (theGradData[0]) \e and the gradient
	* of log density with respect to the model parameters (other components in theGradData)
	* \param theDate uint: time at which gradient is computed
	* \param theValue const cRegArchValue&: value of the random variable
	* \param theGradData cRegArchGradient&: concatenation of derivatives with respect to the random variable and the model parameters
	*/
	void cGedResiduals::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData) const
	{
	uint myBegIndex = theGradData.GetNMeanParam() + theGradData.GetNVarParam();
		GradLogDensity(theValue.mEpst[theDate], theGradData.mCurrentGradLogDens, *mDistrParameter, myBegIndex);
		theGradData.mCurrentDiffLogDensity = DiffLogDensity(theValue.mEpst[theDate]);
	}

	void cGedResiduals::RegArchParamToVector(cDVector& theDestVect, uint theIndex) const
	{
		if (theDestVect.GetSize() < theIndex - 1)
			throw cError("Wrong size");
		theDestVect[theIndex] =(*mDistrParameter)[0];
	}

	void cGedResiduals::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
		if (1 + theIndex > theSrcVect.GetSize())
			throw cError("Wrong size");
		mDistrParameter[0] = theSrcVect[theIndex];
	}

	double cGedResiduals::ComputeEspAbsEps(void)
	{
	double myBeta =(*mDistrParameter)[0];
		return gsl_sf_gamma(2/myBeta)/sqrt(gsl_sf_gamma(1/myBeta)*gsl_sf_gamma(3/myBeta));
	}

	void cGedResiduals::ComputeGradBetaEspAbsEps(cDVector& theGrad)
	{
	double myBeta =(*mDistrParameter)[0];
	double myGamma1 = gsl_sf_gamma(1.0 / myBeta);
	double myGamma2 = gsl_sf_gamma(2.0 / myBeta);
	double myGamma3 = gsl_sf_gamma(3.0 / myBeta);
	double myPsi01 = gsl_sf_psi(1.0 / myBeta);
	double myPsi02 = gsl_sf_psi(2.0 / myBeta);
	double myPsi03 = gsl_sf_psi(3.0 / myBeta);
		theGrad[0] = myGamma2 * (3 * myPsi03 - 4 * myPsi02 + myPsi01) / (2 * sqrt(myGamma1*myGamma3)*myBeta*myBeta);
	}

	void cGedResiduals::ComputeHessBetaEspAbsEps(cDMatrix &theHess)
	{
	double myBeta =(*mDistrParameter)[0];
	double myGamma1 = gsl_sf_gamma(1.0 / myBeta);
	double myGamma2 = gsl_sf_gamma(2.0 / myBeta);
	double myGamma3 = gsl_sf_gamma(3.0 / myBeta);
	double myPsi01 = gsl_sf_psi(1.0 / myBeta);
	double myPsi02 = gsl_sf_psi(2.0 / myBeta);
	double myPsi03 = gsl_sf_psi(3.0 / myBeta);
	double myPsi11 = gsl_sf_psi_1(1.0 / myBeta);
	double myPsi12 = gsl_sf_psi_1(2.0 / myBeta);
	double myPsi13 = gsl_sf_psi_1(3.0 / myBeta);
	double myAux = ((12 * myPsi03 - 16 * myPsi01 + 4 * myPsi01)*myBeta + 18 * myPsi13 - 9 * myPsi03 * myPsi03 + (24 * myPsi01 - 6 * myPsi01) * myPsi03 - 16 * myPsi12 - 16 * myPsi01 * myPsi01 + 8 * myPsi01 * myPsi02 + 2 * myPsi11 - myPsi01 * myPsi01);
		myAux *= -myGamma2 / (4 * sqrt(myGamma1 * myGamma3) * pow(myBeta, 4));
		theHess[0][0] = myAux;
	}

	double cGedResiduals::Diff2LogDensity(double theX) const
	{
	double myBeta =(*mDistrParameter)[0];
		return (myBeta - 1) / theX * DiffLogDensity(theX);
	}

	void cGedResiduals::GradDiffLogDensity(double theX, const cDVector& theDistrParam, cDVector& theGrad)
	{
	double myBeta =(*mDistrParameter)[0];
	double myGamma1 = gsl_sf_gamma(1.0 / myBeta);
	double myGamma3 = gsl_sf_gamma(3.0 / myBeta);
	double myPsi01 = gsl_sf_psi(1.0 / myBeta);
	double myPsi03 = gsl_sf_psi(3.0 / myBeta);
	double myAlpha = sqrt(myGamma1 / myGamma3);
	double myXPBeta = pow(fabs(theX), myBeta);
	double myAux = myXPBeta * (2 * (log(myAlpha) - log(fabs(theX))) * myBeta + 3 * myPsi03 - myPsi01 - 2);
		myAux /= 2 * theX * pow(myAlpha, myBeta);
		theGrad[0] = myAux;
	}

	static void HessLogDensity(double theX, cDMatrix& theHess, const cDVector& theDistrParam, uint theBegIndex)
	{
	double myBeta = theDistrParam[0];
	double myGamma1 = gsl_sf_gamma(1.0 / myBeta);
	double myGamma3 = gsl_sf_gamma(3.0 / myBeta);
	double myPsi01 = gsl_sf_psi(1.0 / myBeta);
	double myPsi03 = gsl_sf_psi(3.0 / myBeta);
	double myPsi11 = gsl_sf_psi_1(1.0 / myBeta);
	double myPsi13 = gsl_sf_psi_1(3.0 / myBeta);
	double myAlpha = sqrt(myGamma1 / myGamma3);
	double myAlpha2 = myAlpha * myAlpha;
	double myXPBeta = pow(fabs(theX), myBeta);
	double myAlphaPBeta = pow(myAlpha, myBeta);
	double myBeta2 = myBeta*myBeta;
	double myBeta3 = myBeta*myBeta2;
	double myBeta4 = myBeta2*myBeta2;
	double myLogX = log(fabs(theX));
	double myLogAlpha2 = log(myAlpha2);
	double myAux = myXPBeta * pow(log(myAlpha2)  - 2 * myLogX, 2) * myBeta4;
		myAux += 4 * myXPBeta * (log(myAlpha) - myLogX) * (3 * myPsi03 -  myPsi01) * myBeta3;
		myAux += (4 * myAlphaPBeta + myXPBeta * pow(3 * myPsi03 - myPsi01, 2)) * myBeta2;
		myAux += (12 * (myPsi01 - myPsi03) * myAlphaPBeta + myXPBeta * (18 * myPsi13 - 2 * myPsi11)) * myBeta;
		myAux += 6 * (myPsi11 - 3 * myPsi13) * myAlphaPBeta;
		myAux /= -4 * myAlphaPBeta * myBeta4;
		theHess[theBegIndex][theBegIndex] = myAux;
	}

	void cGedResiduals::ComputeHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	double myX = theData.mEpst[theDate];
	cDVector myGradDiffLogDens = cDVector(1);
		GradDiffLogDensity(myX, *mDistrParameter, myGradDiffLogDens);
	uint myBegIndex = theHessData.GetNMeanParam() + theHessData.GetNVarParam();
		theHessData.mCurrentGradDiffLogDensity[myBegIndex] = myGradDiffLogDens[0];
		HessLogDensity(myX, theHessData.mCurrentHessDens, *mDistrParameter, myBegIndex);
	}

	void cGedResiduals::GetParamName(uint theIndex, char** theName)
	{
		uint myIndex = theIndex;
		sprintf(theName[myIndex++], "GED BETA");
	}

	void cGedResiduals::GetParamName(uint theIndex, string theName[])
	{
	uint myIndex = theIndex;
	char myChar[100];
		sprintf(myChar, "GED BETA");
		theName[myIndex++] = myChar;
	}

}//namespace
