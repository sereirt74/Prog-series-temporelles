#include "StdAfxRegArchLib.h"

namespace RegArchLib {
	
	double StudentLogDensity(double theX, double theDof)
	{
		return log(gsl_ran_tdist_pdf(theX, theDof));
	}

	void StudentGradLogDensity(double theX, double theDof, cDVector& theGrad)
	{
	double myt2 = theX * theX;
	double myAux1 = myt2 + theDof;
		theGrad[0] = -(theDof + 1)*theX / myAux1;
	double myAux2 = log(myAux1) - gsl_sf_psi((theDof + 1) / 2.0) - log(theDof) + gsl_sf_psi(theDof / 2.0);
		theGrad[1] = -0.5*(myAux2 + (1 - myt2) / myAux1);
	}

	double SkewtLogDensity(double theX, double theDof, double theGamma)
	{
	double myAux = 2 * theGamma / (theGamma * theGamma + 1);
		if (theX < 0)
			return log(myAux) + StudentLogDensity(theX * theGamma, theDof);
		else
			return -log(myAux) + StudentLogDensity(theX / theGamma, theDof);

	}

	void SkewtGradLogDensity(double theX, double theDof, double theGamma, cDVector& theGrad)
	{
	cDVector myGradStudent(2);
	double myAux = 2 * theGamma / (theGamma*theGamma + 1);
		if (theX < 0)
		{
			StudentGradLogDensity(theX * theGamma, theDof, myGradStudent);
			theGrad[0] = myGradStudent[0] * theGamma;
			theGrad[1] = myGradStudent[1];
			theGrad[2] = 1.0 / theGamma - myAux + theX * myGradStudent[0];
		}
		else
		{
		StudentGradLogDensity(theX / theGamma, theDof, myGradStudent);
		theGrad[0] = myGradStudent[0] / theGamma;
		theGrad[1] = myGradStudent[1];
		theGrad[2] = 1.0 / theGamma - myAux - theX / (theGamma * theGamma) * myGradStudent[0];
		}
	}

	double SkewtDiffLogDensity(double theX, double theDof, double theGamma)
	{
	cDVector myGradStudent(2);

		if (theX < 0)
		{
			StudentGradLogDensity(theX * theGamma, theDof, myGradStudent);
			return(myGradStudent[0] * theGamma);
		}
		else
		{
			StudentGradLogDensity(theX / theGamma, theDof, myGradStudent);
			return(myGradStudent[0] / theGamma);
		}

	}

	double SkewtExpect(double theDof, double theGamma)
	{
		return 2 * (theGamma * theGamma - 1)*sqrt(theDof)*gsl_sf_gamma((theDof + 1) / 2) / (SQRT_PI*theGamma*(theDof - 1)*gsl_sf_gamma(theDof / 2));
	}

	double SkewtVar(double theDof, double theGamma)
	{
	double myMean = SkewtExpect(theDof, theGamma);
		return theDof * (pow(theGamma, 4) - theGamma * theGamma + 1) / ((theDof - 2)*theGamma*theGamma) - myMean * myMean;
	}

	void GradSkewtExpect(double theDof, double theGamma, cDVector& theGrad)
	{
	double myGamma0 = gsl_sf_gamma(theDof / 2);
	double myGamma1 = gsl_sf_gamma((theDof + 1) / 2);
	double myPsi00 = gsl_sf_psi(theDof / 2);
	double myPsi01 = gsl_sf_psi((theDof + 1) / 2);
	double myGammaC = theGamma * theGamma;
	double mykm1 = theDof - 1;
	double mySqrtk = sqrt(theDof);
	double myNum1 = (theDof * mykm1 * (myPsi01 - myPsi00) - theDof - 1) * myGamma1 * (myGammaC - 1);
	double myDen1 = SQRT_PI * myGamma0 * mykm1 * mykm1 * mySqrtk * theGamma;
		theGrad[0] = myNum1 / myDen1;  // d/dDOF
	double myNum2 = 2 * mySqrtk * myGamma1 * (myGammaC + 1);
	double myDen2 = SQRT_PI * myGamma0 * mykm1 * myGammaC;
		theGrad[1] = myNum2 / myDen2; // d/dGamma
	}

	static void GradSkewtExpectedSquare(double theDof, double theGamma, cDVector& theGrad)
	{
		double myGamma2 = theGamma * theGamma;
		double myGamma3 = theGamma * myGamma2;
		double myGamma4 = theGamma * myGamma3;
		double mykm2 = theDof - 2;
		theGrad[0] = -(2 * (myGamma4 - myGamma2 + 1)) / (mykm2 * mykm2 * myGamma2);
		theGrad[1] = 2 * theDof * (myGamma4 - 1) / (myGamma3 * mykm2);
	}

	void GradSkewtVar(double theDof, double theGamma, cDVector& theGrad, cDVector* theGradExpVal)
	{
	cDVector myGradE2(2);
	double myExpVal = SkewtExpect(theDof, theGamma);
		GradSkewtExpectedSquare(theDof, theGamma, myGradE2);
	cDVector myGradExpVal(2);
		if (theGradExpVal == NULL)
			GradSkewtExpect(theDof, theGamma, myGradExpVal);
		else
			myGradExpVal = *theGradExpVal;
		theGrad = myGradE2 - 2 * myExpVal * myGradExpVal;
	}

	double MixNormSigma(double thep, double theSigma1, double theSigma2)
	{
		return sqrt(thep * theSigma1 * theSigma1 + (1 - thep) * theSigma2 * theSigma2);
	}

	double MixNormSigma(const cDVector& theDistrParam)
	{
		return sqrt(theDistrParam[0] * theDistrParam[1] * theDistrParam[1] + (1 - theDistrParam[0])*theDistrParam[2] * theDistrParam[2]);
	}

	void MixNormGradSigma(double thep, double theSigma1, double theSigma2, cDVector& theGrad)
	{
	double myStd = MixNormSigma(thep, theSigma1, theSigma2);
		theGrad[0] = (theSigma1 * theSigma1 - theSigma2 * theSigma2) / (2 * myStd);
		theGrad[1] = thep * theSigma1 / myStd;
		theGrad[2] = (1 - thep) * theSigma2 / myStd;
	}

	void MixNormGradSigma(const cDVector& theDistrParam, cDVector& theGrad)
	{
		MixNormGradSigma(theDistrParam[0], theDistrParam[1], theDistrParam[2], theGrad);
	}

	void MixNormHessSigma(double thep, double theSigma1, double theSigma2, cDMatrix& theHess)
	{
	double myVar1 = theSigma1 * theSigma1;
	double myVar2 = theSigma2 * theSigma2;
	double myStd = MixNormSigma(thep, theSigma1, theSigma2);
	double myStd2 = myStd * myStd;
	double myStd3 = myStd * myStd2;
	double myq = 1 - thep;
		theHess[0][0] = -pow(myVar1 - myVar2, 2) / (4 * myStd3);
		theHess[0][1] = theHess[1][0] = theSigma1 * (myVar2 + myStd2) / (2 * myStd3);
		theHess[0][2] = theHess[2][0] = -theSigma2 * (myVar1 + myStd2) / (2 * myStd3);
		theHess[1][1] = thep * myq * myVar2 / myStd3;
		theHess[1][2] = theHess[2][1] = -thep * myq * theSigma1 * theSigma2 / myStd3;
		theHess[2][2] = thep * myq * myVar1 / myStd3;
	}

	void MixNormHessSigma(const cDVector& theDistrParam, cDMatrix& theHess)
	{
		MixNormHessSigma(theDistrParam[0], theDistrParam[1], theDistrParam[2], theHess);
	}
		
	double MixNormLogDensity(double theX, double thep, double theSigma1, double theSigma2)
	{
		return log(thep*gsl_ran_gaussian_pdf(theX, theSigma1) + (1 - thep)*gsl_ran_gaussian_pdf(theX, theSigma2));

	}

	double MixNormLogDensity(double theX, const cDVector &theDistrParam)
	{
		return log(theDistrParam[0]*gsl_ran_gaussian_pdf(theX, theDistrParam[1]) + (1 - theDistrParam[0])*gsl_ran_gaussian_pdf(theX, theDistrParam[2]));

	}

	double MixNormDiffLogDensity(double theX, double thep, double theSigma1, double theSigma2)
	{
	double myf1 = gsl_ran_gaussian_pdf(theX, theSigma1);
	double myf2 = gsl_ran_gaussian_pdf(theX, theSigma2);
	double myq = 1 - thep;
	double myDen = thep*myf1 + myq*myf2;
		return -theX * (thep*myf1 / (theSigma1*theSigma1) + myq*myf2 / (theSigma2*theSigma2)) / myDen;
	}

	double MixNormDiffLogDensity(double theX, const cDVector& theDistrParameter)
	{
		return MixNormDiffLogDensity(theX, theDistrParameter[0], theDistrParameter[1], theDistrParameter[2]);
	}

	double MixNormDiff2LogDensity(double theX, double thep, double theSigma1, double theSigma2)
	{
		double myf1 = gsl_ran_gaussian_pdf(theX, theSigma1);
		double myf2 = gsl_ran_gaussian_pdf(theX, theSigma2);
		double myNum1 = (thep *myf1 / pow(theSigma1, 2) + (1 - thep)*myf2 / pow(theSigma2, 2));
		double myDen1 = thep * myf1 + (1 - thep)*myf2;
		double myNum1P = -theX * (thep*myf1 / pow(theSigma1, 4) + (1 - thep)*myf2 / pow(theSigma2, 4));
		double myDen1P = -theX * myNum1;
		double myDiff1LogDensSurx = myNum1 / myDen1;
		return -myDiff1LogDensSurx - theX * (myNum1P*myDen1 - myNum1 * myDen1P) / pow(myDen1, 2);
	}

	double MixNormDiff2LogDensity(double theX, const cDVector& theDistrParameter)
	{
		return MixNormDiff2LogDensity(theX, theDistrParameter[0], theDistrParameter[1], theDistrParameter[2]);
	}

	void MixNormGradDiffLogDens(double theX, double thep, double theSigma1, double theSigma2, cDVector& theGrad)
	{
	double myf1 = gsl_ran_gaussian_pdf(theX, theSigma1);
	double myf2 = gsl_ran_gaussian_pdf(theX, theSigma2);
	double myq = 1 - thep;
	double myVar1 = theSigma1 * theSigma1;
	double myVar2 = theSigma2 * theSigma2;
/*
	double myNum = theX * (myVar2 - myVar1);
	double myDen = theSigma1 * theSigma2*(2 * thep*myq*theSigma1*theSigma2 + thep * thep*myVar2*myf1 + myq*myq*myVar1*myf2);
		theGrad[0] = myNum / myDen;
		myDen = (thep*myf1 + myq*myf2);
		myNum = thep * myf1 / theSigma1 * (theX*theX / myVar1 - 1);
		myNum *= myq*myf2*(1 / myVar1 - 1 / myVar2);
		myNum -= 2 * myDen / myVar1;
		myDen *= myDen;
		theGrad[1] = myNum / myDen;
		myDen = (thep*myf1 + myq*myf2);
		myNum = myq*myf2 / theSigma2 * (theX*theX / myVar2 - 1);
		myNum *= thep*myf1*(1 / myVar2 - 1 / myVar1);
		myNum -= 2 * myDen / myVar2;
		myDen *= myDen;
		theGrad[2] = myNum / myDen;
*/
	double myNum = -theX * (thep * myf1 / myVar1 + myq * myf2 / myVar2);
	double myDen = thep * myf1 + myq * myf2;
	double myDen2 = myDen * myDen;
	double myNumP = -theX * (myf1 / myVar1 - myf2 / myVar2);
	double myDenP = myf1 - myf2;
		theGrad[0] = (myNumP * myDen - myNum * myDenP) / myDen2;
	double myX2 = theX * theX;
		myNumP = -theX * thep * (myX2 - 3 * myVar1) * myf1 / pow(theSigma1, 5);
		myDenP = thep * myf1 * (myX2 - myVar1) / pow(theSigma1, 3);
		theGrad[1] = (myNumP * myDen - myNum * myDenP) / myDen2;
		myNumP = -theX * myq * (myX2 - 3 * myVar2) * myf2 / pow(theSigma2, 5);
		myDenP = myq * myf2 * (myX2 - myVar2) / pow(theSigma2, 3);
		theGrad[2] = (myNumP * myDen - myNum * myDenP) / myDen2;
	}

	void MixNormGradDiffLogDens(double theX, const cDVector& theDistrParam, cDVector& theGrad)
	{
		MixNormGradDiffLogDens(theX, theDistrParam[0], theDistrParam[1], theDistrParam[2], theGrad);
		
	}

	void MixNormGradLogDensity(double theX, double thep, double theSigma1, double theSigma2, cDVector& theGrad)
	{
	double myf1 = gsl_ran_gaussian_pdf(theX, theSigma1);
	double myf2 = gsl_ran_gaussian_pdf(theX, theSigma2);
	double myDen = thep * myf1 + (1 - thep)*myf2;
		theGrad[0] = (myf1 - myf2) / myDen;
		theGrad[1] = -thep * (theSigma1*theSigma1 - theX * theX)*myf1 / (pow(theSigma1, 3)*myDen);
		theGrad[2] = -(1-thep) * (theSigma2*theSigma2 - theX * theX)*myf2 / (pow(theSigma2, 3)*myDen);
	}

	void MixNormGradLogDensity(double theX, const cDVector& theDistrParam, cDVector& theGrad)
	{
		MixNormGradLogDensity(theX, theDistrParam[0], theDistrParam[1], theDistrParam[2], theGrad);

	}
		
	void MixNormHessLogDensity(double theX, double thep, double theSigma1, double theSigma2, cDMatrix& theHess)
	{
	double myf1 = gsl_ran_gaussian_pdf(theX, theSigma1);
	double myf2 = gsl_ran_gaussian_pdf(theX, theSigma2);
	double myq = 1 - thep;
	double myDen1 = thep * myf1 + myq * myf2;
	double myDen1c = myDen1 * myDen1;
	double myVar1 = theSigma1 * theSigma1;
	double myVar2 = theSigma2 * theSigma2;
	double myX2 = theX * theX;
		theHess[0][0] = -(myf1 - myf2)*(myf1 - myf2) / myDen1c;
		theHess[1][0] = theHess[0][1] = -(myVar1 - myX2) * myf1 * myf2 / myDen1c / pow(theSigma1, 3);
		theHess[2][0] = theHess[0][2] = (myVar2 - myX2) * myf1 * myf2 / myDen1c / pow(theSigma2, 3);
	double myNum2 = -thep * (myVar1 - myX2)*myf1;
	double myDen2 = pow(theSigma1, 3)*myDen1;
	double myNum2P = -thep * (myVar1 * myVar1 + 2 * myVar1 * myX2 - myX2 * myX2) * myf1 / pow(theSigma1, 3);
	double myDen2P = thep * (2 * myVar1 + myX2) * myf1 + 3 * myq * myVar1 * myf2;
		theHess[1][1] = (myNum2P * myDen2 - myNum2 * myDen2P) / (myDen2 * myDen2);
		myNum2P = 0 ;
		myDen2P = pow(theSigma1 / theSigma2, 3) * myq * (myX2 - myVar2) * myf2;
		theHess[1][2] = theHess[2][1] = (myNum2P * myDen2 - myNum2 * myDen2P) / (myDen2 * myDen2);
	double myNum3 = -myq * (myVar2 - myX2) * myf2;
	double myDen3 = pow(theSigma2, 3) * myDen1;
	double myNum3P = -myq * (myVar2 * myVar2 + 2 * myVar2 * myX2 - myX2 * myX2) * myf2 / pow(theSigma2, 3);
	double myDen3P = myq * (2 * myVar2 + myX2) * myf2 + 3 * thep * myVar2 * myf1;
		theHess[2][2] = (myNum3P * myDen3 - myNum3 * myDen3P) / (myDen3 * myDen3);
	}

	void MixNormHessLogDensity(double theX, const cDVector& theDistrParam, cDMatrix& theHess)
	{
		MixNormHessLogDensity(theX, theDistrParam[0], theDistrParam[1], theDistrParam[2], theHess);
	}

} //namespace

