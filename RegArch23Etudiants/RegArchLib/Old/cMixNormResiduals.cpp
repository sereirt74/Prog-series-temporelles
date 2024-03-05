#include "StdAfxRegArchLib.h"

/*!
 \file cMixNormResiduals.cpp
 \brief implementation of the class for N(0, 1) conditional distribution.
 \author Jean-Baptiste DURAND, Ollivier TARAMASCO
 \date dec-18-2006 - Last change feb-18-2011
*/
namespace RegArchLib {
	/*!
	 * \fn cMixNormResiduals::cMixNormResiduals(cDVector* theDistrParameter, bool theSimulFlag):cAbstResiduals(eNormal, NULL, theSimulFlag)
	 * \param bool theSimulFlag: true if created for simulation
	 * \details: mvBool is initialised by cAbstResiduals constructor and theDistrParameter is never used.
	 */
	cMixNormResiduals::cMixNormResiduals(cDVector* theDistrParameter, bool theSimulFlag) : cAbstResiduals(eMixNorm, theDistrParameter, theSimulFlag)
	{
		MESS_CREAT("cMixNormResiduals")
	}

	cMixNormResiduals::cMixNormResiduals(double thep, double theSigma1, double theSigma2, bool theSimulFlag) : cAbstResiduals(eMixNorm, NULL, theSimulFlag)
	{
		mDistrParameter=new cDVector(3);
		(*mDistrParameter)[0] = thep;
		(*mDistrParameter)[1] = theSigma1;
		(*mDistrParameter)[2] = theSigma2;
	}

	cMixNormResiduals::cMixNormResiduals(const cMixNormResiduals& theSrc): cAbstResiduals(eMixNorm, theSrc.mDistrParameter, (theSrc.mtR != NULL))
	{
		MESS_CREAT("cMixNormResiduals")
	}


	/*!
	 * \fn cMixNormResiduals::~cMixNormResiduals
	 * \details: nothing to do here
	 */
	cMixNormResiduals::~cMixNormResiduals()
	{
		MESS_DESTR("cMixNormResiduals")
	}

	/*!
	 * \fn cAbstCondVar* cMixNormResiduals::::PtrCopy()
	 */
/*
	cAbstResiduals* cMixNormResiduals::PtrCopy() const
	{
	cMixNormResiduals *myMixNormResiduals = NULL;
	cDVector* myDistrParameter = new cDVector(mDistrParameter);
	bool mySimulFlag = (mtR != NULL);

		myMixNormResiduals = new cMixNormResiduals(myDistrParameter, mySimulFlag);
		delete myDistrParameter;
		return myMixNormResiduals;
	}
*/
	/*!
	 * \fn void cMixNormResiduals::Print(ostream& theOut) const
	 * \param ostream& theOut: the output stream, default cout.
	 */
#ifndef _RDLL_
	void cMixNormResiduals::Print(ostream& theOut) const
	{
		theOut << "Conditional mixture of gaussian distribution with:" << endl;
		theOut << " p=" << (*mDistrParameter)[0] << endl;
		theOut << " Sigma1=" << (*mDistrParameter)[1] << endl;
		theOut << " Sigma2=" << (*mDistrParameter)[2] << endl;
	}
#else
	void cMixNormResiduals::Print(void)
	{
		Rprintf("Conditional mixture of gaussian distribution with:\n");
		Rprintf(" p=%f\n", (*mDistrParameter)[0]);
		Rprintf(" Sigma1=%f\n", (*mDistrParameter)[1]);
		Rprintf(" Sigma2=%f\n", (*mDistrParameter)[2]);

	}
#endif // _RDLL_

	void cMixNormResiduals::SetDefaultInitPoint(void)
	{
		(*mDistrParameter)[0] = 0.5;
		(*mDistrParameter)[1] = 1;
		(*mDistrParameter)[2] = 3;
	}

	void cMixNormResiduals::Generate(uint theNSample, cDVector& theYt) const
	{
		theYt.ReAlloc(theNSample);

	double myStd = MixNormSigma(*mDistrParameter);
		for (uint t = 0; t < theNSample; t++)
		{
		double myp = gsl_ran_bernoulli(mtR, (*mDistrParameter)[0]);
			theYt[t] = (myp * gsl_ran_gaussian(mtR, (*mDistrParameter)[1]) + (1 - myp)*gsl_ran_gaussian(mtR, (*mDistrParameter)[2]))/myStd;

		}
	}

	/*!
	 * \fn double cMixNormResiduals::LogDensity(double theX) const
	 * \param double theX: the point where density is computed.
	 * \brief Compute the log density of N(0, 1)
	 */
	double cMixNormResiduals::LogDensity(double theX) const
	{
	// double myStd = MixNormSigma(mDistrParameter);
	//	return MixNormLogDensity(theX*myStd, (*mDistrParameter)[0], (*mDistrParameter)[1], (*mDistrParameter)[2]) + log(myStd);
		return TrueLogDensity(theX, *mDistrParameter, &MixNormLogDensity, &MixNormSigma);
	}

	/*!
	 * \fn double cMixNormResiduals::GetNParam(void)
	 * \param void.
	 * \brief return 0: no parameter for N(0,1) residuals.
	 */
	uint cMixNormResiduals::GetNParam(void) const
	{
		return 3;
	}

	/*!
	 * \fn double cMixNormResiduals::DiffLogDensity( theX)
	 * \param TheX. when called, theX is theValue.mEpst[theDate]
	 * lnfe(x; beta) = ln( sigma_nc ) + ln( fnc(sigma_nc*x ; beta) )
	 * fnc(x; beta) = p * fn( sigma_nc*x ; sigma1) + (1-p) * fn(sigma_nc*x ; sigma2)
	 * d/dx ( ln(fnc(x ; beta) ) = fnc'(x)/fnc(x)
	 * \return d/dx( ln(fe( theX ;beta) ) )
	 */
	double cMixNormResiduals::DiffLogDensity(double theX) const
	{	
	//double myStd = MixNormSigma(mDistrParameter);
			
	//	return myStd * MixNormDiffLogDensity(theX*myStd, (*mDistrParameter)[0], (*mDistrParameter)[1], (*mDistrParameter)[2]);
		return DiffTrueDensity(theX, *mDistrParameter, &MixNormDiffLogDensity, &MixNormSigma);
	}

	/*!
	 * \fn static void GradLogDensity(double theX, cDVector& theGrad)
	 * \brief Compute the derivative of log density of a Gaussian distribution with respect to the random variable (theGrad[0])
	 * \e and the gradient of log density with respect to the model parameters (other components in theGrad)
	 * \param theX double: value of the random variable
	 * \param theGrad cDVector&: concatenation of derivatives with respect to the random variable and the model parameters
	 * lnfe(x; beta) = ln( sigma_nc ) + ln( fnc(sigma_nc*x ; beta) )
	 * fnc(x; beta) = p * fn( sigma_nc*x ; sigma1) + (1-p) * fn(sigma_nc*x ; sigma2)
	 */
	static void GradLogDensity(double theX, cDVector& theGrad, const cDVector& theDistrParam, uint theBegIndex)
	{
	double myp = theDistrParam[0];
	double mySigma1 = theDistrParam[1];
	double mySigma2 = theDistrParam[2];
	double myStd = MixNormSigma(theDistrParam);
	double myDiffLogDensity= MixNormDiffLogDensity(theX*myStd, myp, mySigma1, mySigma2);
	cDVector myMixNormGradSigma(3);
		MixNormGradSigma(myp, mySigma1, mySigma2, myMixNormGradSigma);
	cDVector myGradMixNorm(3);
		MixNormGradLogDensity(theX*myStd, myp, mySigma1, mySigma2, myGradMixNorm);
	double myFact = (theX*myDiffLogDensity + 1 / myStd);
		for (uint i = 0; i < 3; i++)
		{
			theGrad[theBegIndex + i] = myMixNormGradSigma[i] * myFact + myGradMixNorm[i];
		}
	}

	/*!
	 * \fn void cMixNormResiduals::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData)
	 * \brief Compute the derivative of log density with respect to the random variable (theGradData[0]) \e and the gradient
	 * of log density with respect to the model parameters (other components in theGradData)
	 * \param theDate uint: time at which gradient is computed
	 * \param theValue const cRegArchValue&: value of the random variable
	 * \param theGradData cRegArchGradient&: concatenation of derivatives with respect to the random variable and the model parameters
	 */
	void cMixNormResiduals::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData) const
	{
	uint myBegIndex = theGradData.GetNMeanParam() + theGradData.GetNVarParam();
		GradLogDensity(theValue.mEpst[theDate], theGradData.mCurrentGradLogDens, *mDistrParameter, myBegIndex);
		theGradData.mCurrentDiffLogDensity = DiffLogDensity(theValue.mEpst[theDate]);
	}

	void cMixNormResiduals::RegArchParamToVector(cDVector& theDestVect, uint theIndex) const
	{
		if (theDestVect.GetSize() < theIndex + 3)
			throw cError("Wrong size");
		for (uint i = 0; i < 3; i++)
			theDestVect[theIndex + i] = (*mDistrParameter)[i];
	}
	
	void cMixNormResiduals::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
		if (theIndex > theSrcVect.GetSize()-3)
			throw cError("Wrong size");
		for (uint i = 0; i < 3; i++)
			(*mDistrParameter)[i] = theSrcVect[theIndex+i];

	}

	double cMixNormResiduals::ComputeEspAbsEps(void)
	{
	double myp = (*mDistrParameter)[0];
	double mySigma1 = (*mDistrParameter)[1];
	double mySigma2 = (*mDistrParameter)[2];
	double myStd = MixNormSigma(*mDistrParameter);

		return 2 * (myp*mySigma1 + (1 - myp)*mySigma2) / (SQRT_2_PI*myStd);
	
	}

	void cMixNormResiduals::ComputeGradBetaEspAbsEps(cDVector& theGrad)
	{
	double myp = (*mDistrParameter)[0];
	double mySigma1 = (*mDistrParameter)[1];
	double mySigma2 = (*mDistrParameter)[2];
	double myStd = MixNormSigma(*mDistrParameter);
	double myDen = SQRT_2_PI * pow(myStd, 3);
	double myDifSig = mySigma2 - mySigma1;

		theGrad[0] = myDifSig * myDifSig*(myp*mySigma1 - (1 - myp)*mySigma2) / myDen;
		theGrad[1] = 2 * (1 - myp) * myp * mySigma2 * myDifSig / myDen;
		theGrad[2] = -2 * (1 - myp) * myp * mySigma1 * myDifSig / myDen;
	}

	void cMixNormResiduals::ComputeHessBetaEspAbsEps(cDMatrix &theHess)
	{
	double myp = (*mDistrParameter)[0];
	double myq = 1 - myp;
	double mySigma1 = (*mDistrParameter)[1];
	double myVar1 = mySigma1 * mySigma1;
	double mySigma2 = (*mDistrParameter)[2];
	double myVar2 = mySigma2 * mySigma2;
	double myStd = MixNormSigma(*mDistrParameter);
	double myStd2 = myStd * myStd;
	double myStd3 = myStd2 * myStd;
	double myDen = SQRT_2_PI * pow(myStd, 3);
	double myDifSig = mySigma2 - mySigma1;
	double myDifSig2 = myDifSig * myDifSig;
	double myNum1 = myDifSig2 * (myp * mySigma1 - myq * mySigma2);
	double mydNum1p = 2 * myDifSig2 * (mySigma1 + mySigma2);
	double mydDenp = 3 * SQRT_2_PI * (myVar1-myVar2)* myStd / 2;
	double myDenc = myDen * myDen;
		theHess[0][0] = (mydNum1p * myDen - myNum1 * mydDenp) / myDenc;
	double mydNum11 = -2 * myDifSig * (myp * mySigma2 - 2 * mySigma2 + 3 * myp * mySigma1);
	double mydDen1 = 3 * SQRT_2_PI * myp * mySigma1 * myStd;
		theHess[0][1] = theHess[1][0] = (mydNum11 * myDen - myNum1 * mydDen1) / myDenc;
	double mydDen2 = 3 * SQRT_2_PI * myq * mySigma2 * myStd;
	double mydNum12 = 2 * myDifSig * ((myp + 1) * mySigma1 - 3 * myq * mySigma2) ;
		theHess[0][2] = theHess[2][0] = (mydNum12 * myDen - myNum1 * mydDen2) / myDenc;
	double myNum2 = 2 * myp * myq * myDifSig * mySigma2 ;
	double mydNum21 = -2 * myp * myq *mySigma2;
		theHess[1][1] = (mydNum21 * myDen - myNum2 * mydDen1) / myDenc;
	double mydNum22 = 2 * myp * myq * (2 * mySigma2 - mySigma1);
		theHess[1][2] = theHess[2][1] = (mydNum22 * myDen - myNum2 * mydDen2) / myDenc;
	double myNum3 = 2 * myp * myq * myDifSig * mySigma1;
	double mydNum32 = -2 * myp * myq * (mySigma2 - 2 * mySigma1);
		theHess[2][2] = (mydNum32 * myDen - myNum3 * mydDen2) / myDenc;
	}

	double cMixNormResiduals::Diff2LogDensity(double theX) const
	{
	double myp = (*mDistrParameter)[0];
	double mySigma1 = (*mDistrParameter)[1];
	double mySigma2 = (*mDistrParameter)[2];
	double myStd = MixNormSigma(*mDistrParameter);

		return myStd * myStd*MixNormDiff2LogDensity(theX*myStd, myp, mySigma1, mySigma2);
	}

	void cMixNormResiduals::GradDiffLogDensity(double theX, const cDVector& theDistrParam, cDVector& theGrad)
	{
		GradTrueDiffLogDensity(theX, theDistrParam, MixNormDiffLogDensity, MixNormDiff2LogDensity, MixNormSigma, MixNormGradSigma, MixNormGradDiffLogDens, theGrad);
	}

	static void HessLogDensity(double theX, cDMatrix& theHess, const cDVector& theDistrParam, uint theBegIndex)
	{
/*
double myp = theDistrParam[0];
	double mySigma1 = theDistrParam[1];
	double mySigma2 = theDistrParam[2];
	double myStd = MixNormSigma(theDistrParam);
	double myStd2 = myStd * myStd;
	double myX = theX * myStd;
	cDVector myGradStd(3);
		MixNormGradSigma(myp, mySigma1, mySigma2, myGradStd);
	cDMatrix myHessSigma(3, 3);
		MixNormHessSigma(myp, mySigma1, mySigma2, myHessSigma);
	cDMatrix myHess(3, 3);
		myHess = (theX * theX*MixNormDiff2LogDensity(myX, myp, mySigma1, mySigma2)+1/myStd2)*myGradStd*Transpose(myGradStd);
		myHess += (theX*MixNormDiffLogDensity(myX, myp, mySigma1, mySigma2) + 1 / myStd)*myHessSigma;
	cDMatrix myHessLogDens(3, 3);
		MixNormHessLogDensity(myX, myp, mySigma1, mySigma2, myHessLogDens);
		myHess += myHessLogDens;
	cDVector myGradDiffLog(3);
		MixNormGradLogDensity(myX, myp, mySigma1, mySigma2, myGradDiffLog);
		myHess += theX * (myGradDiffLog*Transpose(myGradStd) + myGradStd * Transpose(myGradDiffLog));
		for (int i = 0; i < 3; i++)
			for (int j = i; j < 3; j++)
				theHess[theBegIndex + i][theBegIndex + j] = theHess[theBegIndex + j][theBegIndex + i] = myHess[i][j];
*/
	cDMatrix myHess(3, 3);
		HessTrueLogDensity(theX, theDistrParam, &MixNormDiffLogDensity, &MixNormDiff2LogDensity, &MixNormSigma, &MixNormGradSigma, &MixNormGradDiffLogDens, &MixNormHessSigma, &MixNormHessLogDensity, myHess);
		for (int i = 0; i < 3; i++)
			for (int j = i; j < 3; j++)
				theHess[theBegIndex + i][theBegIndex + j] = theHess[theBegIndex + j][theBegIndex + i] = myHess[i][j];

	}

	void cMixNormResiduals::ComputeHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	double myX = theData.mEpst[theDate];
	cDVector myGradDiffLogDens = cDVector(3);
		GradDiffLogDensity(myX, *mDistrParameter, myGradDiffLogDens);
	uint myBegIndex = theHessData.GetNMeanParam() + theHessData.GetNVarParam();
		for (int i = 0 ; i < 3 ; i++)
			theHessData.mCurrentGradDiffLogDensity[myBegIndex+i] = myGradDiffLogDens[i];
		HessLogDensity(myX, theHessData.mCurrentHessDens, *mDistrParameter, myBegIndex);
	}

	void cMixNormResiduals::GetParamName(uint theIndex, char** theName)
	{
	uint myIndex = theIndex;
		sprintf(theName[myIndex++], "Bernouilli parameter p");
		sprintf(theName[myIndex++], "Sigma1");
		sprintf(theName[myIndex++], "Sigma2");
	}

	void  cMixNormResiduals::GetParamName(uint theIndex, string theName[])
	{
	uint myIndex = theIndex;
		char myChar[100];
		sprintf(myChar, "Bernouilli parameter p");
		theName[myIndex++] = myChar;
		sprintf(myChar, "Sigma1");
		theName[myIndex++] = myChar;
		sprintf(myChar, "Sigma2");
		theName[myIndex++] = myChar;
	}

}//namespace
