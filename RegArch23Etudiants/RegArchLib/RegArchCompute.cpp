#include "StdAfxRegArchLib.h"
#include <gsl/gsl_cdf.h>
/*!
	\file RegArchCompute.cpp
	\brief implementation of the Simulation / Estimation procedures for general RegArchModel
	
	\author Jean-Baptiste DURAND, Ollivier TARAMASCO 
	\date dec-18-2006 - last change feb-18-2011
*/


namespace RegArchLib {
	using namespace VectorAndMatrixNameSpace;


	/*!
	 * \fn void RegArchSimul(const uint theNSample, const cRegArchModel& theModel, cRegArchValue& theData)
	 * \param const uint theNSample: size of the sample
	 * \param const cRegArchModel& theModel: the RegArch model
	 * \param cRegArchValue& theData: output parameter. The Y(t) values are stored in theData.mYt
	 */
	void RegArchSimul(const uint theNSample, const cRegArchModel& theModel, cRegArchValue& theData)
	{
#ifdef _CLI_
//		cout << "Debut de RegArchSimul" << endl;
#endif //_CLI_

		theData.ReAlloc(theNSample) ;
		theModel.mResids->Generate(theNSample, theData.mEpst) ;
#ifdef _CLI_
//		cout << "lo1" << endl;
#endif //_CLI_
		theModel.mVar->UpdateProxyVarParameters();
#ifdef _CLI_
//		cout << "lo2" << endl;
#endif //_CLI_
		if (theModel.mMean != NULL)
			theModel.mMean->UpdateProxyMeanParameters();
#ifdef _CLI_
//		cout << "lo3" << endl;
#endif //_CLI_
		for (uint t = 0 ; t < theNSample ; t++)
		{	theData.mHt[t] = theModel.mVar->ComputeVar(t, theData) ;
			if (theModel.mMean != NULL)
				theData.mMt[t] = theModel.mMean->ComputeMean(t, theData) ;
			theData.mUt[t] = sqrt(theData.mHt[t])*theData.mEpst[t] ;
			theData.mYt[t] = theData.mMt[t] + theData.mUt[t] ;
		}
	}

   /*!
	 * \fn void RegArchSimul(const uint theNSample, const cRegArchModel& theModel, cDVector& theYt)
	 * \param const uint theNSample: size of the sample
	 * \param const cRegArchModel& theModel: the RegArch model
	 * \param cDVector& theYt: output parameter. 
	 */
	void RegArchSimul(const uint theNSample, const cRegArchModel& theModel, cDVector& theYt, cDMatrix* theXt, cDMatrix* theXvt)
	{
	cRegArchValue myValue(theNSample, theXt, theXvt) ;
		RegArchSimul(theNSample, theModel, myValue) ;
		theYt = myValue.mYt ;

	}

	void FillValue(uint theDate, const cRegArchModel& theModel, cRegArchValue& theValue)
	{
		theValue.mHt[theDate] = theModel.mVar->ComputeVar(theDate, theValue);
		if (theModel.mMean != NULL)
			theValue.mMt[theDate] = theModel.mMean->ComputeMean(theDate, theValue);
		theValue.mUt[theDate] = theValue.mYt[theDate] - theValue.mMt[theDate];
		theValue.mEpst[theDate] = theValue.mUt[theDate] / sqrt(theValue.mHt[theDate]);
	}

	void FillValueForNumericGrad(uint theDate, const cRegArchModel& theModel, cRegArchValue& theValue, cNumericDerivative& theNumDeriv)
	{
		FillValue(theDate, theModel, theValue);
	cRegArchModel myModel(theModel);
	uint myNParam = theModel.GetNParam();
	uint myNDistrParam = myModel.mResids->GetNParam();
	uint myNMeanVarParam = myNParam - myNDistrParam;
	cDVector myVectParam(myNParam);
	cEgarch* myEgarch;
		if (theModel.mVar->GetCondVarType() == eEgarch)
		{
			myEgarch = static_cast<cEgarch*>(myModel.mVar);
		}
	double myh = theNumDeriv.Geth();
		myModel.RegArchParamToVector(myVectParam);
		for (uint j = 0; j < myNParam; j++)
		{
			theNumDeriv.mh1[j] = fabs(myh*myVectParam[j]);
			if (theNumDeriv.mh1[j] < 1e-16)
				theNumDeriv.mh1[j] = fabs(myh);
			myVectParam[j] += theNumDeriv.mh1[j];
			myModel.VectorToRegArchParam(myVectParam);
			myModel.mVar->UpdateProxyVarParameters();
			if (myModel.mMean != NULL)
				myModel.mMean->UpdateProxyMeanParameters();
			if (j >= myNMeanVarParam)
			{
				if (myModel.mVar->GetCondVarType() == eEgarch)
				{
					myEgarch->SetEspAbsEps(myModel.mResids->ComputeEspAbsEps());
					myModel.SetVar(*myEgarch);
				}
			}
			FillValue(theDate, myModel, theNumDeriv.mValueForGrad[j]);
			myVectParam[j] -= theNumDeriv.mh1[j];
		}
		if (myNDistrParam > 0)
		{
			for (uint j = 0; j < myNDistrParam; j++)
			{
			uint k = myNMeanVarParam + j;
				myVectParam[k] += theNumDeriv.mh1[k];
				myModel.VectorToRegArchParam(myVectParam);
				if (myModel.mVar->GetCondVarType() == eEgarch)
				{
					myEgarch->SetEspAbsEps(myModel.mResids->ComputeEspAbsEps());
					myModel.SetVar(*myEgarch);
				}
				theNumDeriv.mLogDensForGrad[j] = myModel.mResids->LogDensity(theValue.mEpst[theDate]);
				myVectParam[j] -= theNumDeriv.mh1[k];

			}
		}
					
	}

	void FillValueForNumericGradAndHess(uint theDate, const cRegArchModel& theModel, cRegArchValue& theValue, cNumericDerivative& theNumDeriv)
	{
		FillValue(theDate, theModel, theValue);
	cRegArchModel myModel(theModel);
	uint myNParam = theModel.GetNParam();
	uint myNDistrParam = myModel.mResids->GetNParam();
	uint myNMeanVarParam = myNParam - myNDistrParam;
	cDVector myVectParam(myNParam);
	double myh = theNumDeriv.Geth();
	cEgarch* myEgarch=NULL;
		if (theModel.mVar->GetCondVarType() == eEgarch)
		{
		cEgarch* myEgarch1 = static_cast<cEgarch *>(theModel.mVar);
			myEgarch = new cEgarch(*myEgarch1);
		}
		myModel.RegArchParamToVector(myVectParam);
		for (uint j = 0; j < myNParam; j++)
		{
		double myh1 = fabs(myh*myVectParam[j]);
			if (myh1 < 1e-16)
				myh1 = fabs(myh);
			myVectParam[j] += myh1;
			theNumDeriv.mh1[j] = myh1;
			myModel.VectorToRegArchParam(myVectParam);
			if (j >= myNMeanVarParam)
			{
				if (myModel.mVar->GetCondVarType() == eEgarch)
				{
					myEgarch->SetEspAbsEps(myModel.mResids->ComputeEspAbsEps());
					myModel.SetVar(*myEgarch);
				}
			}
			myModel.mVar->UpdateProxyVarParameters();
			if (myModel.mMean != NULL)
				myModel.mMean->UpdateProxyMeanParameters();
			FillValue(theDate, myModel, theNumDeriv.mValueForGrad[j]);
			for (uint k = j; k < myNParam; k++)  
			{
			double myh2;
				if (k == j)
					myh2 = myh1 ;
				else
				{
					myh2 = fabs(myh*myVectParam[k]);
					if (myh2 < 1e-16)
						myh2 = fabs(myh);
				}
				theNumDeriv.mh2[j][k] = theNumDeriv.mh2[k][j] = myh2;
				myVectParam[k] += myh2;
				myModel.VectorToRegArchParam(myVectParam);
				if (k >= myNMeanVarParam)
				{
					if (myModel.mVar->GetCondVarType() == eEgarch)
					{
						myEgarch->SetEspAbsEps(myModel.mResids->ComputeEspAbsEps());
						myModel.SetVar(*myEgarch);
					}
				}
				myModel.mVar->UpdateProxyVarParameters();
				if (myModel.mMean != NULL)
					myModel.mMean->UpdateProxyMeanParameters();
				FillValue(theDate, myModel, theNumDeriv.mValueForHess[j][k]);
//				theNumDeriv.mValueForHess[k][j].mEpst[theDate] = theNumDeriv.mValueForHess[j][k].mEpst[theDate];
//				theNumDeriv.mValueForHess[k][j].mMt[theDate] = theNumDeriv.mValueForHess[j][k].mMt[theDate];
//				theNumDeriv.mValueForHess[k][j].mHt[theDate] = theNumDeriv.mValueForHess[j][k].mHt[theDate];
//				theNumDeriv.mValueForHess[k][j].mUt[theDate] = theNumDeriv.mValueForHess[j][k].mUt[theDate];
				myVectParam[k] -= myh2;
			}
			myVectParam[j] -= myh1;
		}
		if (myNDistrParam > 0)
		{
			for (uint j = 0; j < myNDistrParam; j++)
			{
			uint k = myNMeanVarParam + j;
				myVectParam[k] += theNumDeriv.mh1[k];
				myModel.VectorToRegArchParam(myVectParam);
				if (myModel.mVar->GetCondVarType() == eEgarch)
				{
					myEgarch->SetEspAbsEps(myModel.mResids->ComputeEspAbsEps());
					myModel.SetVar(*myEgarch);
				}
				theNumDeriv.mLogDensForGrad[j] = myModel.mResids->LogDensity(theValue.mEpst[theDate]);
				for (uint i = 0; i < myNParam; i++)
				{
					myVectParam[i] += theNumDeriv.mh2[k][i];
					myModel.VectorToRegArchParam(myVectParam);
					if (i >= myNMeanVarParam)
					{
						if (myModel.mVar->GetCondVarType() == eEgarch)
						{
							myEgarch->SetEspAbsEps(myModel.mResids->ComputeEspAbsEps());
							myModel.SetVar(*myEgarch);
						}
					}
					theNumDeriv.mLogDensForHess[j][i] = myModel.mResids->LogDensity(theValue.mEpst[theDate]);
					myVectParam[i] -= theNumDeriv.mh2[k][i];
					myModel.VectorToRegArchParam(myVectParam);
				}
			double myh3 = fabs(myh*theValue.mEpst[theDate]);
				if (myh3 < 1e-16)
					myh3 = myh;
				theNumDeriv.mGradDiffForGrad[j] = myModel.mResids->LogDensity(theValue.mEpst[theDate] + myh3);
				myVectParam[k] -= theNumDeriv.mh1[k];

			}
		}
		if (myEgarch != NULL)
			delete myEgarch;
	}

   /*!
	* \fn static double ComputeLt(uint theDate, const cRegArchModel& theModel, cRegArchValue& theValue)
	* \param uint theDate: the date
	* \param const cRegArchModel& theModel: the model
	* \param cRegArchValue& theValue: output parameter.
	* \brief return the lt(theta ; Yt)
	*/
	static double ComputeLt(uint theDate, const cRegArchModel& theModel, cRegArchValue& theValue)
	{
		FillValue(theDate, theModel, theValue);
	double myRes  = -0.5*log(theValue.mHt[theDate]) + theModel.mResids->LogDensity(theValue.mEpst[theDate]);
		return myRes;
	}

	/*!
	 * \fn double RegArchLLH(const cRegArchModel& theModel, cDVector* theYt, cDMatrix* theXt)
	 * \param const cRegArchModel& theModel: the model
	 * \param cDVector* theYt: the observations
	 * \param cDMatrix* theXt: the regressors matrix. Default NULL
	 * \brief return the log-likelihood value
	 */
	double RegArchLLH(const cRegArchModel& theModel, cDVector* theYt, cDMatrix* theXt)
	{
	cRegArchValue myValue(theYt, theXt) ;
		return RegArchLLH(theModel, myValue) ;
	}

	/*!
	 * \fn double RegArchLLH(const cRegArchModel& theModel,cRegArchValue& theData)
	 * \param const cRegArchModel& theModel: the model
	 * \param cRegArchValue& theData: theData.mYt contains the observations.
	 * \brief return the log-likelihood value
	 */
	double RegArchLLH(const cRegArchModel& theModel, cRegArchValue& theData)
	{
	int mySize = (int)theData.mYt.GetSize() ;
	double myRes = 0 ;
		theData.mEpst = theData.mHt = theData.mMt = theData.mUt = 0.0 ;
		theModel.mVar->UpdateProxyVarParameters();
		if (theModel.mMean != NULL)
			theModel.mMean->UpdateProxyMeanParameters();
		for(int t=0 ; t < mySize ; t++)
		{
			FillValue(t, theModel, theData);
			myRes += -0.5*log(theData.mHt[t]) + theModel.mResids->LogDensity(theData.mEpst[t]) ;
		}
		return myRes ;
	}

	void RegArchGradLt(int theDate, cRegArchModel& theModel, cRegArchValue& theValue, cRegArchGradient& theGradData, cDVector& theGradlt)
	{	theModel.ComputeGrad(theDate, theValue, theGradData);
	double mySigmat = sqrt(theValue.mHt[theDate]);
		theGradlt = (-1.0 / mySigmat) * theGradData.mCurrentGradSigma + theGradData.mCurrentDiffLogDensity * theGradData.mCurrentGradEps + theGradData.mCurrentGradLogDens;
	}
	
	void NumericRegArchGradLt(uint theDate, cRegArchModel& theModel, cRegArchValue* theValue, cDVector& theGradlt, double theh)
	{
	uint myNParam = theModel.GetNParam();
	cDVector myParam(myNParam);
	cDVector myParam1(myNParam);
		theModel.RegArchParamToVector(myParam);
	cRegArchModel myModel = cRegArchModel(theModel) ;
		theGradlt.ReAlloc(myNParam);
	double mylt0 = ComputeLt(theDate, theModel, theValue[myNParam]);
		myModel.VectorToRegArchParam(myParam);
		myParam1 = myParam;
		for (uint i = 0; i < myNParam; i++)
		{
		double myh1 = abs(myParam[i] * theh);
		myParam1[i] += myh1;
		myModel.VectorToRegArchParam(myParam1);
		double mylt1 = ComputeLt(theDate, myModel, theValue[i]);
		theGradlt[i] = (mylt1 - mylt0) / myh1;
		myParam1[i] -= myh1;
		}
	}

	void RegArchLtAndGradLt(int theDate, cRegArchModel& theModel, cRegArchValue& theValue, cRegArchGradient& theGradData, double& theLt, cDVector& theGradlt)
	{
		theModel.ComputeGrad(theDate, theValue, theGradData);
	double mySigmat = sqrt(theValue.mHt[theDate]);
		theGradlt = (-1.0 / mySigmat) * theGradData.mCurrentGradSigma + theGradData.mCurrentDiffLogDensity * theGradData.mCurrentGradEps + theGradData.mCurrentGradLogDens;
		theLt = -0.5*log(theValue.mHt[theDate]) + theModel.mResids->LogDensity(theValue.mEpst[theDate]) ;
	}

	void RegArchGradLLH(cRegArchModel& theModel, cRegArchValue& theData, cDVector& theGradLLH)
	{
	cRegArchGradient myGradData=cRegArchGradient(&theModel) ;
	cDVector myGradlt(myGradData.GetNParam()) ;
		theGradLLH = 0.0L ;
		theModel.mVar->UpdateProxyVarParameters();
		if (theModel.mMean != NULL)
			theModel.mMean->UpdateProxyMeanParameters();
		for (int t = 0 ; t < (int)theData.mYt.GetSize() ; t++)
		{	RegArchGradLt(t, theModel, theData, myGradData, myGradlt) ;
			theGradLLH += myGradlt ;
			myGradData.Update();
		}
	}

	void RegArchLLHAndGradLLH(cRegArchModel& theModel, cRegArchValue& theValue, double& theLLH, cDVector& theGradLLH)
	{
	cRegArchGradient myGradData(&theModel) ;
	cDVector myGradlt(myGradData.GetNParam()) ;
	double myLt ;
		theGradLLH = 0.0L ;
		theLLH = 0.0 ;
		theModel.mVar->UpdateProxyVarParameters();
		if (theModel.mMean != NULL)
			theModel.mMean->UpdateProxyMeanParameters();
		theModel.mVar->UpdateProxyVarParameters();
		for (int t = 0 ; t < (int)theValue.mYt.GetSize() ; t++)
		{	RegArchLtAndGradLt(t, theModel, theValue, myGradData, myLt, myGradlt) ;
			theGradLLH += myGradlt ;
			theLLH += myLt ;
			myGradData.Update();
		}
	}

	void RegArchHessLt(int theDate, cRegArchModel& theModel, cRegArchValue& theValue, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cDMatrix& theHesslt)
	{
	cDVector myGradlt(theGradData.GetNParam());
		theHesslt = 0.0;
		RegArchGradLt(theDate, theModel, theValue, theGradData, myGradlt);
		theModel.mVar->ComputeHess(theDate, theValue, theGradData, theHessData, theModel.mResids);
		theHessData.ComputeHessSigmatFromHessVart(theDate, theValue, theGradData);
		if (theModel.mMean != NULL)
			theModel.mMean->ComputeHess(theDate, theValue, theGradData, theHessData);
		theModel.mResids->ComputeHess(theDate, theValue, theGradData, theHessData);

	double myHt = theValue.mHt[theDate];
	double mySigmat = sqrt(myHt);
		/* Hessien Eps(t, theta) */
		theHessData.mCurrentHessEps = theGradData.mCurrentGradSigma * VectorAndMatrixNameSpace::Transpose(theGradData.mCurrentGradMu) / myHt;
		theHessData.mCurrentHessEps += theGradData.mCurrentGradMu * Transpose(theGradData.mCurrentGradSigma) / myHt;
		theHessData.mCurrentHessEps -= theHessData.mCurrentHessMu / mySigmat;
		theHessData.mCurrentHessEps += 2 * theValue.mUt[theDate] / (myHt * mySigmat) * theGradData.mCurrentGradSigma*Transpose(theGradData.mCurrentGradSigma);
		theHessData.mCurrentHessEps -= theValue.mUt[theDate] / myHt *theHessData.mCurrentHessSigma;
		/* Hessien lt*/
		theHesslt = theGradData.mCurrentGradSigma * Transpose(theGradData.mCurrentGradSigma) / myHt;
		theHesslt -= theHessData.mCurrentHessSigma / mySigmat;
	double myDiff2LogDensity = theModel.mResids->Diff2LogDensity(theValue.mEpst[theDate]);
		theHesslt += theGradData.mCurrentGradEps * Transpose(theGradData.mCurrentGradEps)*myDiff2LogDensity;
		theHesslt += theGradData.mCurrentDiffLogDensity * theHessData.mCurrentHessEps;
		theHesslt += theHessData.mCurrentGradDiffLogDensity * Transpose(theGradData.mCurrentGradEps);
		theHesslt += theGradData.mCurrentGradEps * Transpose(theHessData.mCurrentGradDiffLogDensity);
		theHesslt += theHessData.mCurrentHessDens;
	}

	void NumericRegArchHessLt(int theDate, cRegArchModel& theModel, cRegArchValue* theValue, cRegArchGradient* theGradData, cDMatrix& theHesslt, double theh)
	{
	uint myNParam = theModel.GetNParam();
	cDVector myGradLt0 = cDVector(myNParam);
	cDVector myGradLt1 = cDVector(myNParam);
	cDVector myParam(myNParam);
		theModel.RegArchParamToVector(myParam);
		RegArchGradLt(theDate, theModel, theValue[myNParam], theGradData[myNParam], myGradLt0);
		theGradData[myNParam].Update();
		for (uint i = 0; i < myNParam; i++)
		{
		double myh = abs(myParam[i] * theh);
			myParam[i] += myh;
			theModel.VectorToRegArchParam(myParam);
			RegArchGradLt(theDate, theModel, theValue[i], theGradData[i], myGradLt1);
			theGradData[i].Update();
			for (uint j = i; j < myNParam; j++)
			theHesslt[i][j] = theHesslt[j][i] = (myGradLt1[j] - myGradLt0[j]) / myh;
				myParam[i] -= myh;
		}
		theModel.VectorToRegArchParam(myParam);
	}

	void RegArchGradAndHessLt(int theDate, cRegArchModel& theModel, cRegArchValue& theValue, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cDVector& theGradlt, cDMatrix& theHesslt)
	{
	double mylt = 0.0;
		theGradlt.ReAlloc(theModel.GetNParam());
		theHesslt = 0.0;
		RegArchLtAndGradLt(theDate, theModel, theValue, theGradData, mylt, theGradlt);
		theModel.mVar->ComputeHess(theDate, theValue, theGradData, theHessData, theModel.mResids);
		theHessData.ComputeHessSigmatFromHessVart(theDate, theValue, theGradData);
		if (theModel.mMean != NULL)
			theModel.mMean->ComputeHess(theDate, theValue, theGradData, theHessData);
		theModel.mResids->ComputeHess(theDate, theValue, theGradData, theHessData);
	double mySigmat = sqrt(theValue.mHt[theDate]);
	/* Hessien Eps(t, theta) */
		theHessData.mCurrentHessEps = theGradData.mCurrentGradSigma * Transpose(theGradData.mCurrentGradMu) / theValue.mHt[theDate];
		theHessData.mCurrentHessEps += theGradData.mCurrentGradMu * Transpose(theGradData.mCurrentGradSigma) / theValue.mHt[theDate];
		theHessData.mCurrentHessEps -= theHessData.mCurrentHessMu / mySigmat;
		theHessData.mCurrentHessEps += 2 * theValue.mUt[theDate] / pow(mySigmat, 3)*theGradData.mCurrentGradSigma*Transpose(theGradData.mCurrentGradSigma);
		theHessData.mCurrentHessEps -= theValue.mUt[theDate] / theValue.mHt[theDate] * theHessData.mCurrentHessSigma;
	/* Hessien lt*/
		theHesslt = theGradData.mCurrentGradSigma * Transpose(theGradData.mCurrentGradSigma) / theValue.mHt[theDate];
		theHesslt -= theHessData.mCurrentHessSigma / mySigmat;
	double myDiff2LogDensity = theModel.mResids->Diff2LogDensity(theValue.mEpst[theDate]);
		theHesslt += theGradData.mCurrentGradEps * Transpose(theGradData.mCurrentGradEps)*myDiff2LogDensity;
		theHesslt += theGradData.mCurrentDiffLogDensity * theHessData.mCurrentHessEps;
		theHesslt += theHessData.mCurrentGradDiffLogDensity * Transpose(theGradData.mCurrentGradEps);
		theHesslt += theGradData.mCurrentGradEps * Transpose(theHessData.mCurrentGradDiffLogDensity);
		theHesslt += theHessData.mCurrentHessDens;
	}

	void RegArchLtGradAndHessLt(int theDate, cRegArchModel& theModel, cRegArchValue& theValue, double& thelt, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cDVector& theGradlt, cDMatrix& theHesslt)
	{
	uint myNParam = theModel.GetNParam();
		theGradlt.ReAlloc(myNParam);
		theHesslt.ReAlloc(myNParam, myNParam);
		RegArchLtAndGradLt(theDate, theModel, theValue, theGradData, thelt, theGradlt);
		theModel.mVar->ComputeHess(theDate, theValue, theGradData, theHessData, theModel.mResids);
		if (theModel.mMean != NULL)
			theModel.mMean->ComputeHess(theDate, theValue, theGradData, theHessData);
		theModel.mResids->ComputeHess(theDate, theValue, theGradData, theHessData);
		theHessData.ComputeHessSigmatFromHessVart(theDate, theValue, theGradData);
	double mySigmat = sqrt(theValue.mHt[theDate]);
		/* Hessien Eps(t, theta) */
		theHessData.mCurrentHessEps = theGradData.mCurrentGradSigma * Transpose(theGradData.mCurrentGradMu) / theValue.mHt[theDate];
		theHessData.mCurrentHessEps += theGradData.mCurrentGradMu * Transpose(theGradData.mCurrentGradSigma) / theValue.mHt[theDate];
		theHessData.mCurrentHessEps -= theHessData.mCurrentHessMu / mySigmat;
		theHessData.mCurrentHessEps += 2 * theValue.mUt[theDate] / pow(mySigmat, 3)*theGradData.mCurrentGradSigma*Transpose(theGradData.mCurrentGradSigma);
		theHessData.mCurrentHessEps -= theValue.mUt[theDate] / theValue.mHt[theDate] * theHessData.mCurrentHessSigma;
		/* Hessien lt*/
		theHesslt = theGradData.mCurrentGradSigma * Transpose(theGradData.mCurrentGradSigma) / theValue.mHt[theDate];
		theHesslt -= theHessData.mCurrentHessSigma / mySigmat;
	double myDiff2LogDensity = theModel.mResids->Diff2LogDensity(theValue.mEpst[theDate]);
		theHesslt += theGradData.mCurrentGradEps * Transpose(theGradData.mCurrentGradEps)*myDiff2LogDensity;
		theHesslt += theGradData.mCurrentDiffLogDensity * theHessData.mCurrentHessEps;
		theHesslt += theHessData.mCurrentGradDiffLogDensity * Transpose(theGradData.mCurrentGradEps);
		theHesslt += theGradData.mCurrentGradEps * Transpose(theHessData.mCurrentGradDiffLogDensity);
		theHesslt += theHessData.mCurrentHessDens;
	}

	void RegArchHessLLH(cRegArchModel& theModel, cRegArchValue& theValue, cDMatrix& theHessLLH)
	{
	cRegArchGradient myGradData = cRegArchGradient(&theModel);
	cRegArchHessien myHessData = cRegArchHessien(&theModel);
	uint myNParam = myHessData.GetNParam();
		cDMatrix myHesslt(myNParam, myNParam);
		theHessLLH.ReAlloc(myNParam, myNParam);
		theModel.mVar->UpdateProxyVarParameters();
		if (theModel.mMean != NULL)
			theModel.mMean->UpdateProxyMeanParameters();
		for (int t = 0; t < (int)theValue.mYt.GetSize(); t++)
		{
			RegArchHessLt(t, theModel, theValue, myGradData, myHessData, myHesslt);
			theHessLLH += myHesslt;
			myGradData.Update();
			myHessData.Update();
		}
	}


	
#ifdef _LOG_


#endif // _LOG_
}
