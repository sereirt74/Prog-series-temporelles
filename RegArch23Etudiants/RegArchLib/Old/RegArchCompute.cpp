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
		
		theData.ReAlloc(theNSample) ;
		theModel.mResids->Generate(theNSample, theData.mEpst) ;
		theModel.mVar->UpdateProxyVarParameters();
		if (theModel.mMean != NULL)
			theModel.mMean->UpdateProxyMeanParameters();
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
			theModel.mMean->ComputeHess(theDate, theValue, theGradData, theHessData, theModel.mResids);
		theModel.mResids->ComputeHess(theDate, theValue, theGradData, theHessData, theModel.mResids);

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
			theModel.mMean->ComputeHess(theDate, theValue, theGradData, theHessData, theModel.mResids);
		theModel.mResids->ComputeHess(theDate, theValue, theGradData, theHessData, theModel.mResids);
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
			theModel.mMean->ComputeHess(theDate, theValue, theGradData, theHessData, theModel.mResids);
		theModel.mResids->ComputeHess(theDate, theValue, theGradData, theHessData, theModel.mResids);
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

	void NumericComputeJ(cRegArchModel& theModel, cRegArchValue& theValue, cDMatrix& theJ, double theh=1e-5)
	{
	//on initialise l'ensemble des variables utiles
	uint myNParam = theModel.GetNParam();
	cDVector myGradij(myNParam), myVect(myNParam), myVect0(myNParam);
		theJ.ReAlloc(myNParam, myNParam);
	uint myT = theValue.mYt.GetSize();
	cRegArchGradient myGradData(&theModel);

	cDVector* myGrad0 = new cDVector[myT];

		for (uint t = 0; t < myT; t++)
		{
			myGrad0[t].ReAlloc(myNParam);
			RegArchGradLt(t, theModel, theValue, myGradData, myGrad0[t]);
		}

		//on met � 0 la hessienne
		theJ = 0.0L;

		theModel.RegArchParamToVector(myVect0);

		myVect = myVect0;
	double myhhi;
		for (uint i = 0; i < myNParam; i++)
		{
		//on fait la somme sur T
			myhhi = fabs(theh*myVect0[i]);
			if (myhhi < 1e-16)
				myhhi = theh;
			myVect[i] += myhhi;
			theModel.VectorToRegArchParam(myVect);
			if (theModel.mMean != NULL)
				theModel.mMean->UpdateProxyMeanParameters();
			theModel.mVar->UpdateProxyVarParameters();
			myGradData.ReInitialize();
			for (uint t = 0; t < myT; t++)
			{
				RegArchGradLt(t, theModel, theValue, myGradData, myGradij);
				for (uint j = i; j < myNParam; j++)
				{
					//on somme dans theHessLLH
					theJ[i][j] += (myGradij[j] - myGrad0[t][j]) / myhhi;
				}
			}
			myVect[i] -= myhhi;
		}
	//on divise par T
		for (uint i = 0; i < myNParam; i++)
			for (uint j = i + 1; j < myNParam; j++)
				theJ[j][i] = theJ[i][j];
		theJ /= (double)myT;
		// On d�salloue
		for (uint t = 0; t < myT; t++)
			myGrad0[t].Delete();
		delete[] myGrad0;
	}

	void NumericComputeCov(cRegArchModel &theModel, cRegArchValue &theData, cDMatrix &theCov)
	{
	uint myNParam = theModel.GetNParam();
	theCov.ReAlloc(myNParam,myNParam);
	cDMatrix myI(myNParam, myNParam);  
	cDMatrix myJ(myNParam, myNParam);
	RegArchComputeI(theModel,theData,myI);
//	cout<< " I : " << endl;
	myI.Print();
//	cout << "J : " << endl;
	NumericComputeJ(theModel, theData, myJ);
	myJ.Print();

	theCov =  Inv(myJ) * myI * Inv(myJ);
	theCov = theCov/(int)theData.mYt.GetSize();
	}

	void RegArchComputeI(cRegArchModel &theModel,cRegArchValue &theData, cDMatrix &theI)
	{
	uint myNParam = theModel.GetNParam();

		theI.ReAlloc(myNParam, myNParam);
	cRegArchGradient myGradData(&theModel);
	cDVector myGradlt(myNParam);
	cRegArchGradient myGradlt2(myNParam);
	cDMatrix myGradltTranspose(1,myNParam);
	
		if (theModel.mMean != NULL)
			theModel.mMean->UpdateProxyMeanParameters();
		theModel.mVar->UpdateProxyVarParameters();
		for (int t = 0 ; t < (int)theData.mYt.GetSize() ; t++)
		{
			//Calcul du gradient de lt pour th�ta chapeau
			RegArchGradLt(t,theModel,theData,myGradData, myGradlt);
			myGradltTranspose = Transpose(myGradlt);
			theI += (myGradlt * myGradltTranspose);
		}
		theI /= (double)(theData.mYt.GetSize());
	}

	void RegArchComputeIAndJ(cRegArchModel &theModel,cRegArchValue &theData, cDMatrix &theI, cDMatrix &theJ)
	{
	uint myNParam = theModel.GetNParam();

	theI.ReAlloc(myNParam, myNParam);
	theJ.ReAlloc(myNParam, myNParam);
	cRegArchGradient myGradData(&theModel);
	cRegArchHessien myHessData(&theModel) ;

	cDVector myGradlt(myNParam);
	cDMatrix myHesslt(myNParam, myNParam) ;

	cRegArchGradient myGradlt2(myNParam);
	
	cDMatrix myGradltTranspose(1,myNParam);
	cDMatrix myGradI(myNParam,myNParam);
		if (theModel.mMean != NULL)
			theModel.mMean->UpdateProxyMeanParameters();
		theModel.mVar->UpdateProxyVarParameters();
		for (int t = 0 ; t < (int)theData.mYt.GetSize() ; t++)
		{
		//Calcul du gradient et du hessien de lt pour th�ta chapeau
			RegArchGradAndHessLt(t, theModel, theData, myGradData, myHessData, myGradlt, myHesslt) ;

			myGradltTranspose = Transpose(myGradlt);
			myGradI += (myGradlt * myGradltTranspose);
			theJ -= myHesslt ;
		}
		theI = myGradI / (double)theData.mYt.GetSize();
		theJ /= (double)theData.mYt.GetSize() ;
	}

	void RegArchComputeCov(cRegArchModel& theModel, cRegArchValue& theValue, cDMatrix& theCov)
	{
	uint myNParam = theModel.GetNParam() ;
	cDMatrix myI(myNParam, myNParam) ;
	cDMatrix myJ(myNParam, myNParam) ;
	cDMatrix myInvJ(myNParam, myNParam) ;

		RegArchComputeIAndJ(theModel, theValue, myI, myJ) ;
//		RegArchComputeI(theModel, theValue, myI);
//		NumericComputeJ(theModel, theValue, myJ);
		myInvJ = Inv(myJ) ;
		theCov = myInvJ * myI * myInvJ  / theValue.mYt.GetSize();
	}

	void RegArchComputeCov(cRegArchModel& theModel, cRegArchValue& theValue, cDMatrix& theCov, int& theError)
	{
	uint myNParam = theModel.GetNParam();
	cDMatrix myI(myNParam, myNParam);
	cDMatrix myJ(myNParam, myNParam);
	cDMatrix myInvJ(myNParam, myNParam);

		RegArchComputeIAndJ(theModel, theValue, myI, myJ);
		//		RegArchComputeI(theModel, theValue, myI);
		//		NumericComputeJ(theModel, theValue, myJ);
		myInvJ = Inv(myJ, theError);
		theCov = myInvJ * myI * myInvJ / theValue.mYt.GetSize();
	}

	void NumericRegArchGradLLH(cRegArchModel& theModel, cRegArchValue& theValue, cDVector& theGradLLH, double theh)
	{
	double myLLH0 = RegArchLLH(theModel, theValue) ;
	int myNParam = (int)theGradLLH.GetSize() ;
	int myNLawParam = (int)theModel.mResids->GetNParam() ;
	eCondVarEnum myVarType = theModel.mVar->GetCondVarType() ;

	cDVector myVectParam(myNParam), myVect0(myNParam) ;
	
		theModel.RegArchParamToVector(myVectParam) ;
		theModel.mVar->UpdateProxyVarParameters();
		myVect0 = myVectParam ;
		for (int i = 0 ; i < myNParam ; i++)
		{
		double myhh = fabs(theh * myVectParam[i]) ;
			if (myhh < 1e-16)
				myhh = theh ;
			myVectParam[i] += myhh ;
			theModel.VectorToRegArchParam(myVectParam) ;
			theModel.mVar->UpdateProxyVarParameters();
	
		double myLLH1 = RegArchLLH(theModel, theValue) ;
			theGradLLH[i] = (myLLH1 - myLLH0)/myhh ;
			myVectParam[i] -= myhh ;
		}
		theModel.VectorToRegArchParam(myVect0) ;
	}

	void NumericRegArchHessLLHold(cRegArchModel& theModel, cRegArchValue& theValue, cDMatrix& theHessLLH, double theh)
	{
	int myNParam = (int)theModel.GetNParam() ;
	cDVector myGrad0(myNParam), myGradij(myNParam), myVect(myNParam), myVect0(myNParam) ;
		theHessLLH.ReAlloc(myNParam, myNParam) ;
	
		theModel.RegArchParamToVector(myVect) ;
		myVect0 = myVect ;
		RegArchGradLLH(theModel, theValue, myGrad0) ; 	

		for (int i = 0 ; i < myNParam ; i++)
		{	
		double myhhi = fabs(theh*myGrad0[i]) ;
			if (myhhi < 1e-16)
				myhhi = theh ;
			myVect[i] += myhhi ;
			theModel.VectorToRegArchParam(myVect) ;
			RegArchGradLLH(theModel, theValue, myGradij) ;
			myVect[i] -= myhhi ;
			for (int j = 0 ; j < myNParam ; j++)
				theHessLLH[i][j] = (myGradij[j] - myGrad0[j])/myhhi ;
		}
		theModel.VectorToRegArchParam(myVect0) ;
		for (int i = 0 ; i < myNParam-1 ; i++)
			for (int j = i+1 ; j < myNParam ; j++)
				theHessLLH[i][j] = theHessLLH[j][i] = (theHessLLH[i][j] + theHessLLH[j][i])/2.0 ;
	}	

	void RegArchEstim(cRegArchModel& theModel, cRegArchValue& theValue, sGSLMultiMinResult& theResStruct, cRegArchModel& theResModel, cDVector* theInitPoint, eGSLMultiMinAlgoEnum theAlgo, double theStopValue, int theMaxIter, bool theVerbose)
	{
	uint myNParam = theModel.GetNParam() ;

		if (theModel.mMean != NULL)
			theResModel.SetMean(*(theModel.mMean)) ;
		else
			theResModel.mMean = NULL ;

		theResModel.SetVar(*(theModel.mVar)) ;
		theResModel.SetResid(*(theModel.mResids)) ;

		if (theInitPoint == NULL)
		{
			theResModel.SetDefaultInitPoint(theValue);		
		}
		else
		{	if (theInitPoint->GetGSLVector() == NULL)
				theResModel.SetDefaultInitPoint(theValue) ;
			else
				theResModel.VectorToRegArchParam(*theInitPoint) ;
		}
	cDVector myX0(myNParam) ;
		theResModel.RegArchParamToVector(myX0) ;
	cGSLMultiMin myMultiMin(myX0, theAlgo) ;

	sParamOptimStruct myOtherParam ;
		myOtherParam.mParam = &theResModel ;
		myOtherParam.mValue = &theValue ;

	gsl_multimin_function_fdf myFunct ;
		myFunct.df = GslGradLLHFunction ;
		myFunct.f = GslLLHFunction ;
		myFunct.fdf = GslLLHAndGradLLHFunction ;
		myFunct.n = myNParam ;
		myFunct.params = &myOtherParam ;

		myMultiMin.SetFunction(&myFunct) ;

	cDVector myX1(myNParam) ;

	sGSLMultiMinResult myGSLRes ;
 		myMultiMin.GSLOptim(myX1, theResStruct, theStopValue, theMaxIter, theVerbose) ;

		theResModel.VectorToRegArchParam(myX1) ;
	
	}

	void RegArchEstim(cRegArchModel& theModel, cRegArchValue& theValue,  sGSLMultiMinResult& theResStruct, cRegArchModel& theResModel, cDVector* theInitPoint, sGSLMultiMinAlgoParam* theAlgoParam)
	{
	uint myNParam = theModel.GetNParam() ;

		if (theInitPoint == NULL)
		{
			theResModel.SetDefaultInitPoint(theValue);
		}
		else
		{	if (theInitPoint->GetGSLVector() == NULL)
				theResModel.SetDefaultInitPoint(theValue) ;
			else
			{
				theResModel.VectorToRegArchParam(*theInitPoint);
			}
		}

	cDVector myX0(myNParam) ;
		theResModel.RegArchParamToVector(myX0) ;

	sGSLMultiMinAlgoParam* myAlgoParam ;
	bool myExist = (theAlgoParam != NULL) ;	
		if (!myExist)
		{
			myAlgoParam = (sGSLMultiMinAlgoParam *)malloc(sizeof(sGSLMultiMinAlgoParam)) ;
			myAlgoParam->mAlgoType = eBFGSTwo ;
			myAlgoParam->mNMaxIter = 200 ;
			myAlgoParam->mStepSize = 0.01 ;
			myAlgoParam->mTol = 0.01 ;
			myAlgoParam->mVerbose = false ;
		}
		else
		{
			myAlgoParam = theAlgoParam ;
		}

	cGSLMultiMin myMultiMin(myX0, myAlgoParam->mAlgoType) ;

	sParamOptimStruct myOtherParam ;
		myOtherParam.mParam = &theResModel ;
		myOtherParam.mValue = &theValue ;

	gsl_multimin_function_fdf myFunct ;
		myFunct.df = GslGradLLHFunction ;
		myFunct.f = GslLLHFunction ;
		myFunct.fdf = GslLLHAndGradLLHFunction ;
		myFunct.n = myNParam ;
		myFunct.params = &myOtherParam ;

		myMultiMin.SetFunction(&myFunct) ;

	cDVector myX1(myX0) ;

		myMultiMin.GSLOptim(myX1, theResStruct, *myAlgoParam) ;

		if (!myExist)
			free(myAlgoParam) ;

		theResModel.VectorToRegArchParam(myX1) ;
	}

	void RegArchEstim(cRegArchModel& theModel, cRegArchValue& theValue, cNLOPTResult& theResStruct, cRegArchModel& theResModel, cDVector* theInitPoint, cNLOPTAlgoParam* theAlgoParam)
	{
	uint myNParam = theModel.GetNParam();

		if (theModel.mMean != NULL)
			theResModel.SetMean(*(theModel.mMean));
		else
			theResModel.mMean = NULL;

		theResModel.SetVar(*(theModel.mVar));
		theResModel.SetResid(*(theModel.mResids));

		if (theInitPoint == NULL)
		{
			theResModel.SetDefaultInitPoint(theValue);
		}
		else
		{
			if (theInitPoint->GetGSLVector() == NULL)
				theResModel.SetDefaultInitPoint(theValue);
			else
			theResModel.VectorToRegArchParam(*theInitPoint);
		}
	cDVector myX0(myNParam);
		theResModel.RegArchParamToVector(myX0);

	cNLOPTAlgoParam myAlgoParam = cNLOPTAlgoParam();
		myAlgoParam.mfTol = theAlgoParam->mfTol;
		myAlgoParam.mMaxComputeTime = theAlgoParam->mMaxComputeTime;
		myAlgoParam.mMaxFuncEval = theAlgoParam->mMaxFuncEval;
		myAlgoParam.mMinimisation = theAlgoParam->mMinimisation;
		myAlgoParam.mStopVal = theAlgoParam->mStopVal;
		myAlgoParam.mxTol = theAlgoParam->mxTol;
		myAlgoParam.mVerbose = theAlgoParam->mVerbose;

	cNloptWrapperCpp myOptim = cNloptWrapperCpp();

	sParamOptimStruct myOtherParam;
		myOtherParam.mParam = &theResModel;
		myOtherParam.mValue = &theValue;
		myOtherParam.mVerbose = myAlgoParam.mVerbose;
		myOtherParam.mNFuncEval = 0;

		myOptim.SetAlgorithm(myAlgoParam.mAlgo, myNParam);
		myOptim.SetObjectiveFunc((nlopt_func)NloptLLHAndGradLLHFunction, false, &myOtherParam);

	double* myX1 = GSLVectorToDoubleStar(myX0);

		myOptim.Optimize(myX1, myNParam, myAlgoParam, theResStruct);

		cDVector myXOptim(myNParam, theResStruct.mXOptim);
			theResModel.VectorToRegArchParam(myXOptim);
	}

	void RegArchEstim(cRegArchModel& theModel, cRegArchValue& theValue, cNLOPTResult& theResStruct, cRegArchModel& theResModel, cDVector* theInitPoint, nlopt_algorithm theAlgo, double theStopValue, double thefTol, double thexTol, double theNMaxSec, int theMaxFuncEval, bool theMinimisation, bool theVerbose)
	{
	cNLOPTAlgoParam myAlgoParam =cNLOPTAlgoParam();
		myAlgoParam.mAlgo = theAlgo;
		myAlgoParam.mfTol = thefTol;
		myAlgoParam.mMaxComputeTime = theNMaxSec;
		myAlgoParam.mMaxFuncEval = theMaxFuncEval;
		myAlgoParam.mMinimisation = theMinimisation;
		myAlgoParam.mStopVal = theStopValue;
		myAlgoParam.mxTol = thexTol;
		myAlgoParam.mVerbose = theVerbose;
	
		RegArchEstim(theModel, theValue, theResStruct, theResModel, theInitPoint, &myAlgoParam);
	}

	void RegArchStatTable(cRegArchModel &theModel, cRegArchValue& theValue, cDMatrix& theTable)
	{
	uint myNParam = theModel.GetNParam();
	theTable.ReAlloc(myNParam, 4);
	cDVector myTeta(myNParam);
		theModel.RegArchParamToVector(myTeta);
		for (uint i = 0; i < myNParam; i++)
			theTable[i][0] = myTeta[i];
	cDMatrix myCov(myNParam, myNParam);
	int myError = 0;
		RegArchComputeCov(theModel, theValue, myCov, myError);

		if (myError == 0)
		{
			for (uint i = 0; i < myNParam; i++)
			{
				theTable[i][1] = sqrt(myCov[i][i]);
				theTable[i][2] = theTable[i][0] / theTable[i][1];
				theTable[i][3] = (1.0 - gsl_cdf_gaussian_P(fabs(theTable[i][2]), 1.0))*2.0;
			}
		}
		else
		{
			for (uint i = 0; i < myNParam; i++)
			{
				theTable[i][1] = -1;
				theTable[i][2] = -1;
				theTable[i][3] = -1;
			}

		}
	}
}
