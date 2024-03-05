#include "StdAfxRegArchLib.h"

namespace RegArchLib {

	/*!
	 * \fn cRegArchModel::cRegArchModel()
	 * \param None
	 * \details A simple constructor
	 */
	cRegArchModel::cRegArchModel()
	{
		mMean = NULL;
		mVar = NULL ;
		mResids = NULL ;
	}

	/*!
	 * \fn cRegArchModel::cRegArchModel(cCondMean* theMean, cAbstCondVar* theVar, cAbstResiduals* theResiduals)
	 * \param cCondMean* theMean: pointer to conditional mean model
	 * \param cAbstCondVar* theVar: pointer to conditional variance model
	 * \param cAbstResiduals* theResiduals: pointer to conditional residuals distribution model
	 * \details A simple constructor
	 */
	cRegArchModel::cRegArchModel(cCondMean& theMean, cAbstCondVar& theVar, cAbstResiduals& theResiduals)
	{
		if (mMean != NULL)
			delete mMean;
		mMean = new cCondMean(theMean);
		if (mVar != NULL)
			delete mVar;
		mVar = CreateOneRealCondVar(theVar);
		if (mResids != NULL)
			delete mResids;
		mResids = CreateRealCondResiduals(theResiduals);
		
 		MESS_CREAT("cRegArchModel") ;
	}

	/*!
	 * \fn cRegArchModel::cRegArchModel(cRegArchModel& theModel)
	 * \param cRegArchModel& theModel: the source
	 * \details recopy constructor
	 */
	cRegArchModel::cRegArchModel(const cRegArchModel& theModel)
	{	
		if (theModel.mMean != NULL)
			mMean = new cCondMean(*(theModel.mMean)) ;
		else
			mMean = new cCondMean() ;
		mVar = CreateOneRealCondVar(*(theModel.mVar));
		mResids = CreateRealCondResiduals(*(theModel.mResids));

		MESS_CREAT("cRegArchModel") ;
	}

	/*!
	 * \fn cRegArchModel::~cRegArchModel()
	 */
	cRegArchModel::~cRegArchModel()
	{	if (mMean != NULL)
		{	//mMean->Delete() ;
			delete mMean;
		}
		if (mVar != NULL)
		{	//mVar->Delete() ;
			delete mVar;
		}
		if (mResids != NULL)
		{	//mResids->Delete() ;
			delete mResids;
		}
		MESS_DESTR("cRegArchModel") ;
	}

	/*!
	 * \fn cRegArchModel& cRegArchModel::operator=(cRegArchModel& theRegArchModel)
	 * \param cRegArchModel& theRegArchModel: the source
	 * \details = operator for cRegArchModel
	 */
	cRegArchModel& cRegArchModel::operator=(cRegArchModel& theRegArchModel)
	{
		
		*mMean =  *(theRegArchModel.mMean) ;
		//mVar = theRegArchModel.mVar->PtrCopy() ;
		*mVar = *(theRegArchModel.mVar);
		//mResids = theRegArchModel.mResids->PtrCopy() ;
		*mResids = *(theRegArchModel.mResids);
		return *this ;
	}

	void cRegArchModel::Delete(void)
	{
		if (mMean != NULL)
		{	//mMean->Delete() ;
			delete mMean;
			mMean = NULL;
		}
		if (mVar != NULL)
		{	//mVar->Delete() ;
			delete mVar;
			mVar = NULL;
		}
		if (mResids != NULL)
		{	//mResids->Delete() ;
			delete mResids;
			mResids = NULL;
		}
	}

	/*
	void cRegArchModel::ReAllocMean(uint theNewSize)
	{
		if (mMean != NULL)
		{
			mMean->Delete() ;
			delete mMean ;
		}
		mMean = new cCondMean(theNewSize) ;
	}
	*/

	/*!
	 * \fn void cRegArchModel::SetMean(cCondMean& theCondMean)
	 * \param cCondMean& theCondMean: the conditional mean model
	 */
	void cRegArchModel::SetMean(cCondMean& theCondMean)
	{
		if (mMean != NULL)
		{
			delete mMean ;
		}

		mMean = new cCondMean(theCondMean) ;
	}

	/*!
	 * \fn void cRegArchModel::GetNMean(void)
	 * \param void
	 * \details return the number of mean components
	 */
	int cRegArchModel::GetNMean(void)
	{
		if (mMean == NULL)
			return 0 ;
		else
			return mMean->GetNMean() ;
	}

	/*!
	 * \fn void cRegArchModel::AddOneMean(cAbstCondMean& theOneMean)
	 * \param cAbstCondMean& theOneMean: the conditional mean component model
	 * \details Add a new mean component
	 */
	void cRegArchModel::AddOneMean(cAbstCondMean& theOneMean)
	{
		if (mMean == NULL)
			mMean = new cCondMean();
		mMean->AddOneMean(theOneMean) ;
	}

	/*!
	 * \fn void cRegArchModel::GetOneMean(int theNumMean)
	 * \param int theNumMean: the index of the mean component
	 * \details Return theNumMean th mean component
	 */
	cAbstCondMean* cRegArchModel::GetOneMean(int theNumMean)
	{	if (mMean != NULL)
			return mMean->GetOneMean(theNumMean) ;
		else
			return NULL ;
	}

	bool cRegArchModel::IsGoodMeanType(eCondMeanEnum theMeanEnum, int theIndex)
	{
		if (mMean != NULL)
		{	cAbstCondMean* myCondMean = mMean->GetOneMean(theIndex) ;
			return (myCondMean->GetCondMeanType() == theMeanEnum) ;
		}
		else
			return false ;
	}

#ifdef _RDLL_
	void cRegArchModel::PrintMean(void)
	{
		if (mMean != NULL)
			mMean->Print();
	}

	void cRegArchModel::PrintVar(void)
	{
		if (mVar != NULL)
			mVar->Print();
	}

	void cRegArchModel::PrintResiduals(void)
	{
		if (mResids != NULL)
			mResids->Print();
	}

	void cRegArchModel::Print(void)
	{
		Rprintf("Regression with ARCH type residuals parameters:\n");
		Rprintf("-----------------------------------------------\n");
		PrintResiduals();
		Rprintf("\n");
		PrintMean();
		Rprintf("\n");
		Rprintf("\nConditional variance parameters:\n");
		Rprintf("--------------------------------\n");
		PrintVar();
		Rprintf("\n");
	}

#else

	/*!
	* \fn void cRegArchModel::PrintMean(ostream& theOut) const
	* \param ostream& theOut: output stream (screen or file). Default: cout
	*/
	void cRegArchModel::PrintMean(ostream& theOut) const
    {	if (mMean != NULL)
    		mMean->Print(theOut) ;
    }

	/*!
	* \fn void cRegArchModel::PrintVar(ostream& theOut) const
	* \param ostream& theOut: output stream (screen or file). Default: cout
	*/
	void cRegArchModel::PrintVar(ostream& theOut) const
	{
		if (mVar != NULL)
			mVar->Print(theOut);
	}
	/*!
	* \fn void cRegArchModel::PrintResiduals(ostream& theOut) const
	* \param ostream& theOut: output stream (screen or file). Default: cout
	*/
	void cRegArchModel::PrintResiduals(ostream& theOut) const
	{
		if (mResids != NULL)
			mResids->Print(theOut);
	}

	/*!
	* \fn void cRegArchModel::Print(ostream& theOut) const
	* \param ostream& theOut: output stream (screen or file). Default: cout
	*/
	void cRegArchModel::Print(ostream& theOut) const
	{
		theOut << "Regression with ARCH type residuals parameters:" << endl;
		theOut << "-----------------------------------------------" << endl;
		PrintResiduals(theOut);
		theOut << endl;
		PrintMean(theOut);
		theOut << endl;
		theOut << "Conditional variance parameters:" << endl;
		theOut << "--------------------------------" << endl;
		PrintVar(theOut);
		theOut << endl;

	}

#endif //_RDLL_
	/*!
	 * \fn void cRegArchModel::SetVar(cAbstCondVar& theCondVar)
	 * \param cCondVar& theCondVar: the conditional variance model
	 */
	void cRegArchModel::SetVar(cAbstCondVar& theCondVar)
	{
		if (mVar != NULL)
		{
			delete mVar ;
			mVar = NULL;
		}
		mVar = CreateOneRealCondVar(theCondVar);

	}

	cAbstCondVar* cRegArchModel::GetVar(void)
	{
		return mVar  ;
	}

	/*!
	 * \fn void cRegArchModel::SetResid(cAbstResiduals& theCondResiduals)
	 * \param cAbstResiduals& theCondResiduals: the conditional residuals model
	 */
	void cRegArchModel::SetResid(cAbstResiduals& theCondResiduals)
	{
		if (mResids != NULL)
		{
			delete mResids ;
		}
	//	mResids = theCondResiduals.PtrCopy() ;
		mResids = CreateRealCondResiduals(theCondResiduals) ;	
	}

	cAbstResiduals*  cRegArchModel::GetResid(void)
	{
			return mResids ;
	}

	void cRegArchModel::SetDefaultInitPoint(cRegArchValue& theValue)
	{
		if (mMean != NULL)
		{
			mMean->SetDefaultInitPoint(theValue);
			for (uint t = 0; t < theValue.mYt.GetSize(); t++)
			{
				theValue.mMt[t] = mMean->ComputeMean(t, theValue);
				theValue.mUt[t] = theValue.mYt[t] - theValue.mMt[t];
			}
		}
		else
		{
			for (uint t = 0; t < theValue.mYt.GetSize(); t++)
			{
				theValue.mUt[t] = theValue.mYt[t];
			}
		}
		mVar->SetDefaultInitPoint(theValue) ;
		mResids->SetDefaultInitPoint() ;
	}

	uint cRegArchModel::GetNParam(void) const
    {
		uint myNParam = 0;
		if (mMean != NULL)
			myNParam = mMean->GetNParam();
		if (mVar != NULL)
			myNParam += mVar->GetNParam();
		if (mResids != NULL)
			myNParam += mResids->GetNParam();

    	return myNParam;
    }

	uint cRegArchModel::GetNLags(void) const
    {
    	if (mMean == NULL)
    		return mVar->GetNLags() ;
    	else
    		return MAX(mMean->GetNLags(), mVar->GetNLags()) ;
    }

	void cRegArchModel::ComputeGrad(uint theDate, cRegArchValue& theValue, cRegArchGradient& theGradData)
    {
		FillValue(theDate, *this, theValue);
		mVar->ComputeGrad(theDate, theValue, theGradData, mResids);
	double mySigmat = sqrt(theValue.mHt[theDate]);
		theGradData.mCurrentGradSigma = theGradData.mCurrentGradVar / (2.0 * mySigmat);
		if (mMean != NULL)
		{
			theGradData.mCurrentGradMu = 0.0;
			mMean->ComputeGrad(theDate, theValue, theGradData);
		}
		mResids->ComputeGrad(theDate, theValue, theGradData);
		theGradData.mCurrentGradEps = -1.0*(theValue.mEpst[theDate] * theGradData.mCurrentGradSigma + theGradData.mCurrentGradMu) / mySigmat;
    }

	void cRegArchModel::NumericComputeGrad(uint theDate, cRegArchValue& theData, cRegArchGradient& theGradData, cAbstResiduals* theResiduals, cNumericDerivative& theNumDeriv)
	{
		FillValueForNumericGrad(theDate, *this, theData, theNumDeriv);
	uint myNParam = GetNParam();
	uint myNDistrParam = mResids->GetNParam();
	uint myNMeanVarParam = myNParam - myNDistrParam;
	cDVector myParam(myNParam);
		RegArchParamToVector(myParam);
	double myF0Mean = 0;
		if (mMean != NULL)
			myF0Mean = theData.mMt[theDate];
	double myF0Var = theData.mHt[theDate];
	double myF0Dens = mResids->LogDensity(theData.mEpst[theDate]);
	double myF0Eps = theData.mEpst[theDate];
		for (uint i = 0; i < myNParam; i++)
		{
		double myh1 = theNumDeriv.mh1[i];
		double myF1Mean = 0;
			if (mMean != NULL)
				myF1Mean = theNumDeriv.mValueForGrad[i].mMt[theDate];
		double myF1Var = theNumDeriv.mValueForGrad[i].mHt[theDate];
		double myF1Eps = theNumDeriv.mValueForGrad[i].mEpst[theDate];
		double myF1Dens;
			if (i >= myNMeanVarParam)
				myF1Dens = theNumDeriv.mLogDensForGrad[i - myNMeanVarParam];
			if (mMean != NULL)
					theGradData.mCurrentGradMu[i] = (myF1Mean - myF0Mean) / myh1;
				theGradData.mCurrentGradVar[i] = (myF1Var - myF0Var) / myh1;
			theGradData.mCurrentGradEps[i] = (myF1Eps - myF0Eps) / myh1;
			if (i >= myNMeanVarParam)
				theGradData.mCurrentGradLogDens[i] = (myF1Dens - myF0Dens) / myh1;
		}
	double mySigmat = sqrt(theData.mHt[theDate]);
		theGradData.mCurrentGradSigma = theGradData.mCurrentGradVar / (2.0 * mySigmat);
	double myh0 = theNumDeriv.Geth();
	double myh1 = fabs(myh0*theData.mEpst[theDate]);
		if (myh1 < 1e-16)
			myh1 = myh0;
	double myF1 = mResids->LogDensity(theData.mEpst[theDate]+myh1);
		theGradData.mCurrentDiffLogDensity = (myF1 - mResids->LogDensity(myF0Eps)) / myh1;	
	}

	void cRegArchModel::ComputeGradAndHess(uint theDate, cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
		ComputeGrad(theDate, theData, theGradData);
		mVar->ComputeHess(theDate, theData, theGradData, theHessData, mResids);
		if (mMean != NULL)
			mMean->ComputeHess(theDate, theData, theGradData, theHessData);
		mResids->ComputeHess(theDate, theData, theGradData, theHessData);
		theHessData.ComputeHessSigmatFromHessVart(theDate, theData, theGradData);
	double mySigmat = sqrt(theData.mHt[theDate]);
		/* Hessien Eps(t, theta) */
		theHessData.mCurrentHessEps = theGradData.mCurrentGradSigma * Transpose(theGradData.mCurrentGradMu) / theData.mHt[theDate];
		theHessData.mCurrentHessEps += theGradData.mCurrentGradMu * Transpose(theGradData.mCurrentGradSigma) / theData.mHt[theDate];
		theHessData.mCurrentHessEps -= theHessData.mCurrentHessMu / mySigmat;
		theHessData.mCurrentHessEps += 2 * theData.mUt[theDate] / pow(mySigmat, 3)*theGradData.mCurrentGradSigma*Transpose(theGradData.mCurrentGradSigma);
		theHessData.mCurrentHessEps -= theData.mUt[theDate] / theData.mHt[theDate] * theHessData.mCurrentHessSigma;
	}

	void cRegArchModel::NumericComputeGradAndHess(uint theDate, cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals, cNumericDerivative& theNumDeriv)
	{
		FillValueForNumericGradAndHess(theDate, *this, theData, theNumDeriv);
	uint myNParam = GetNParam();
	uint myNDistrParam = mResids->GetNParam();
	uint myNMeanVarParam = myNParam - myNDistrParam;
	cDVector myParam(myNParam);
		RegArchParamToVector(myParam);
	double myF0Mean = theData.mMt[theDate];
	double myF0Var = theData.mHt[theDate];
	double myF0Dens = mResids->LogDensity(theData.mEpst[theDate]);
	double myF0Eps = theData.mEpst[theDate];
	double* myF1Mean = new double[myNParam];
	double* myF1Var = new double[myNParam];
	double* myF1Dens = new double[myNParam];
	double* myF1Eps = new double[myNParam];
		for (uint i = 0; i < myNParam; i++)
		{
			myF1Mean[i] = theNumDeriv.mValueForGrad[i].mMt[theDate];
			myF1Var[i] = theNumDeriv.mValueForGrad[i].mHt[theDate];
			myF1Eps[i] = theNumDeriv.mValueForGrad[i].mEpst[theDate];
			if (i >= myNMeanVarParam)
				myF1Dens[i] = theNumDeriv.mLogDensForGrad[i - myNMeanVarParam];
			theGradData.mCurrentGradMu[i] = (myF1Mean[i] - myF0Mean) / theNumDeriv.mh1[i];
			theGradData.mCurrentGradVar[i] = (myF1Var[i] - myF0Var) / theNumDeriv.mh1[i];
			theGradData.mCurrentGradEps[i] = (myF1Eps[i] - myF0Eps) / theNumDeriv.mh1[i];
			if (i >= myNMeanVarParam)
				theGradData.mCurrentGradLogDens[i] = (myF1Dens[i] - myF0Dens) / theNumDeriv.mh1[i];
		}
		for (uint i = 0; i < myNParam; i++)
		{
			for (uint j = i; j < myNParam; j++)
			{
			double myF2 = theNumDeriv.mValueForHess[i][j].mMt[theDate];
			double myh12 = theNumDeriv.mh1[i] * theNumDeriv.mh2[i][j];
				theHessData.mCurrentHessMu[i][j] = theHessData.mCurrentHessMu[j][i] = (myF2 - myF1Mean[i] - myF1Mean[j] + myF0Mean) / myh12;
				myF2 = theNumDeriv.mValueForHess[i][j].mHt[theDate];
				theHessData.mCurrentHessVar[i][j] = theHessData.mCurrentHessVar[j][i] = (myF2 - myF1Var[i] - myF1Var[j] + myF0Var) / myh12;
				myF2 = theNumDeriv.mValueForHess[i][j].mEpst[theDate];
				theHessData.mCurrentHessEps[i][j] = theHessData.mCurrentHessEps[j][i] = (myF2 - myF1Eps[i] - myF1Eps[j] + myF0Eps) / myh12;
				if (i >= myNMeanVarParam)
				{
					myF2 = theNumDeriv.mLogDensForHess[i-myNMeanVarParam][j];
					theHessData.mCurrentHessDens[i][j] = theHessData.mCurrentHessDens[j][i] = (myF2 - myF1Dens[i] - myF1Dens[j] + myF0Dens) / myh12;
				}
			}
		}
	double mySigmat = sqrt(theData.mHt[theDate]);
		theGradData.mCurrentGradSigma = theGradData.mCurrentGradVar / (2.0 * mySigmat);
	double myh0 = theNumDeriv.Geth();
	double myh1 = fabs(myh0*theData.mEpst[theDate]);
		if (myh1 < 1e-16)
			myh1 = myh0;
	double myF1 = mResids->LogDensity(theData.mEpst[theDate] + myh1);
		theGradData.mCurrentDiffLogDensity = (myF1 - myF0Dens) / myh1;
		for (uint i = 0; i < myNDistrParam; i++)
		{
		uint k = myNMeanVarParam + i;
			theHessData.mCurrentGradDiffLogDensity[k] = (theNumDeriv.mGradDiffForGrad[i] - myF1 - theNumDeriv.mLogDensForGrad[i] + myF0Dens) / (myh1*theNumDeriv.mh1[k]);
		}
		theHessData.ComputeHessSigmatFromHessVart(theDate, theData, theGradData);
		delete[] myF1Mean;
		delete[] myF1Var;
		delete[] myF1Dens;
		delete[] myF1Eps;
	}

	void cRegArchModel::RegArchParamToVector(cDVector& theDestVect) const
    {
    uint myIndex = 0 ;
    	if (mMean != NULL)
    	{	mMean->RegArchParamToVector(theDestVect, myIndex) ;
    		myIndex += mMean->GetNParam() ;
    	}
    	mVar->RegArchParamToVector(theDestVect, myIndex) ;
     	myIndex += mVar->GetNParam() ;
    	mResids->RegArchParamToVector(theDestVect, myIndex) ;
    }

	void cRegArchModel::VectorToRegArchParam(const cDVector& theSrcParam)
    {
    uint myIndex = 0 ;
    	if (mMean != NULL)
    	{	mMean->VectorToRegArchParam(theSrcParam, myIndex) ;
    		myIndex += mMean->GetNParam() ;
    	}
    	mVar->VectorToRegArchParam(theSrcParam, myIndex) ;
    	myIndex += mVar->GetNParam() ;
    	mResids->VectorToRegArchParam(theSrcParam, myIndex) ;
    
    	if ((mVar->GetCondVarType() == eEgarch) && (mResids->GetNParam() > 0))
    	{	
    	cEgarch* myVar = dynamic_cast<cEgarch *>(mVar) ;
    		if (myVar != NULL)
    		{	
    		double myAux = mResids->ComputeEspAbsEps() ;
    			myVar->SetEspAbsEps(myAux) ;
    		}
    	}
    
    }

	void cRegArchModel::GetParamName(char** theName)
	{
	uint myNParam = GetNParam();
		if (mMean != NULL)
			mMean->GetParamName(theName);
	uint myIndex = mMean->GetNParam();
		mVar->GetParamName(myIndex, theName);
		myIndex += mVar->GetNParam();
		mResids->GetParamName(myIndex, theName);	
	}

	void cRegArchModel::GetParamName(string theName[])
	{
	uint myNParam = GetNParam();
		if (theName == NULL)
			theName = new string[myNParam];
		if (mMean != NULL)
			mMean->GetParamName(theName);
	uint myIndex = mMean->GetNParam();
		mVar->GetParamName(myIndex, theName);
		myIndex += mVar->GetNParam();
		mResids->GetParamName(myIndex, theName);
	}
	
	void cRegArchModel::GetParamName(cTabOfString& theName)
	{
	uint myNParam = GetNParam();
		theName.ReAlloc(myNParam);
		uint myIndex = 0;
		if (mMean != NULL)
		{
			mMean->GetParamName(theName.mTab);
			myIndex = mMean->GetNParam();
		}
		mVar->GetParamName(myIndex, theName.mTab);
		myIndex += mVar->GetNParam();
		mResids->GetParamName(myIndex, theName.mTab);
	}
	
}//namespace

