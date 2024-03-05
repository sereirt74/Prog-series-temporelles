#include "StdAfxRegArchLib.h"

/*!
 \file cAbstCondMean.cpp
 \brief sources for abstract class cConst methods.

 \author Jean-Baptiste Durand, Ollivier TARAMASCO
 \date dec-18-2006 - Last change feb-18-2011
*/
namespace RegArchLib {
    /*!
     * \fn cRegArchGradient::cRegArchGradient(uint theNPast, uint theNMeanParam, uint theNVarParam, uint theNDistrParameter)
     *
     * \param uint theNPast: maximum number of past lags used in gradient computations
     * \param uint theNMeanParam: number of conditional mean parameters
     * \param uint theNVarParam:  number of conditional var parameters
     * \param uint theNDistrParameter:  number of conditional distribution parameters
     *
     */
    cRegArchGradient::cRegArchGradient(uint theNPast, uint theNMeanParam, uint theNVarParam, uint theNDistrParameter)
    {
    	mvNPast = theNPast ;
    	mGradHt = new cDVector[theNPast] ;
    	mGradMt = new cDVector[theNPast] ;
    	mGradEpst = new cDVector[theNPast] ;
    	mvNParam = theNMeanParam + theNVarParam + theNDistrParameter ;
    	mvNMeanParam = theNMeanParam ;
    	mvNVarParam = theNVarParam ;
    	mvNDistrParameter = theNDistrParameter ;
    	for (uint i = 0 ; i < theNPast ; i++)
    	{	mGradHt[i].ReAlloc(mvNParam) ;
    		mGradMt[i].ReAlloc(mvNParam) ;
    		mGradEpst[i].ReAlloc(mvNParam) ;
    	}
    	mCurrentGradVar.ReAlloc(mvNParam) ;
    	mCurrentGradSigma.ReAlloc(mvNParam) ;
    	mCurrentGradMu.ReAlloc(mvNParam) ;
		mCurrentGradLogDens.ReAlloc(mvNParam);
    	mCurrentGradEps.ReAlloc(mvNParam) ;
		mCurrentDiffLogDensity = 0.0;
    }
    
    cRegArchGradient::cRegArchGradient(cRegArchModel* theModel)
    {
    	mvNPast = mvNParam = 0 ;
	    if (theModel != NULL)
    	{	mvNPast = theModel->GetNLags() ;
			if (theModel->mMean != NULL)
			{
				mvNMeanParam = theModel->mMean->GetNParam();
			}
			else
    			mvNMeanParam = 0 ;
    		mvNVarParam = theModel->mVar->GetNParam() ;
    		mvNDistrParameter = theModel->mResids->GetNParam() ;
    		mvNParam = mvNMeanParam + mvNVarParam + mvNDistrParameter ;
    	}
    	mGradHt = new cDVector[mvNPast] ;
    	mGradMt = new cDVector[mvNPast] ;
    	mGradEpst = new cDVector[mvNPast] ;
    	for (uint i = 0 ; i < mvNPast ; i++)
    	{	mGradHt[i].ReAlloc(mvNParam) ;
    		mGradMt[i].ReAlloc(mvNParam) ;
    		mGradEpst[i].ReAlloc(mvNParam) ;
    	}
    	mCurrentGradVar.ReAlloc(mvNParam) ;
    	mCurrentGradSigma.ReAlloc(mvNParam) ;
    	mCurrentGradMu.ReAlloc(mvNParam) ;
		mCurrentGradLogDens.ReAlloc(mvNParam);
    	mCurrentGradEps.ReAlloc(mvNParam) ;
		mCurrentDiffLogDensity = 0.0;
	}
    
    cRegArchGradient::~cRegArchGradient()
    {
    	for (uint i = 0 ; i < mvNPast ; i++)
    	{	mGradHt[i].Delete() ;
    		mGradMt[i].Delete() ;
    		mGradEpst[i].Delete() ;
    	}
        if (mvNPast > 0)
        {
            delete[] mGradHt;
            delete[] mGradMt;
            delete[] mGradEpst;
        }
        mvNPast = 0;
    	mCurrentGradVar.Delete() ;
    	mCurrentGradSigma.Delete() ;
    	mCurrentGradMu.Delete() ;
    	mCurrentGradLogDens.Delete() ;
    	mCurrentGradEps.Delete() ;
    	mvNMeanParam =  mvNVarParam = mvNDistrParameter = mvNParam = 0 ;
    }
    
    void cRegArchGradient::Delete(void)
    {
    	for (uint i = 0 ; i < mvNPast ; i++)
    	{	mGradHt[i].Delete() ;
    		mGradMt[i].Delete() ;
    		mGradEpst[i].Delete() ;
    	}
    	delete [] mGradHt ;
    	delete [] mGradMt ;
    	delete [] mGradEpst ;
    	mvNPast = 0 ;
    	mCurrentGradVar.Delete() ;
    	mCurrentGradSigma.Delete() ;
    	mCurrentGradMu.Delete() ;
    	mCurrentGradLogDens.Delete() ;
    	mCurrentGradEps.Delete() ;
    	mvNMeanParam =  mvNVarParam = mvNDistrParameter = mvNParam = 0 ;
    }
    
    void cRegArchGradient::ReAlloc(uint theNPast, uint theNMeanParam, uint theNVarParam, uint theNDistrParameter) 
    {
    	Delete() ;
    	mvNPast = theNPast ;
    	mGradHt = new cDVector[theNPast] ;
    	mGradMt = new cDVector[theNPast] ;
    	mGradEpst = new cDVector[theNPast] ;
    	mvNParam = theNMeanParam + theNVarParam + theNDistrParameter ;
    	mvNMeanParam = theNMeanParam ;
    	mvNVarParam = theNVarParam ;
    	mvNDistrParameter = theNDistrParameter ;
    	for (uint i = 0 ; i < theNPast ; i++)
    	{	mGradHt[i].ReAlloc(mvNParam) ;
    		mGradMt[i].ReAlloc(mvNParam) ;
    		mGradEpst[i].ReAlloc(mvNParam) ;
    	}
    	mCurrentGradVar.ReAlloc(mvNParam) ;
    	mCurrentGradSigma.ReAlloc(mvNParam) ;
    	mCurrentGradMu.ReAlloc(mvNParam) ;
		mCurrentGradLogDens.ReAlloc(mvNParam);
    	mCurrentGradEps.ReAlloc(mvNParam) ;
    }
    
    void cRegArchGradient::ReAlloc(cRegArchModel* theModel)
    {
    	mvNPast = mvNParam = 0 ;
    	if (theModel != NULL)
    	{	mvNPast = theModel->GetNLags() ;
			if (theModel->mMean == NULL)
				mvNMeanParam = 0;
			else
				mvNMeanParam = theModel->mMean->GetNParam();
    		mvNVarParam = theModel->mVar->GetNParam() ;
    		mvNDistrParameter = theModel->mResids->GetNParam() ;
    		mvNParam = mvNMeanParam + mvNVarParam + mvNDistrParameter ;
    	}
    	mGradHt = new cDVector[mvNPast] ;
    	mGradMt = new cDVector[mvNPast] ;
    	mGradEpst = new cDVector[mvNPast] ;
    	for (uint i = 0 ; i < mvNPast ; i++)
    	{	mGradHt[i].ReAlloc(mvNParam) ;
    		mGradMt[i].ReAlloc(mvNParam) ;
    		mGradEpst[i].ReAlloc(mvNParam) ;
    	}
    	mCurrentGradVar.ReAlloc(mvNParam) ;
    	mCurrentGradSigma.ReAlloc(mvNParam) ;
    	mCurrentGradMu.ReAlloc(mvNParam) ;
		mCurrentGradLogDens.ReAlloc(mvNParam);
    	mCurrentGradEps.ReAlloc(mvNParam) ;
    }
    
	void cRegArchGradient::ReInitialize(void)
	{
		for (uint i = 0; i < mvNPast; i++)
		{
			mGradHt[i] = 0.0 ;
			mGradMt[i] = 0.0 ;
			mGradEpst[i] = 0.0 ;
		}
		mCurrentGradVar = 0.0 ;
		mCurrentGradSigma = 0.0 ;
		mCurrentGradMu = 0.0 ;
		mCurrentGradLogDens = 0.0 ;
		mCurrentGradEps = 0.0 ;
		mCurrentDiffLogDensity = 0.0;
	}

    uint cRegArchGradient::GetNPast(void) const
    {
    	return mvNPast ;
    }
    
    uint cRegArchGradient::GetNParam(void)
    {
    	return mvNParam ;
    }
    
    uint cRegArchGradient::GetNMeanParam(void)
    {
    	return mvNMeanParam ;
    }
    
    uint cRegArchGradient::GetNVarParam(void)
    {
    	return mvNVarParam ;
    }
    
    uint cRegArchGradient::GetNDistrParameter(void)
    {
    	return mvNDistrParameter ;
    }
    
    void cRegArchGradient::Update(void)
    {
    	if (mvNPast > 0)
    	{	for (int i = (int)mvNPast - 1 ; i > 0 ; i--)
    		{	mGradHt[i] = mGradHt[i-1] ;
    			mGradMt[i] = mGradMt[i-1] ;
    			mGradEpst[i] = mGradEpst[i-1] ;
    		}
    		mGradMt[0] = mCurrentGradMu  ;
    		mGradHt[0] = mCurrentGradVar ;
    		mGradEpst[0] = mCurrentGradEps;
    	}
    }

	cRegArchGradient& cRegArchGradient::operator =(cRegArchGradient& theSrc)
	{

		mvNPast = theSrc.GetNPast(); 
		mvNParam = theSrc.GetNParam();
		mvNMeanParam = theSrc.GetNMeanParam();
		mvNVarParam = theSrc.GetNVarParam();
		mvNDistrParameter = theSrc.GetNDistrParameter();

		ReAlloc(mvNPast, mvNMeanParam, mvNVarParam, mvNDistrParameter);

		mCurrentGradVar = theSrc.mCurrentGradVar;
		mCurrentGradSigma = theSrc.mCurrentGradSigma;
		mCurrentGradMu = theSrc.mCurrentGradMu;
		mCurrentGradLogDens = theSrc.mCurrentGradLogDens;
		mCurrentGradEps = theSrc.mCurrentGradEps;
		mCurrentDiffLogDensity = theSrc.mCurrentDiffLogDensity;
		for (uint i = 0; i < mvNPast; i++)
		{
			mGradHt[i] = theSrc.mGradHt[i];
			mGradMt[i] = theSrc.mGradMt[i];
			mGradEpst[i] = theSrc.mGradEpst[i];
		}

		return *this;
	}

}//namespace
