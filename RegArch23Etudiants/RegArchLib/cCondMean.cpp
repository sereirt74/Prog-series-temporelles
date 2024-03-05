#include "StdAfxRegArchLib.h"
/*!
 \file cAbstCondVar.cpp
 \brief sources for abstract class cAbstCondVar methods.

 \author Jean-Baptiste DURAND, Ollivier TARAMASCO 
 \date dec-18-2006 - Last change feb-18-2011
*/

namespace RegArchLib {
	/*!
	 * \fn cCondMean::cCondMean(uint theNCondMean)
	 * \param uint theNCondMean: number of conditional means
	 * \details mvCondMean = theNCondMean
	 */
	cCondMean::cCondMean()
	{
		for (uint i = 0; i < MAX_COND_MEAN; i++)
			mvCondMean[i] = NULL;
		MESS_CREAT("cCondMean") ;
	}

	/*!
	 * \fn cCondMean::cCondMean(const cCondMean& theCondMean)
	 * \param cCondMean& theCondMean: conditional mean
	 * \details Recopy constructor
	 */
	cCondMean::cCondMean(const cCondMean& theCondMean)
	{
		for (uint i = 0; i < MAX_COND_MEAN; i++)
			mvCondMean[i] = NULL;
		for (uint i = 1; i < MAX_COND_MEAN; i++)
			if (theCondMean.mvCondMean[i] != NULL)
			{
				mvCondMean[i] = CreateOneRealCondMean(*(theCondMean.mvCondMean[i]));
			}
		MESS_CREAT("cCondMean") ;
	}

	/*!
	 * \fn cCondMean::~cCondMean()
	 * \param None
	 * \details simple destructor
	 */
	cCondMean::~cCondMean()
	{	Delete() ;
		MESS_DESTR("cCondMean") ;
	}

	/*!
	 * \fn void cCondMean::Delete(void)
	 * \param void
	 * \details free memory used par the cCondMean class
	 */
	void cCondMean::Delete(void)
	{	for (uint i = 1; i < MAX_COND_MEAN; i++)
			if (mvCondMean[i] != NULL)
			{
				delete mvCondMean[i];
				mvCondMean[i] = NULL;
			}
	}

	/*!
	 * \fn inline uint cCondMean::GetNMean(void) const
	 * \param void
	 * \brief return mvNCondMean
	 */
	uint cCondMean::GetNMean(void) const
	{
	uint myAux = 0;
		for (uint i = 1; i < MAX_COND_MEAN; i++)
			myAux += 1 * (mvCondMean[i] != NULL);
		return myAux;
	}

	/*!
	* \fn  void cCondMean::AddOneMean(cAbstCondMean* theAbstCondMean)
	* \param cAbstCondMean* theAbstCondMean: conditional mean component to be copied in the mCondMean array.
	* \brief *mvCondMean[theWhatMean] = *theAbstCondMean
	*/
	void cCondMean::AddOneMean(cAbstCondMean& theAbstCondMean)
	{
	int myNum = (int)theAbstCondMean.GetCondMeanType();
		if (mvCondMean[myNum] != NULL)
			mvCondMean[myNum]->Delete();
		mvCondMean[myNum] = CreateOneRealCondMean(theAbstCondMean);
	}

	/*!
	 * \fn cAbstCondMean** cCondMean::GetCondMean(void)
	 * \\details return mvCondMean
	 */
	cAbstCondMean** cCondMean::GetCondMean(void) const
	{	return (cAbstCondMean **)mvCondMean ;
	}

	/*!
	 * \fn cAbstCondMean* cCondMean::GetOneMean(uint theIndex) const
	 * \param uint theIndex: index of component to be returned
	 * \\details return mvCondMean[theWhateMean] 
	 */
	cAbstCondMean* cCondMean::GetOneMean(uint theIndex) const
	{
		return mvCondMean[theIndex] ;
	}

	/*!
	 * \fn void cCondMean::GetCondMeanType(eCondMeanEnum* theCodeType) const
	 * \param eCondMeanEnum* theCodeType: array of all conditional mean component codes
	 * \details fill theCodeType array
	 */
	void cCondMean::GetCondMeanType(eCondMeanEnum* theCodeType) const
	{	
	uint j = 0;
		for (uint i = 1 ; i < MAX_COND_MEAN ; i++)
			if (mvCondMean[i] != NULL)
				theCodeType[j++] = (eCondMeanEnum)i  ;
	}

	/*!
	 * \fn void cCondMean::Print(ostream& theOut) const
	 * \param ostream& theOut: output stream (file or screen). Default cout.
	 */
#ifndef _RDLL_
	void cCondMean::Print(ostream& theOut) const
	{	theOut << "Conditional mean parameters:" << endl ;
		theOut << "----------------------------" << endl ;
		for (uint i = 0 ; i < MAX_COND_MEAN ; i++)
			if (mvCondMean[i] != NULL)
				mvCondMean[i]->Print(theOut) ;
	}
#else
	void cCondMean::Print(void)
	{
		Rprintf("Conditional mean parameters:\n");
		Rprintf("----------------------------\n");
		for (uint i = 1; i < MAX_COND_MEAN; i++)
			if (mvCondMean[i] != NULL)
				mvCondMean[i]->Print();
	}
#endif // _RDLL_

	/*!
	 * \fn double cCondMean::Get(uint theNumMean, uint theIndex, uint theNumParam)
	 * \param uint theNumMean: index of conditional mean
	 * \param uint theNumParam: index of parameter
	 */
	double cCondMean::Get(uint theNumMean, uint theIndex, uint theNumParam)
	{
		if (theNumMean < MAX_COND_MEAN)
		{
			if (mvCondMean[theNumMean] != NULL)
				return mvCondMean[theNumMean]->Get(theIndex, theNumParam);
			else
				throw cError("cCondMean::Get bad index");
		}
		else
			throw cError("cCondMean::Get bad index");
	}
	
	void cCondMean::SetDefaultInitPoint(double theMean, double theVar)
	{
		for (uint i = 1; i < MAX_COND_MEAN; i++)
			if (mvCondMean[i] != NULL)
				mvCondMean[i]->SetDefaultInitPoint(theMean, theVar) ;
	}

	void cCondMean::SetDefaultInitPoint(cRegArchValue& theValue)
	{
		for (uint i = 1; i < MAX_COND_MEAN; i++)
			if (mvCondMean[i] != NULL)
				mvCondMean[i]->SetDefaultInitPoint(theValue);
	}

	/*!
	 * \fn ostream& operator <<(ostream& theOut, const cCondMean& theCondMean)
	 * \param ostream& theOut: output (file or screen).
	 * \param const cCondMean& theCondMean: the conditional mean class to be printed.
	 */
#ifndef _RDLL_
	ostream& operator <<(ostream& theOut, const cCondMean& theCondMean)
	{	theOut << "Conditional mean parameters:" << endl ;
		theOut << "----------------------------" << endl ;
		for (uint i = 0 ; i < MAX_COND_MEAN ; i++)
		{
			if (theCondMean.mvCondMean[i] != NULL)
			{
				theCondMean.mvCondMean[i]->Print(theOut);
				theOut << endl;
			}
		}
		return theOut ;
	}
#endif //_RDLL_

	/*!
	 * \fn cCondMean& cCondMean::operator =(cCondMean& theSrc)
	 * \param cCondMean& theSrc: source class
	 */
	cCondMean& cCondMean::operator =(cCondMean& theSrc)
	{
		Delete() ;
		for (uint i = 0; i < MAX_COND_MEAN; i++)
			mvCondMean[i] = NULL;
		for (uint i = 1; i < MAX_COND_MEAN; i++)
			if (theSrc.mvCondMean[i] != NULL)
			{
				mvCondMean[i] = CreateOneRealCondMean(*(theSrc.mvCondMean[i]));
			}
		return *this ;
	}

	void cCondMean::UpdateProxyMeanParameters(void)
	{
		for (uint i = 1; i < MAX_COND_MEAN; i++)
		{
			if (mvCondMean[i] != NULL)
				mvCondMean[i]->UpdateProxyMeanParameters();
		}
	}

	/*!
	 * \fn double cCondMean::ComputeMean(uint theDate, const cRegArchValue& theData) const
	 * \param int theDate: date of computation
	 * \param const cRegArchValue& theData: past datas.
	 * \details Compute the value of the conditional mean at date theDate. 
	 * theData is not updated here.
	 */
	double cCondMean::ComputeMean(uint theDate, const cRegArchValue& theData) const
	{
	double myMean = 0.0 ;
		for (uint i = 1 ; i < MAX_COND_MEAN ; i++)
			if (mvCondMean[i] != NULL)
				myMean += mvCondMean[i]->ComputeMean(theDate, theData) ;
		return myMean ;
	}

	uint cCondMean::GetNParam(void) const
	{
	uint myNParam = 0 ;
		for (uint i = 1 ; i < MAX_COND_MEAN ; i++)
			if (mvCondMean[i] != NULL)
				myNParam += mvCondMean[i]->GetNParam() ;
		return myNParam ;
	}

	uint cCondMean::GetNLags(void) const
	{
	uint myNLags = 0 ;
		for (uint i = 0 ; i < MAX_COND_MEAN ; i++)
			if (mvCondMean[i] != NULL)
				myNLags = MAX(myNLags, mvCondMean[i]->GetNLags());
		return myNLags ;
	}

	void cCondMean::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData)
	{
 	uint myIndex = 0 ;
		theGradData.mCurrentGradMu = 0.0L ;
		for (uint i = 0 ; i < MAX_COND_MEAN; i++)
			if (mvCondMean[i] != NULL)
			{	mvCondMean[i]->ComputeGrad(theDate, theValue, theGradData, myIndex) ;
				myIndex += mvCondMean[i]->GetNParam() ;
			}
	}

/*
	void cCondMean::NumericComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cAbstResiduals* theResiduals, cNumericDerivative& theNumDeriv)
	{
	uint myIndex = 0;
		theGradData.mCurrentGradMu = 0.0L;
		for (uint i = 0; i < GetNMean(); i++)
		{
			mvCondMean[i]->NumericComputeGrad(theDate, theData, theGradData, myIndex, theResiduals, theNumDeriv);
			myIndex += mvCondMean[i]->GetNParam();
		}
	}


	void cCondMean::NumericComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals, cNumericDerivative& theNumDeriv)
	{
	uint myIndex = 0;
		theGradData.mCurrentGradMu = 0.0L;
		for (uint i = 0; i < GetNMean(); i++)
		{
			mvCondMean[i]->NumericComputeGradAndHess(theDate, theData, theGradData, theHessData, myIndex, theResiduals, theNumDeriv);
			myIndex += mvCondMean[i]->GetNParam();
		}
	}
*/

	void cCondMean::RegArchParamToVector(cDVector& theDestVect, uint theIndex) const
	{
	uint myIndexCour = theIndex ;
		for (uint i = 0; i < MAX_COND_MEAN; i++)
			if (mvCondMean[i] != NULL)
			{	mvCondMean[i]->RegArchParamToVector(theDestVect, myIndexCour) ;
				myIndexCour += mvCondMean[i]->GetNParam() ;
			}
	}

	void cCondMean::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	uint myIndexCour = theIndex ;
		for (uint i = 0; i <MAX_COND_MEAN; i++)
			if (mvCondMean[i] != NULL)
			{	mvCondMean[i]->VectorToRegArchParam(theSrcVect, myIndexCour) ;
				myIndexCour += mvCondMean[i]->GetNParam() ;
			}
	}

	void cCondMean::ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData)
	{
	uint myIndex = 0;
		theHessData.mCurrentHessMu = 0.0;
		for (uint i = 0; i <MAX_COND_MEAN; i++)
			if (mvCondMean[i] != NULL)
			{
				mvCondMean[i]->ComputeHess(theDate, theData, theGradData, theHessData, myIndex);
				myIndex += mvCondMean[i]->GetNParam();
			}
	}

	void cCondMean::ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData)
	{
		uint myIndex = 0;
		theHessData.mCurrentHessMu = 0.0;
		for (uint i = 0; i < MAX_COND_MEAN; i++)
			if (mvCondMean[i] != NULL)
			{
				mvCondMean[i]->ComputeGradAndHess(theDate, theData, theGradData, theHessData, myIndex);
				myIndex += mvCondMean[i]->GetNParam();
			}
	}

	void cCondMean::GetParamName(char** theName)
	{
	uint myIndex = 0;
		for (uint i = 0; i < MAX_COND_MEAN; i++)
			if (mvCondMean[i] != NULL)
			{
				mvCondMean[i]->GetParamName(myIndex, theName);
				myIndex += mvCondMean[i]->GetNParam();	
			}
	
	}

	void cCondMean::GetParamName(string theName[])
	{
	uint myIndex = 0;
		for (uint i = 0; i < MAX_COND_MEAN; i++)
			if (mvCondMean[i] != NULL)
			{
				mvCondMean[i]->GetParamName(myIndex, theName);
				myIndex += mvCondMean[i]->GetNParam();
			}
	}



}//namespace
