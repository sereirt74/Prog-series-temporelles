#include "StdAfxRegArchLib.h"
/*!
\file cNumericDervative.cpp
\brief implementation of the cNumericDerivative methods

\author Jean-Baptiste DURAND, Ollivier TARAMASCO
\date jan-31-2018 - last change jan-31-2018 
*/

namespace RegArchLib {

	cNumericDerivative::cNumericDerivative(uint theNParam, uint theNDistrParam, double theh, uint theNObs)
	{
		mvNParam = theNParam;
		mvNDistrParam = theNDistrParam;
		mvh = theh;
		if (theNParam > 0)
		{
			mh1 = new double[theNParam];
			mh2 = new double*[theNParam];
			for (uint i = 0; i < theNParam; i++)
				mh2[i] = new double[theNParam];
			if (theNDistrParam > 0)
			{
				mLogDensForGrad.ReAlloc(theNDistrParam);
				mLogDensForHess.ReAlloc(theNDistrParam, theNParam);
				mGradDiffForGrad.ReAlloc(theNDistrParam);

			}
		}
		else
		{
			mh1 = NULL;
			mh2 = NULL;
		}
		if (theNParam == 0)
		{
			mValueForGrad = NULL;
			mValueForHess = NULL;
		}
		else
		{
			mValueForGrad = new cRegArchValue[theNParam];
			mValueForHess = new cRegArchValue*[theNParam];
			for (uint i = 0; i < theNParam; i++)
			{
				mValueForHess[i] = new cRegArchValue[theNParam];
				mValueForGrad[i].ReAlloc(theNObs);
				for (uint j = i; j < theNParam; j++)
					mValueForHess[i][j].ReAlloc(theNObs);
			}
		}
	}

	cNumericDerivative::cNumericDerivative(uint theNParam, uint theNDistrParam, double theh, cDVector& theYt, cDMatrix* theXt, cDMatrix *theXvt)
	{
		mvNParam = theNParam;
		mvNDistrParam = theNDistrParam;
		mvh = theh;
		if (theNParam > 0)
		{
			mh1 = new double[theNParam];
			mh2 = new double*[theNParam];
			for (uint i = 0; i < theNParam; i++)
				mh2[i] = new double[theNParam];
			if (theNDistrParam > 0)
			{
				mLogDensForGrad.ReAlloc(theNDistrParam);
				mLogDensForHess.ReAlloc(theNDistrParam, theNParam);
				mGradDiffForGrad.ReAlloc(theNDistrParam);
			}
		}
		else
		{
			mh1 = NULL;
			mh2 = NULL;
		}
		if (theNParam == 0)
		{
			mValueForGrad = NULL;
			mValueForHess = NULL;
		}
		else
		{
			mValueForGrad = new cRegArchValue[theNParam];
			mValueForHess = new cRegArchValue*[theNParam];
			for (uint i = 0; i < theNParam; i++)
			{
				mValueForHess[i] = new cRegArchValue[theNParam];
				mValueForGrad[i].ReAlloc(theYt);
				if (theXt != NULL)
					if (theXt->GetNCol() > 0)
						mValueForGrad[i].mXt = *theXt;
				if (theXvt != NULL)
					if (theXvt->GetNCol() > 0)
						mValueForGrad[i].mXvt = *theXvt;


				for (uint j = i; j < theNParam; j++)
				{
					mValueForHess[i][j].ReAlloc(theYt);
					if (theXt != NULL)
						if (theXt->GetNCol() > 0)
							mValueForHess[i][j].mXt = *theXt;
					if (theXvt != NULL)
						if (theXvt->GetNCol() > 0)
							mValueForHess[i][j].mXvt = *theXvt;

				}
			}
		}
	}

	cNumericDerivative::~cNumericDerivative()
	{
		if (mvNParam > 0)
		{
			for (uint i = 0 ; i < mvNParam; i++)
			{	mValueForGrad[i].Delete();
				for (uint j = 0; j < mvNParam; j++)
				{	mValueForHess[i][j].Delete();
				}
				delete[] mh2[i];
				delete[] mValueForHess[i];
			}
			delete[] mValueForGrad;
			delete[] mValueForHess;
			delete[] mh1;
			delete[] mh2;
			if (mvNDistrParam > 0)
			{
				mLogDensForGrad.Delete();
				mLogDensForHess.Delete();
				mGradDiffForGrad.Delete();

			}
		}
	}

	uint cNumericDerivative::GetNParam(void)
	{
		return mvNParam;

	}

	uint cNumericDerivative::GetNDistrParam(void)
	{
		return mvNDistrParam;

	}

	double cNumericDerivative::Geth(void)
	{
		return mvh;

	}



} // namespace