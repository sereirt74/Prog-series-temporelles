#include "StdAfxRegArchLib.h"
/*!
	\file cFuncMeanAndVar.cpp
	\brief Sources for cFuncMeanAndVar methods.
	\author Ollivier TARAMASCO
	\date apr-25-2023 - Last change feb-18-2011
*/

namespace RegArchLib {
/*
		cDMatrix md2Fxx;
		cDVector* md2Fxu;
		cDMatrix md2Fuu;
		cDVector* md2Fxh;
		cDVector* md2Fuh;
		cDMatrix md2Fhh;
*/
	cFuncMeanAndVar::cFuncMeanAndVar()
	{
		mNtheta = mNu = mNh = 0;
		md2Fxu = NULL;
		md2Fxh = NULL;
		md2Fuh = NULL;
		mCondMeanFuncBool = true;
	}

	cFuncMeanAndVar::cFuncMeanAndVar(uint theNtheta, bool theCondMeanBool, uint theNu, uint theNh)
	{
		mNtheta = theNtheta;
		mCondMeanFuncBool = theCondMeanBool;
		mNu = theNu;
		mNh = theNh;
		mdFx.ReAlloc(mNtheta);
		md2Fxx.ReAlloc(mNtheta, mNtheta);
		if (mNu > 0)
		{
			mdFu.ReAlloc(mNu);
			md2Fuu.ReAlloc(mNu, mNu);
			md2Fxu.ReAlloc(theNtheta, theNu);
			if (mNh > 0)
			{
				md2Fuh.ReAlloc(mNu, mNh) ;
			}
		}
		if (mNh > 0)
		{
			mdFh.ReAlloc(mNh);
			md2Fhh.ReAlloc(mNh, mNh);
			md2Fxh.ReAlloc(mNtheta, mNh);
		}
	}

	cFuncMeanAndVar::~cFuncMeanAndVar()
	{
		mdFu.Delete();
		mdFx.Delete();
		mdFh.Delete();
		md2Fxx.Delete();
		if (mNu > 0)
		{
			md2Fuu.Delete();
			if (mNh > 0)
			{
				md2Fuh.Delete();
			}
		}
		if (mNh > 0)
		{
			md2Fhh.Delete();
			md2Fxh.Delete();
		}
	}

	void cFuncMeanAndVar::ComputeGradFunc(cRegArchGradient& theGradData)
	{
		cDVector myRes;
		if (mCondMeanFuncBool)
			myRes = theGradData.mCurrentGradMu;
		else
			myRes = theGradData.mCurrentGradVar;
		myRes += mdFx;
		if (mNu > 0)
		{
			for (uint k = 0; k < mNu; k++)
				myRes -= mdFu[k] * theGradData.mGradMt[k];
		}

		if (mNh > 0)
		{
			if (mCondMeanFuncBool)
			{
				myRes += mdFh[0] * theGradData.mCurrentGradVar;
			}
			else
			{
				for (uint l = 0; l < mNh; l++)
				{
					myRes += mdFh[l] * theGradData.mGradHt[l];
				}
			}
		}

		if (mCondMeanFuncBool)
			theGradData.mCurrentGradMu = myRes;
		else
			theGradData.mCurrentGradVar = myRes;
	}

	void cFuncMeanAndVar::ComputeHessFunc(const cRegArchGradient& theGradData, cRegArchHessien& theHessData)
	{
		cDMatrix myRes = md2Fxx;
		
		if (mNu > 0)
			for (uint k = 0; k < mNu; k++)
				myRes -= mdFu[k] * theHessData.mHessMt[k];
		if (mNh > 0)
		{
			if (mCondMeanFuncBool)
				myRes += mdFh[0] * theHessData.mCurrentHessVar;
			else
				for (uint m = 0 ; m < mNh ; m++)
					myRes += mdFh[m] * theHessData.mHessHt[m];
		}

		for (uint i = 0 ; i < mNtheta ; i++)
			for (uint j = i; j < mNtheta; j++)
			{
				
				if (mNu > 0)
					for (uint k = 0; k < mNu; k++)
					{
						myRes[i][j] -= md2Fxu[j][k] * theGradData.mGradMt[k][i];
						myRes[i][j] -= md2Fxu[i][k] * theGradData.mGradMt[k][j];
					}

				if (mNh > 0)
				{
					if (mCondMeanFuncBool)
					{
						myRes[i][j] += md2Fxh[j][0] * theGradData.mCurrentGradVar[i];
						myRes[i][j] += md2Fxh[i][0] * theGradData.mCurrentGradVar[j];
					}

					else
						for (uint m = 0; m < mNh; m++)
						{
							myRes[i][j] += md2Fxh[j][m] * theGradData.mGradHt[m][i];
							myRes[i][j] += md2Fxh[i][m] * theGradData.mGradHt[m][j];
						}
				}
				if (mNu > 0)
				{
					for (uint k = 0; k < mNu; k++)
						for (uint l = 0; l < mNu; l++)
							myRes[i][j] += md2Fuu[k][l] * theGradData.mGradMt[k][i] * theGradData.mGradMt[l][j];
				}
				if (mNh > 0)
				{
					if (mCondMeanFuncBool)
						myRes[i][j] += md2Fhh[0][0] * theGradData.mCurrentGradVar[i] * theGradData.mCurrentGradVar[j];
					else
					{
						for (uint m = 0; m < mNh; m++)
							for (uint n = 0; n < mNh; n++)
								myRes[i][j] += md2Fhh[m][m] * theGradData.mGradHt[m][i] * theGradData.mGradHt[n][j];

					}
				}
				if (mNu > 0 && mNh > 0)
				{
					for (uint k = 0; k < mNu; k++)
					{
						if (mCondMeanFuncBool)
							myRes[i][j] -= md2Fuh[k][0] * theGradData.mGradMt[k][i] * theGradData.mCurrentGradVar[j];
						else
						{
							for (uint n = 0 ; n < mNh ; n++)
								myRes[i][j] -= md2Fuh[k][n] * theGradData.mGradMt[k][i] * theGradData.mGradHt[n][j];

						}
					}
					if (mCondMeanFuncBool)
					{
						for (uint l = 0; l < mNu; l++)
							myRes[i][j] -= md2Fuh[l][0] * theGradData.mGradMt[l][j] * theGradData.mCurrentGradVar[i];
					}
					else
					{
						for (uint l = 0; l < mNu; l++)
							for (uint n = 0 ; n < mNh ; n++)
								myRes[i][j] -= md2Fuh[l][n] * theGradData.mGradMt[l][j] * theGradData.mGradHt[n][i];

					}
				}
			}
		for (uint i = 0; i < mNtheta - 1; i++)
			for (uint j = i + 1; j < mNtheta; j++)
				myRes[j][i] = myRes[i][j];

		if (mCondMeanFuncBool)
			theHessData.mCurrentHessMu = myRes;
		else
			theHessData.mCurrentHessVar = myRes;

/*
		for (uint i = 0; i < mNtheta; i++)
		{
			for (uint j = i; j < mNtheta; j++)
			{
				myRes[i][j] = md2Fxx[i][j];
				for (uint m = 0; m < mNu; m++)
					myRes[i][j] -= md2Fxu[m][i] * theGradData.mGradMt[m][j];
				for (uint n = 0; n < mNh; n++)
				{
					if (mCondMeanFuncBool)
						myRes[i][j] += md2Fxh[n][i] * theGradData.mCurrentGradVar[j];
					else
						myRes[i][j] += md2Fxh[n][i] * theGradData.mGradHt[n][j];
				}
				for (uint k = 0; k < mNu; k++)
					myRes[i][j] -= mdFu[k] * theHessData.mHessMt[k][i][j];
				for (uint l = 0; l < mNh; l++)
				{
					if (mCondMeanFuncBool)
						myRes[i][j] += mdFh[l] * theHessData.mCurrentHessVar[i][j];
					else
						myRes[i][j] += mdFh[l] * theHessData.mHessHt[l][i][j];
				}
				for (uint k = 0; k < mNu; k++)
				{
					double myAux = md2Fxu[k][j];
					for (uint m = 0; m < mNu; m++)
						myAux -= md2Fuu[m][k] * theGradData.mGradMt[m][j];
					for (uint n = 0; n < mNh; n++)
					{
						if (mCondMeanFuncBool)
							myAux += md2Fuh[k][n] * theGradData.mCurrentGradVar[j];
						else
							myAux += md2Fuh[k][n] * theGradData.mGradHt[n][j];
					}
					myAux *= theGradData.mGradMt[k][i];
					myRes[i][j] -= myAux;
				}
				for (uint l = 0; l < mNh; l++)
				{
					double myAux = md2Fxh[l][j];
					for (uint m = 0; m < mNu; m++)
						myAux -= md2Fuh[m][l] * theGradData.mGradMt[m][j];
					for (uint n = 0; n < mNh; n++)
					{
						if (mCondMeanFuncBool)
							myAux += md2Fhh[l][n] * theGradData.mCurrentGradVar[j];
						else
							myAux += md2Fhh[l][n] * theGradData.mGradHt[n][j];
					}
					if (mCondMeanFuncBool)
						myAux *= theGradData.mCurrentGradVar[i]; 
					else
						myAux *= theGradData.mGradHt[l][i];
					myRes[i][j] += myAux;
				}
			}
		}
		for (uint i = 0; i < mNtheta - 1; i++)
			for (uint j = i + 1; j < mNtheta; j++)
				myRes[j][i] = myRes[i][j];

		if (mCondMeanFuncBool)
			theHessData.mCurrentHessMu = myRes;
		else
			theHessData.mCurrentHessVar = myRes;

*/
	}

	void cFuncMeanAndVar::ComputeGradAndHessFunc(cRegArchGradient& theGradData, cRegArchHessien& theHessData)
	{
		ComputeGradFunc(theGradData);
		ComputeHessFunc(theGradData, theHessData);
	}

} //namespace