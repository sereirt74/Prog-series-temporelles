#include "StdAfxRegArchLib.h"
/*!
	\file cAbstCondMean.cpp
	\brief Sources for abstract class cAbstCondMean methods.
	\author Jean-Baptiste DURAND, Ollivier TARAMASCO
	\date dec-18-2006 - Last change feb-18-2011
*/

namespace RegArchLib {

	/*!
	 * \fn cAbstCondMean::cAbstCondMean(const eCondMeanEnum theType)
	 * \param const eCondMeanEnum theType. Code of the conditional mean type. Default value: eUnknown
	 * \details Set the real type code of conditional mean component
	 */
	cAbstCondMean::cAbstCondMean(const eCondMeanEnum theType)
	{
		mvCondMeanType = theType ;
  		MESS_CREAT("cAbstCondMean")
	}

	/*!
	 * \fn cAbstCondMean::~cAbstCondMean()
	 */
	cAbstCondMean::~cAbstCondMean()
	{
		mvCondMeanType = eUnknown ;
  		MESS_DESTR("cAbstCondMean")
	}

	/*!
	 * \fn  inline eCondMeanEnum cAbstCondMean::GetCondMeanType(void)
	 * \param void
	 * \return a eCondMeanEnum value.
	 */
	eCondMeanEnum cAbstCondMean::GetCondMeanType(void) const
	{
		return mvCondMeanType ;
	}

	/*!
	 * \fn  void cAbstCondMean::SetCondMeanType(eCondMeanEnum theType)
	 * \param eCondMeanEnum theType
	 * \details mvCondMeanType = theType
	 */
	void cAbstCondMean::SetCondMeanType(eCondMeanEnum theType)
	{	mvCondMeanType = theType ;
	}

	void cAbstCondMean::ComputeGradWithM(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, uint theBegIndex, cAbstResiduals* theResiduals)
	{
/*		double myGradH = 0;
		uint myNu = 0;
		uint myNh = 0;
		uint myNAllParam = theGradData.mCurrentGradMu.GetSize();
		cDVector myGradTheta = cDVector(myNAllParam);
		uint myNPast = GetNLags();
		cDVector myGradU = NULL;
		if (myNPast > 0)
			myGradU.ReAlloc(myNPast);
		ComputeGradForM(theDate, theValue, theBegIndex, myGradTheta, myGradU, myGradH, myNu, myNh);
		theGradData.mCurrentGradMu += myGradTheta;
		for (uint k = 0; k < myNu; k++)
			theGradData.mCurrentGradMu -= myGradU[k] * theGradData.mGradMt[k];
		if (myNh > 0)
			theGradData.mCurrentGradMu += myGradH * theGradData.mCurrentGradVar;
*/
		uint myNu = GetNu();
		uint myNh = GetNh();
		uint myNtheta = theGradData.mCurrentGradMu.GetSize();
		cFuncMeanAndVar myDerivM = cFuncMeanAndVar(myNtheta, true, myNu, myNh);
		ComputeGradForM(theDate, theValue, theBegIndex, myDerivM);
		myDerivM.ComputeGradFunc(theGradData);
/*
		theGradData.mCurrentGradMu += myDerivM.mdFx;
		for (uint k = 0; k < myNu; k++)
			theGradData.mCurrentGradMu -= myDerivM.mdFu[k] * theGradData.mGradMt[k];
		if (myNh > 0)
			theGradData.mCurrentGradMu += myDerivM.mdFh[0] * theGradData.mCurrentGradVar;
*/
	}

	void cAbstCondMean::ComputeGradWithF(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, uint theBegIndex, cAbstResiduals* theResiduals)
	{
	uint myNParam[3];
		GetNParamF(myNParam);
	cDVector myGradF(myNParam[0] + myNParam[1] + myNParam[2]);
		ComputeGradF(theDate, theData, myGradF);
		for (uint i = 0; i < myNParam[0]; i++)
			theGradData.mCurrentGradMu[theBegIndex + i] += myGradF[i];
		for (uint k = 0; k < MIN(myNParam[1], theDate); k++)
			theGradData.mCurrentGradMu -= myGradF[myNParam[0] + k] * theGradData.mGradMt[k];
		if (myNParam[2] > 0)
			theGradData.mCurrentGradMu += myGradF[myNParam[0] + myNParam[1]] * theGradData.mCurrentGradSigma;
	}

	void cAbstCondMean::ComputeHessWithF(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex, cAbstResiduals* theResiduals)
	{
	uint myNParam[3];
		GetNParamF(myNParam);
	uint myN = myNParam[0] + myNParam[1] + myNParam[2];
	uint myNm1 = myNParam[0] + myNParam[1];
	uint myNParamTot = theGradData.GetNParam();
	cDVector myGradF(myN);
 		ComputeGradF(theDate, theData, myGradF);
	cDMatrix myHessF(myN, myN);
		ComputeHessF(theDate, theData, myHessF);
	cDMatrix myD2Fxx(myNParamTot, myNParamTot);
		for (uint i = 0; i < myNParam[0]; i++)
			for (uint j = 0; j < myNParam[0]; j++)
				myD2Fxx[i + theBegIndex][j + theBegIndex] = myHessF[i][j];
	cDVector* myD2Fxy = NULL;
		if (myNParam[1] > 0)
			myD2Fxy = new cDVector[myNParam[1]];
		for (uint i = 0; i < myNParam[1]; i++)
		{
			myD2Fxy[i].ReAlloc(myNParamTot);
			for (uint j = 0; j < myNParam[0]; j++)
				myD2Fxy[i][j + theBegIndex] = myHessF[i + myNParam[0]][j];
		}
	cDVector myD2Fxz;
		if (myNParam[2] > 0)
		{
			myD2Fxz.ReAlloc(myNParamTot);
			for (uint j = 0; j < myNParam[0]; j++)
				myD2Fxz[j + theBegIndex] = myHessF[myNm1][j];
		}
	cDMatrix myD2Fyy = cDMatrix(myNParam[1], myNParam[1]);
		for (uint i = 0; i < myNParam[1]; i++)
			for (uint j = 0; j < myNParam[1]; j++)
				myD2Fyy[i][j] = myHessF[i + myNParam[0]][j + myNParam[0]];
	double* myD2Fyz = NULL;
		if (myNParam[1] > 0)
		{
			myD2Fyz = new double[myNParam[1]];
			for (uint k = 0; k < myNParam[1]; k++)
			{
				myD2Fyz[k] = 0;
				if (myNParam[2] > 0)
					myD2Fyz[k] = myHessF[myNm1][k + myNParam[0]];
			}
		}
		theHessData.mCurrentHessMu += myD2Fxx;
		for (uint k = 0; k < myNParam[1]; k++)
			theHessData.mCurrentHessMu -= myD2Fxy[k] * Transpose(theGradData.mGradMt[k]);
		if (myNParam[2] > 0)
			theHessData.mCurrentHessMu += myD2Fxz * Transpose(theGradData.mCurrentGradSigma);

		for (uint k = 0; k < myNParam[1]; k++)
		{
		cDVector myAux = myD2Fxy[k];
			for (uint l = 0; l < myNParam[1]; l++)
				myAux -= myD2Fyy[k][l] * theGradData.mGradMt[l];
			if (myNParam[2] > 0)
				myAux += myD2Fyz[k] * theGradData.mCurrentGradSigma;
			theHessData.mCurrentHessMu -= myAux * Transpose(theGradData.mGradMt[k]);
			theHessData.mCurrentHessMu -= myGradF[myNParam[0] + k] * theHessData.mHessMt[k];
		}
		if (myNParam[2] > 0)
		{
		cDVector myAux = myD2Fxz;
			for (uint l = 0; l < myNParam[1]; l++)
				myAux -= myHessF[myNm1][myNParam[0]+l] * theGradData.mGradMt[l];
			myAux += myHessF[myNm1][myNm1] * theGradData.mCurrentGradSigma;
			theHessData.mCurrentHessMu += myAux * Transpose(theGradData.mCurrentGradSigma);
			theHessData.mCurrentHessMu += myHessF[myNm1][myNm1] * theHessData.mCurrentHessSigma;
		}
		if (myNParam[1] > 0)
		{
			for (uint i = 0; i < myNParam[1]; i++)
				myD2Fxy[i].Delete();
			delete myD2Fxy;
			delete myD2Fyz;
		}
	}

	void cAbstCondMean::ComputeHessWithM(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex, cAbstResiduals* theResiduals)
	{
/*
		double myGradH = 0;
		uint myNu = 0;
		uint myNh = 0;
		uint myNAllParam = theGradData.mCurrentGradMu.GetSize();
		cDVector myGradTheta = cDVector(myNAllParam);
		uint myNPast = GetNLags();
		cDVector myGradU = NULL;
		if (myNPast > 0)
			myGradU.ReAlloc(myNPast);
		ComputeGradForM(theDate, theData, theBegIndex, myGradTheta, myGradU, myGradH, myNu, myNh);
		cDMatrix myHessTheta2 = cDMatrix(myNAllParam, myNAllParam);
		cDVector* myHessThetaU = NULL;
		cDVector myHessThetaH = NULL;
		cDVector* myHessU2 = NULL;
		cDVector myHessUH = NULL;
		double myHessH2 = 0;
		if (myNPast > 0)
		{
			myHessThetaU = new cDVector[myNPast];
			myHessU2 = new cDVector[myNPast];
			for (uint i = 0; i < myNPast; i++)
			{
				myHessThetaU[i].ReAlloc(myNAllParam);
				myHessU2[i].ReAlloc(myNPast);
			}
			if (myNh > 0)
				myHessUH.ReAlloc(myNPast);
		}
		if (myNh > 0)
		{
			myHessThetaH.ReAlloc(myNAllParam);
		}

		ComputeHessForM(theDate, theData, theBegIndex, myGradTheta, myGradU, myGradH, myNu, myNh, myHessTheta2, myHessThetaU, myHessThetaH, myHessU2, myHessUH, myHessH2);
		theHessData.mCurrentHessMu += myHessTheta2;
		if (myNu > 0)
		{
			for (uint k = 0; k < myNu; k++)
			{
				theHessData.mCurrentHessMu -= myHessThetaU[k] * Transpose(theGradData.mGradMt[k]);
				for (uint l = 0 ; l < myNu ; l++)
					theHessData.mCurrentHessMu -= myHessU2[k][l] * Transpose(theGradData.mGradMt[k]);
				if (myNh > 0)
					theHessData.mCurrentHessMu -= myHessUH * Transpose(theGradData.mGradMt[k]);
				theHessData.mCurrentHessMu -= myGradU[k] * theHessData.mHessMt[k];
			}
		}

		if (myNh > 0)
		{
			theHessData.mCurrentHessMu += myHessThetaH * Transpose(theGradData.mCurrentGradVar);
			for (uint k = 0; k < myNu; k++)
				theHessData.mCurrentHessMu -= theGradData.mGradMt[k] * Transpose(myHessUH);

			theHessData.mCurrentHessMu += theGradData.mCurrentGradVar * Transpose(myHessThetaH);
			theHessData.mCurrentHessMu += myGradH * theHessData.mCurrentHessVar;
		}
		if (myNPast > 0)
		{
			for (uint i = 0; i < myNPast; i++)
			{
				myHessThetaU[i].Delete();
				myHessU2[i].Delete();
			}
			delete[] myHessThetaU;
			delete[] myHessU2;
		}
*/
		uint myNu = GetNu();
		uint myNh = GetNh();
		uint myNtheta = theGradData.mCurrentGradMu.GetSize();
		cFuncMeanAndVar myDerivM = cFuncMeanAndVar(myNtheta, true, myNu, myNh);
		ComputeGradForM(theDate, theData, theBegIndex, myDerivM);
		ComputeHessForM(theDate, theData, theBegIndex, myDerivM);
		myDerivM.ComputeHessFunc(theGradData, theHessData);
	}

	void cAbstCondMean::ComputeGradAndHessWithM(uint theDate, uint theBegIndex, const cRegArchValue& theValue, cRegArchGradient& theGradData, cRegArchHessien& theHessData)
	{
/*
		uint myNAllParam = theGradData.mCurrentGradMu.GetSize();
		uint myNPast = theGradData.GetNPast();
		cDVector myGradTheta = cDVector(myNAllParam);
		cDVector myGradU = cDVector(myNPast);
		double myGradH = 0;
		uint myNu = 0, myNh = 0;
		ComputeGradForM(theDate, theValue, theBegIndex, myGradTheta, myGradU, myGradH, myNu, myNh);
		theGradData.mCurrentGradMu += myGradTheta;
		for (uint k = 0; k < myNu; k++)
			theGradData.mCurrentGradMu -= myGradU[k] * theGradData.mGradMt[k];
		if (myNh > 0)
			theGradData.mCurrentGradMu += myGradH * theGradData.mCurrentGradVar;
		cDMatrix myHessTheta2 = cDMatrix(myNAllParam, myNAllParam);
		cDVector* myHessThetaU = NULL;
		cDVector myHessThetaH = NULL;
		cDVector* myHessU2 = NULL;
		cDVector myHessUH = NULL;
		double myHessH2 = 0;
		if (myNPast > 0)
		{
			myHessThetaU = new cDVector[myNPast];
			myHessU2 = new cDVector[myNPast];
			for (uint i = 0; i < myNPast; i++)
			{
				myHessThetaU[i].ReAlloc(myNAllParam);
				myHessU2[i].ReAlloc(myNPast);
			}
			if (myNh > 0)
				myHessUH.ReAlloc(myNPast);
		}
		if (myNh > 0)
		{
			myHessThetaH.ReAlloc(myNAllParam);
		}

		ComputeHessForM(theDate, theValue, theBegIndex, myGradTheta, myGradU, myGradH, myNu, myNh, myHessTheta2, myHessThetaU, myHessThetaH, myHessU2, myHessUH, myHessH2);
		theHessData.mCurrentHessMu += myHessTheta2;
		if (myNu > 0)
		{
			for (uint k = 0; k < myNu; k++)
			{
				theHessData.mCurrentHessMu -= myHessThetaU[k] * Transpose(theGradData.mGradMt[k]);
				theHessData.mCurrentHessMu -= theGradData.mGradMt[k] * Transpose(myHessThetaU[k]);
				theHessData.mCurrentHessMu -= myGradU[k] * theHessData.mHessMt[k];
			}
		}

		if (myNh > 0)
		{
			theHessData.mCurrentHessMu += myHessThetaH * Transpose(theGradData.mCurrentGradVar);
			theHessData.mCurrentHessMu += theGradData.mCurrentGradVar * Transpose(myHessThetaH);
			theHessData.mCurrentHessMu += myGradH * theHessData.mCurrentHessVar;
		}

		if (myNPast > 0)
		{
			for (uint i = 0; i < myNPast; i++)
			{
				myHessThetaU[i].Delete();
				myHessU2[i].Delete();
			}
			delete[] myHessThetaU;
			delete[] myHessU2;
		}
*/
		uint myNu = GetNu();
		uint myNh = GetNh();
		uint myNtheta = theGradData.mCurrentGradMu.GetSize();
		cFuncMeanAndVar myDerivM = cFuncMeanAndVar(myNtheta, true, myNu, myNh);
		ComputeGradForM(theDate, theValue, theBegIndex, myDerivM);
		myDerivM.ComputeGradFunc(theGradData);
		ComputeHessForM(theDate, theValue, theBegIndex, myDerivM);
		myDerivM.ComputeHessFunc(theGradData, theHessData);

	}

	void cAbstCondMean::ComputeGradAndHessWithF(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex, cAbstResiduals* theResiduals)
	{
	uint myNParam[3];
		GetNParamF(myNParam);
	uint myN = myNParam[0] + myNParam[1] + myNParam[2];
	uint myNm1 = myNParam[0] + myNParam[1];
	uint myNParamTot = theGradData.GetNParam();

	cDVector myGradF(myN);
		ComputeGradF(theDate, theData, myGradF);
		for (uint i = 0; i < myNParam[0]; i++)
			theGradData.mCurrentGradMu[theBegIndex + i] += myGradF[i];
		for (uint k = 0; k < MIN(myNParam[1], theDate); k++)
			theGradData.mCurrentGradMu -= myGradF[myNParam[0] + k] * theGradData.mGradMt[k];
		if (myNParam[2] > 0)
			theGradData.mCurrentGradMu += myGradF[myNParam[0] + myNParam[1]] * theGradData.mCurrentGradSigma;
		
	cDMatrix myHessF(myN, myN);
		ComputeHessF(theDate, theData, myHessF);
	cDMatrix myD2Fxx(myNParamTot, myNParamTot);
	for (uint i = 0; i < myNParam[0]; i++)
		for (uint j = 0; j < myNParam[0]; j++)
			myD2Fxx[i + theBegIndex][j + theBegIndex] = myHessF[i][j];
	cDVector* myD2Fxy = NULL;
		if (myNParam[1] > 0)
			myD2Fxy = new cDVector[myNParam[1]];
		for (uint i = 0; i < myNParam[1]; i++)
		{
			myD2Fxy[i].ReAlloc(myNParam[0]);
			for (uint j = 0; j < myNParam[0]; j++)
				myD2Fxx[i][j + theBegIndex] = myHessF[i + myNParam[0]][j];
		}
	cDVector myD2Fxz ;
		if (myNParam[2] > 0)
		{
			myD2Fxz.ReAlloc(myNParamTot);
			for (uint j = 0; j < myNParam[0]; j++)
				myD2Fxz[j+theBegIndex] = myHessF[myNm1][j];
		}
	cDMatrix myD2Fyy = cDMatrix(myNParam[1], myNParam[1]);
		for (uint i = 0; i < myNParam[1]; i++)
			for (uint j = 0; j < myNParam[1]; j++)		
				myD2Fyy[i][j] = myHessF[i + myNParam[0]][j + myNParam[0]];
	double* myD2Fyz = NULL;
		if (myNParam[1] > 0)
		{
			myD2Fyz = new double[myNParam[1]];
			for (uint k = 0; k < myNParam[1]; k++)
			{
				myD2Fyz[k] = 0;
				if (myNParam[2] > 0)
					myD2Fyz[k] = myHessF[myNm1][k + myNParam[0]];
			}
		}
		theHessData.mCurrentHessMu += myD2Fxx;
		for (uint k = 0; k < myNParam[1]; k++)
			theHessData.mCurrentHessMu -= myD2Fxy[k] * Transpose(theGradData.mGradMt[k]);
		if (myNParam[2] > 0)
			theHessData.mCurrentHessMu += myD2Fxz * Transpose(theGradData.mCurrentGradSigma);
		
		for (uint k = 0; k < myNParam[1]; k++)
		{
		cDVector myAux = myD2Fxy[k];
			for (uint l = 0; l < myNParam[1]; l++)
				myAux -= myD2Fyy[k][l] * theGradData.mGradMt[l];
			if (myNParam[2] > 0)
				myAux += myD2Fyz[k] * theGradData.mCurrentGradSigma;
			theHessData.mCurrentHessMu -= myAux * Transpose(theGradData.mGradMt[k]);
			theHessData.mCurrentHessMu -= myGradF[myNParam[0] + k] * theHessData.mHessMt[k];
		}
		if (myNParam[2] > 0)
		{
		cDVector myAux = myD2Fxz;
			for (uint l = 0; l < myNParam[1]; l++)
				myAux -= myHessF[myNm1][myNParam[0]] * theGradData.mGradMt[l];
			myAux += myHessF[myNm1][myNm1] * theGradData.mCurrentGradSigma;
			theHessData.mCurrentHessMu += myAux * Transpose(theGradData.mCurrentGradSigma);
			theHessData.mCurrentHessMu += myHessF[myNm1][myNm1] * theHessData.mCurrentHessSigma;
		}
		if (myNParam[1] > 0)
		{
			for (uint i = 0; i < myNParam[1]; i++)
				myD2Fxy[i].Delete();
			delete myD2Fxy;
			delete myD2Fyz;
		}

	}

/*
	void cAbstCondMean::NumericComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, uint theBegIndex, cAbstResiduals* theResiduals, cNumericDerivative& theNumDeriv)
	{
	uint myNParam = GetNParam();
	cDVector myParam(myNParam);
		RegArchParamToVector(myParam);
	double myh0 = theNumDeriv.Geth();
	double myF0 = ComputeMean(theDate, theData);
		for (uint i = 0; i < myNParam; i++)
		{
		double myh = fabs(myh0*myParam[i]);
			if (myh < 1e-16)
				myh = myh0;
			myParam[i] += myh;
			VectorToRegArchParam(myParam);
		double myF1 = ComputeMean(theDate, theNumDeriv.mValueForGrad[i+theBegIndex]);
			theGradData.mCurrentGradMu[theBegIndex + i] = (myF1 - myF0) / myh;
			myParam[i] -= myh;
		}
	}

	void cAbstCondMean::NumericComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex, cAbstResiduals* theResiduals, cNumericDerivative& theNumDeriv)
	{

		uint myNParamTot = theGradData.GetNParam();
		uint myNParam = GetNParam();
		cDVector myParam(myNParamTot);
			RegArchParamToVector(myParam, theBegIndex);
		double myh0 = theNumDeriv.Geth();
		double myF0 = theData.mMt[theDate];
		double* myF1 = new double[myNParam];
		double* myh1 = new double[myNParam];
		for (uint i = 0; i < myNParam; i++)
		{
			myh1[i] = fabs(myh0*myParam[i+theBegIndex]);
			if (myh1[i] < 1e-16)
				myh1[i] = myh0;
			myF1[i] = theNumDeriv.mValueForGrad[i + theBegIndex].mMt[theDate];
			theGradData.mCurrentGradMu[theBegIndex + i] = (myF1[i] - myF0) / myh1[i];
		}
		for (uint i = 0 ; i < myNParam ; i++)
		{	for (uint j = i; j < myNParam; j++)
			{
			double myh2;
				if (j > i)
					myh2 = fabs(myh0*myParam[j]);
				else
					myh2 = fabs(myh0*(myParam[j]-myh1[j]));
				if (myh2 < 1e-16)
					myh2 = myh0;
				double myF2 = theNumDeriv.mValueForHess[i+theBegIndex][j+theBegIndex].mMt[theDate];
				theHessData.mCurrentHessMu[theBegIndex+i][theBegIndex+j] = theHessData.mCurrentHessMu[theBegIndex+j][theBegIndex+i] = (myF2 - myF1[i] - myF1[j] + myF0) / (myh1[i]*myh2);
			}
		}
		delete[] myh1;
		delete[] myF1;
	}
*/

#ifndef _RDLL_
	/*!
	 * \fn ostream& operator <<(ostream& theOut, const cAbstCondMean& theAbstCondMean)
	 * \param ostream& theOut: output stream (file or screen). Default cout.
	 * \param const cAbstCondMean& theAbstCondMean: the mean component model.
	 * \details Uses cAbstCondMean::Print method.
	 */
	ostream& operator <<(ostream& theOut, const cAbstCondMean& theAbstCondMean)
	{
		theAbstCondMean.Print(theOut) ;
		return theOut ;
	}
#endif // _RDLL_

	template<class T>
	T* TemplateCreateOneRealCondMean(void)
	{
//	static std::unique_ptr<T> myCondMean(new T());
//		return &*myCondMean;
		return new T();
	}

	template<class T>
	T* TemplateCreateOneRealCondMean(cAbstCondMean* theAbstCondMean)
	{
	T* mySrc = static_cast<T *>(theAbstCondMean);
		if (mySrc)
		{
		T* myCondMean;
			myCondMean = new T(*mySrc);
			return myCondMean;
	
//	static std::shared_ptr<T> myCondMean(new T(*mySrc));
//			return &*myCondMean;
		}
		else
			throw cError("Wrong contional mean in TemplateCreateOneRealCondMean");
	}

	/*!
	* \fn cAbstCondMean* CreateOneRealCondMean(const eCondMeanEnum theType)
	* \param const eCondMeanEnum theType: code of the real condtional mean component.
	* \details This function has to be changed when adding a new conditional mean type.
	*/
	cAbstCondMean* CreateOneRealCondMean(eCondMeanEnum theType)
	{
		switch (theType)
		{
		case eConst:
			return TemplateCreateOneRealCondMean<cConst>();
			break;
		case eAr:
			return TemplateCreateOneRealCondMean<cAr>();
			break;
		case eMa:
			return TemplateCreateOneRealCondMean<cMa>();
			break;
		case eLinReg:
			return TemplateCreateOneRealCondMean<cLinReg>();
			break;
		case eStdDevInMean:
			return TemplateCreateOneRealCondMean<cStdDevInMean>();
			break;
		case eVarInMean:
			return TemplateCreateOneRealCondMean<cVarInMean>();
			break;
		case eArfima:
			return TemplateCreateOneRealCondMean<cArfima>();
			break;
		default:
			throw cError("CreateOneRealCondMean: unknown conditional mean type");
			break;
		}
	}

	cAbstCondMean* CreateOneRealCondMean(cAbstCondMean& theAbstCondMean)
	{
		switch (theAbstCondMean.GetCondMeanType())
		{
		case eConst:
			return TemplateCreateOneRealCondMean<cConst>(&theAbstCondMean);
			break;
		case eAr:
			return TemplateCreateOneRealCondMean<cAr>(&theAbstCondMean);	
			break;
		case eMa:
			return TemplateCreateOneRealCondMean<cMa>(&theAbstCondMean);
			break;
		case eLinReg:
			return TemplateCreateOneRealCondMean<cLinReg>(&theAbstCondMean);
			break;
		case eStdDevInMean:		
			return TemplateCreateOneRealCondMean<cStdDevInMean>(&theAbstCondMean);
			break;
		case eVarInMean:
			return TemplateCreateOneRealCondMean<cVarInMean>(&theAbstCondMean);
			break;
		case eArfima:
			return TemplateCreateOneRealCondMean<cArfima>(&theAbstCondMean);
			break;
		default:
			throw cError("CreateOneRealCondMean: unknown conditional mean type");
			break;
		}
	}
	
	cAbstCondMean* cAbstCondMean::PtrCopy(void) 	
	{
		return CreateOneRealCondMean(*this);
	}
	
}//namespace
