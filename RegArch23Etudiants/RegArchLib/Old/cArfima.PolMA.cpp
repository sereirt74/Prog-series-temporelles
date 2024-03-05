#include "StdAfxRegArchLib.h"
/*!
\file cArfima.cpp
\brief sources for class cAr methods.

\author Jean-Baptiste Durand, Ollivier TARAMASCO
\date feb-18-2016 - Last change feb-18-2016
*/

namespace RegArchLib {
	/*!
	* \fn cArfimar(uint theNAr = 0, uint theNMa = 0, double theFracD = 0, uint theNTruncLag = 100):cAbstCondMean(eArfima)
	* \param const uint theNAr: number of AR lags.
	*/
	
	cArfima::cArfima(uint theNAr, uint theNMa, double theFracD , uint theNTruncLag) :cAbstCondMean(eArfima)
	{
		mvAr.ReAlloc(theNAr);
		mvMa.ReAlloc(theNMa);
		mvFracD = theFracD;
		mvNTruncLag = theNTruncLag;
		mvPolMa.Resize(mvNTruncLag);
	uint myNParam = theNAr + theNMa + 1;
		
		mvGradPolMa = new cPolynome[myNParam];
		mvHessPolMa = new cPolynome*[myNParam];
		for (uint i = 0; i < myNParam; i++)
		{
			mvGradPolMa[i].Resize(mvNTruncLag);
			mvHessPolMa[i] = new cPolynome[myNParam];
			for (uint j = 0; j < myNParam; j++)
				mvHessPolMa[i][j].Resize(mvNTruncLag);
		}
		MESS_CREAT("cArfima")
	}

	cArfima::cArfima(const cDVector& theAr, const cDVector& theMa, double theFracD, uint theNTruncLag) :cAbstCondMean(eArfima)
	{
		mvAr = theAr;
		mvMa = theMa;
		mvFracD = theFracD;
		mvNTruncLag = theNTruncLag;
		mvPolMa.Resize(mvNTruncLag);
		uint myNParam = theAr.GetSize() + theMa.GetSize() + 1;

		mvGradPolMa = new cPolynome[myNParam];
		mvHessPolMa = new cPolynome*[myNParam];
		for (uint i = 0; i < myNParam; i++)
		{
			mvGradPolMa[i].Resize(mvNTruncLag);
			mvHessPolMa[i] = new cPolynome[myNParam];
			for (uint j = 0; j < myNParam; j++)
				mvHessPolMa[i][j].Resize(mvNTruncLag);
		}
		MESS_CREAT("cArfima")
	}

	cArfima::cArfima(const cArfima& theArfima) :cAbstCondMean(eUnknown)
	{	
		*this = theArfima;
		MESS_CREAT("cAr")
	}

	cArfima:: ~cArfima()
	{
	uint myNParam = mvAr.GetSize() + mvMa.GetSize() + 1;
		mvAr.Delete();
		mvMa.Delete();
		mvPolMa.Delete();
		for (uint i = 0; i < myNParam; i++)
		{
			mvGradPolMa[i].Delete();
			for (uint j = 0; j < myNParam; j++)
				mvHessPolMa[i][j].Delete();
			delete[] mvHessPolMa[i];
		}
		delete[] mvGradPolMa;
		delete[] mvHessPolMa;

			MESS_DESTR("cArfima")
	}

/*
cAbstCondMean* cArfima::PtrCopy(void)
	{
//	cArfima* myArfima = new cArfima(*this);

//		return myArfima;
		return cAbstCondMeanPtrCopy<cArfima>();
	}
*/
	void cArfima::Delete(void)
	{
	uint myNParam = mvAr.GetSize() + mvMa.GetSize() + 1;
		mvAr.Delete();
		mvMa.Delete();
		mvPolMa.Delete();
		for (uint i = 0; i < myNParam; i++)
		{
			mvGradPolMa[i].Delete();
			for (uint j = 0; j < myNParam; j++)
				mvHessPolMa[i][j].Delete();
			delete[] mvHessPolMa[i];
		}
		delete[] mvGradPolMa;
		delete[] mvHessPolMa;
	}

#ifndef _RDLL_
	void cArfima::Print(ostream& theOut) const
	{
	uint myNAr = mvAr.GetSize();
	uint myNMa = mvMa.GetSize();
		theOut << "ARFIMA(" << myNAr << ", d, " << myNMa << ") model with:" << endl;
		for (uint i = 0; i < myNAr; i++)
			theOut << "\tAR[" << i + 1 << "]=" << mvAr[i] << endl;
		theOut << "\td=" << mvFracD << endl;
		for (uint i = 0; i < myNMa; i++)
			theOut << "\tMA[" << i + 1 << "]=" << mvMa[i] << endl;
	}
#else
	void cArfima::Print(void)
	{
	uint myNAr = mvAr.GetSize();
	uint myNMa = mvMa.GetSize();
		Rprintf("ARFIMA(%d, d, ) model with:\n", myNAr, myNMa);
		for (uint i = 0; i < myNAr; i++)
			Rprintf("\tAR[%d]=%f\n", i + 1, mvAr[i]);
		Rprintf("\td=%f\n", mvFracD);
		for (uint i = 0; i < myNMa; i++)
			Rprintf("\tMA[%d]=%f\n", i + 1, mvMa[i]);
	}
#endif // _RDLL_
	void  cArfima::SetDefaultInitPoint(double theMean, double theVar)
	{
		mvAr = 0.1;
		mvMa = 0.1;
		mvFracD = 0.5;

	}

	void  cArfima::SetDefaultInitPoint(cRegArchValue& theValue)
	{
		mvAr = 0.1;
		mvMa = 0.1;
		mvFracD = 0.5;

	}

	void  cArfima::ReAlloc(const uint theSize, const uint theNumParam)
	{
		switch (theNumParam)
		{
		case 0:
			mvAr.ReAlloc(theSize);
			break;
		case 1:
			mvMa.ReAlloc(theSize);
			break;
		default:
			throw cError("cArfima::ReAlloc - theNumParam must be in 0 or 1");
			break;
		}
	}
	
	void  cArfima::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{
		switch (theNumParam)
		{
		case 0: // mvAr
			mvAr = theVectParam;
			break;
		case 1: // mvMa
			mvMa = theVectParam;
			break;
		default:
			throw cError("cArfima::ReAlloc - theNumParam must be in 0 or 1");
			break;
		}

	}
	
	void  cArfima::Set(const double theValue, const uint theIndex, const uint theNumParam)
	{
		switch (theNumParam)
		{
		case 0: // AR
			if (theIndex < mvAr.GetSize())
				mvAr[theIndex] = theValue;
			else
				throw cError("cArfima::Set - wrong index");
			break;
		case 1:
			if (theIndex < mvMa.GetSize())
				mvMa[theIndex] = theValue;
			else
				throw cError("cArfima::Set - wrong index");
			break;
		case 2:
			mvFracD = theValue;
			break;

		default:
			throw cError("cArfima::Set - theNumParam must be in 0, 1, 2");
			break;

		}
	}
	
	void  cArfima::Set(const cDVector& theVectParam, const uint theNumParam)
	{
		switch (theNumParam)
		{
		case 0:
			mvAr = theVectParam;
			break;
		case 1:
			mvMa = theVectParam;
			break;
		case 2:
			if (theVectParam.GetSize() > 0)
				mvFracD = theVectParam[0];
			else
				throw cError("cArfima::Set - Size of theVectParam must be > 0");
			break;
		default:
			throw cError("cArfima::Set - theNumParam must be in 0, 1, 2");
			break;
		}
	}
	
	double  cArfima::Get(const uint theIndex, const uint theNumParam)
	{
		switch (theNumParam)
		{
		case 0:
			return mvAr[theIndex];
			break;
		case 1:
			return mvMa[theIndex];
			break;
		case 2:
			return mvFracD;
			break;
		}

	}

	cDVector& cArfima::Get(const uint theNumParam)
	{
	cDVector* myAux;
		switch (theNumParam)
		{
		case 0:
			return mvAr;
			break;
		case 1:
			return mvMa;
			break;
		case 2:
			myAux = new cDVector(1, mvFracD);
			return *myAux;
			break;
		}

	}

	cArfima& cArfima::operator =(const cArfima& theSrc)
	{
		if (GetCondMeanType() == eUnknown)
		{
			mvAr = cDVector(theSrc.mvAr);
			mvMa = cDVector(theSrc.mvMa);
			mvPolMa = cPolynome(theSrc.mvPolMa);
			SetCondMeanType(eArfima);
		}
		else
		{
			mvAr = theSrc.mvAr;
			mvMa = theSrc.mvMa ;
			mvPolMa = theSrc.mvPolMa;
		}
		mvFracD = theSrc.mvFracD;
		mvNTruncLag = theSrc.mvNTruncLag;
		
	uint myNParam = mvAr.GetSize() + mvMa.GetSize() + 1;
		mvGradPolMa = new cPolynome[myNParam];
		mvHessPolMa = new cPolynome*[myNParam];
		for (uint i = 0; i < myNParam; i++)
		{
			mvGradPolMa[i].Resize(mvNTruncLag);
			mvHessPolMa[i] = new cPolynome[myNParam];
			for (uint j = 0; j < myNParam; j++)
				mvHessPolMa[i][j].Resize(mvNTruncLag);
		}
		return *this;
	}
	
/*
	void  cArfima::UpdateProxyMeanParameters(void)
	{
	uint myNAr = GetNAr();
	uint myNMa = GetNMa();
	uint myNParam = myNAr + myNMa + 1;
	
	cPolynome myPhi(myNAr);
	cPolynome myTeta(myNMa);
		myPhi[0] = myTeta[0] = 1.0;
		for (uint i = 1; i <= myNAr; i++)
			myPhi[i] = -mvAr[i-1];
		for (uint i = 1; i <= myNMa; i++)
			myTeta[i] = mvMa[i-1];
	cPolynome myDelta(0);
		ComputeDeltaPowD(-mvFracD, mvNTruncLag, myDelta);
	cPolynome myPol1(0), myRest(0);
		IncrPowDiv(myTeta, myPhi, mvNTruncLag+myNAr, myPol1, myRest);
	cPolynome myPol2(0);
		IncrPowDiv(myPol1, myPhi, mvNTruncLag + myNAr, myPol2, myRest);
	cPolynome myPol3(0);
			myPol3 = TrunkMult(myDelta, myPol2, mvNTruncLag);
		mvPolMa = TrunkMult(myDelta, myPol1, mvNTruncLag);
		mvPolMa[0] = 0.0;
	cPolynome myPol4(0);
		IncrPowDiv(myDelta, myPhi, mvNTruncLag + myNAr, myPol4, myRest);
	cPolynome myPol5(mvNTruncLag);
		for (uint i = 1; i <= mvNTruncLag; i++)
			myPol5[i] = 1 / (double)i;
	cPolynome myPol6(0);
		myPol6 = TrunkMult(myPol5, myDelta, mvNTruncLag);
	cPolynome myPol7(0);
		myPol7 = TrunkMult(myPol6, myPol1, mvNTruncLag);
	cPolynome myPolX(1);
		myPolX[1] = 1;
		for (uint i = 0; i < myNAr; i++)
		{
			myPol3 = TrunkMult(myPol3, myPolX, mvNTruncLag);
			mvGradPolMa[i] = myPol3;
		}
		for (uint i = myNAr; i < myNAr + myNMa; i++)
		{
			myPol4 = TrunkMult(myPol4, myPolX, mvNTruncLag);
			mvGradPolMa[i] = myPol4;
		}
		mvGradPolMa[myNAr + myNMa] = myPol7;
	cDVector* myDPhi = new cDVector[myNAr];
	cDVector* myDTeta = new cDVector[myNMa];
	cDVector myDFracD(myNParam, 0.0);
		myDFracD[myNParam - 1] = 1.0;
		for (uint i = 0; i < myNAr ; i++)
		{
			myDPhi[i].ReAlloc(myNParam, 0.0);
			myDPhi[i] = 1.0;
		}
		for (uint j = 0; j < myNMa ; j++)
		{
			myDPhi[j].ReAlloc(myNParam, 0.0);
			myDPhi[j+myNAr] = 1.0;
		}
	
	}	
*/
	void  cArfima::ReAllocProxyMeanParameters(uint theOldNParam)
	{
		for (uint i = 0; i < theOldNParam; i++)
		{
			mvGradPolMa[i].Delete();
			for (uint j = 0; j < theOldNParam; j++)
				mvHessPolMa[i][j].Delete();
			delete[] mvHessPolMa[i];
		}
		delete[] mvGradPolMa;
		delete[] mvHessPolMa;
	uint myNParam = GetNParam();

		mvPolMa.Resize(mvNTruncLag);
		mvGradPolMa = new cPolynome[myNParam];
		mvHessPolMa = new cPolynome*[myNParam];
		for (uint i = 0; i < myNParam; i++)
		{
			mvGradPolMa[i].Resize(mvNTruncLag);
			mvHessPolMa[i] = new cPolynome[myNParam];
			for (uint j = 0; j < myNParam; j++)
				mvHessPolMa[i][j].Resize(mvNTruncLag);
		}

	}

	void  cArfima::UpdateProxyMeanParameters(void)
	{
		uint myNAr = GetNAr();
		uint myNMa = GetNMa();
		uint myNParam = myNAr + myNMa + 1;
		uint myNTrunk = mvNTruncLag ;
		cPolynome myP(myNAr);
		cPolynome myQ(myNMa);
		myP[0] = myQ[0] = 1.0;
		for (uint i = 1; i <= myNAr; i++)
			myP[i] = -mvAr[i - 1];
		for (uint i = 1; i <= myNMa; i++)
			myQ[i] = mvMa[i - 1];
		cPolynome myDeltaPmD(0);
		ComputeDeltaPowD(-mvFracD, myNTrunk, myDeltaPmD);
		cPolynome myLogDelta(0);
		ComputeLogDelta(myNTrunk, myLogDelta);
		cPolynome myUn(0);
		myUn[0] = 1.0;
		cPolynome myReste(0);
		cPolynome myPm1(0);
		IncrPowDiv(myUn, myP, myNTrunk, myPm1, myReste);
		cPolynome myPm1Q(0);
		myPm1Q = TrunkMult(myPm1, myQ, myNTrunk);
		cPolynome myPm1DeltaPmD(0);
		myPm1DeltaPmD = TrunkMult(myPm1, myDeltaPmD, myNTrunk);
		cPolynome myPm1QDeltaPmD(0);
		myPm1QDeltaPmD = TrunkMult(myPm1Q, myDeltaPmD, myNTrunk);		
		mvPolMa = myPm1QDeltaPmD;
		mvPolMa[0] = 0.0;

		cPolynome myPm2QDeltaPmD(0);
		myPm2QDeltaPmD = TrunkMult(myPm1, myPm1QDeltaPmD, myNTrunk);
		cPolynome myPm1QDeltaPmDLogDelta(0);
		myPm1QDeltaPmDLogDelta = TrunkMult(myLogDelta, myPm1QDeltaPmD, myNTrunk);
		cPolynome myX(1);
		myX[0] = 0.0;
		myX[1] = 1.0;
		for (uint i = 0; i < myNParam; i++)
			mvGradPolMa[i] = cPolynome(0);				
		// Derivée AR
		if (myNAr > 0)
		{
			mvGradPolMa[0] = TrunkMult(myX, myPm2QDeltaPmD, mvNTruncLag);
			for (uint i = 1; i < myNAr; i++)
				mvGradPolMa[i] = TrunkMult(myX, mvGradPolMa[i - 1], mvNTruncLag);
		}
		// Derivée MA
		if (myNMa > 0)
		{
			mvGradPolMa[myNAr] = TrunkMult(myX, myPm1DeltaPmD, mvNTruncLag);
			for (uint j = 1; j < myNMa; j++)
				mvGradPolMa[j+myNAr] = TrunkMult(myX, mvGradPolMa[j + myNAr - 1], mvNTruncLag);
		}
		// Derivee FracD
		mvGradPolMa[myNParam - 1] = TrunkPoly(myPm1QDeltaPmDLogDelta, mvNTruncLag);
		mvGradPolMa[myNParam - 1] *= -1.0;

		for (uint i = 0; i < myNParam; i++)
			(mvGradPolMa[i])[0] = 0.0;

		cPolynome myPm2DeltaPmD(0);
		myPm2DeltaPmD = TrunkMult(myPm1, myPm1DeltaPmD, myNTrunk);
		cPolynome myPm1DeltaPmdDLogDelta(0);
		myPm1DeltaPmdDLogDelta = TrunkMult(myLogDelta, myPm1DeltaPmD, myNTrunk);
		cPolynome myPm2QDeltaPmDLogDelta(0);
		myPm2QDeltaPmDLogDelta = TrunkMult(myLogDelta, myPm2QDeltaPmD, myNTrunk);
		cPolynome myPm3QDeltaPmD(0);
		myPm3QDeltaPmD = TrunkMult(myPm1, myPm2QDeltaPmD, myNTrunk);
		cPolynome myPm1QDeltaPmDLogDeltaP2(0);
		myPm1QDeltaPmDLogDeltaP2 = TrunkMult(myLogDelta, myPm1QDeltaPmDLogDelta, myNTrunk);

		for (uint i = 0; i < myNParam; i++)
			for (uint j = 0; j < myNParam; j++)
				mvHessPolMa[i][j] = cPolynome(0);

		// Hessien AR x AR		
		for (uint i = 1; i <= myNAr; i++)
		{
			for (uint j = i ; j <= myNAr; j++)
			{
				cPolynome myXpipj = cPolynome(i+j);
				myXpipj[i+j] = 2.0;
				mvHessPolMa[i - 1][j - 1] = TrunkMult(myPm3QDeltaPmD, myXpipj, mvNTruncLag);
				mvHessPolMa[j - 1][i - 1] = mvHessPolMa[i - 1][j - 1] ;
			}
		}
		// Hessien AR x MA
		for (uint i = 1; i <= myNAr; i++)
		{
			for (uint j = i ; j <= myNMa; j++)
			{
				cPolynome myXpipj = cPolynome(i+j);
				myXpipj[i+j] = 1.0;
				mvHessPolMa[i - 1][j - 1 + myNAr] = TrunkMult(myPm2DeltaPmD, myXpipj, mvNTruncLag);
				mvHessPolMa[j - 1 + myNAr][i - 1] = mvHessPolMa[i - 1][j - 1 + myNAr];
			}
		}
		// Hessien AR x FracD
		for (uint i = 1; i <= myNAr; i++)
		{
			cPolynome myXpi = cPolynome(i);
			myXpi[i] = -1.0;
			mvHessPolMa[i - 1][myNAr + myNMa] = TrunkMult(myPm2QDeltaPmDLogDelta, myXpi, mvNTruncLag);
			mvHessPolMa[myNAr + myNMa][i - 1] = mvHessPolMa[i - 1][myNAr + myNMa];
		}
		// Hessien MA x MA - Y'en n'a pas

		// Hessien MA x FracD
		for (uint i = 1; i <= myNMa; i++)
		{
			cPolynome myXpi = cPolynome(i);
			myXpi[i] = -1.0;
			mvHessPolMa[i + myNAr - 1][myNAr + myNMa] =  TrunkMult(myPm1DeltaPmdDLogDelta, myXpi, mvNTruncLag);
			mvHessPolMa[myNAr + myNMa][i + myNAr - 1] = mvHessPolMa[i + myNAr - 1][myNAr + myNMa];
		}
		//Hessien FracD x FracD
		mvHessPolMa[myNParam - 1][myNParam - 1] = TrunkPoly(myPm1QDeltaPmDLogDeltaP2, mvNTruncLag);

		for (uint i = 0; i < myNParam; i++)
			for (uint j = 0; j < myNParam; j++)
				(mvHessPolMa[i][j])[0] = 0.0;
//		delete[] myMult;

	}

	double cArfima::ComputeMean(uint theDate, const cRegArchValue& theData) const
	{
	
		return mvPolMa.BackwardPolOp(theData.mUt, theDate);

	}
	
	uint cArfima::GetNParam(void) const
	{
		return mvAr.GetSize() + mvMa.GetSize() + 1;
	}
	
	uint cArfima::GetNAr(void) const
	{
		return mvAr.GetSize();
	}

	uint cArfima::GetNMa(void) const
	{
		return mvMa.GetSize();
	}
		
	uint cArfima::GetNLags(void) const
	{
		return MAX(MAX(mvMa.GetSize(), mvAr.GetSize()), mvNTruncLag);
	}

	void cArfima::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, uint theBegIndex)
	{
	uint myNAr = GetNAr();
	uint myNMa = GetNMa();
	uint myNParam = myNAr + myNMa + 1;
		
	cDVector myRes(theGradData.GetNParam());
		for (uint i = 0; i < myNParam; i++)
			myRes[theBegIndex + i] = mvGradPolMa[i].BackwardPolOp(theValue.mUt, theDate);
		for (uint t = 1; t < MIN(theDate, mvNTruncLag + 1); t++)
			myRes -= mvPolMa[t] * theGradData.mGradMt[t - 1];

		theGradData.mCurrentGradMu += myRes;
	}

	void cArfima::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{
	uint mySize = GetNParam();
		if (theDestVect.GetSize() < mySize + theIndex)
				throw cError("Wrong size");
		mvAr.SetSubVectorWithThis(theDestVect, theIndex);
		mvMa.SetSubVectorWithThis(theDestVect, theIndex + mvAr.GetSize());
		theDestVect[theIndex + mvAr.GetSize() + mvMa.GetSize()] = mvFracD;

	}
	
	void cArfima::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	uint mySize = theSrcVect.GetSize();
		if (GetNParam() + theIndex > mySize)
			throw cError("Wrong size");
		mvAr.SetThisWithSubVector(theSrcVect, theIndex);
		mvMa.SetThisWithSubVector(theSrcVect, theIndex + mvAr.GetSize());
		mvFracD = theSrcVect[theIndex + mvAr.GetSize() + mvMa.GetSize()];


	}

	void cArfima::ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex)
	{
	uint myNAr = GetNAr();
	uint myNMa = GetNMa();
	uint myNParam = myNAr + myNMa + 1;
	cDMatrix myRes = cDMatrix(theGradData.mCurrentGradMu.GetSize(), theGradData.mCurrentGradMu.GetSize());

/*
		for (uint i = 0; i < myNParam; i++)
			for (uint j = i; j < myNParam; j++)
				theHessData.mCurrentHessVar[myBegIndex + i][myBegIndex + j] = mvHessPolAr[i][j].BackwardPolOp(theData.mUt, theDate, 2);
		for (uint k = 0; k < myNParam; k++)
			for (uint l = k+1; l < myNParam; l++)
				theHessData.mCurrentHessVar[myBegIndex + l][myBegIndex + k] = theHessData.mCurrentHessVar[myBegIndex + k][myBegIndex + l];

		cDVector myGradPol(theGradData.GetNParam());
		for (uint n = 1; n < MIN(theDate, mvNTruncLag + 1); n++)
		{
			for (uint i = 0; i < myNParam; i++)
				myGradPol[i+myBegIndex] = mvGradPolAr[i].mCoeff[n];
			theHessData.mCurrentHessVar -= 2 * theData.mUt[theDate - n] * (myGradPol*Transpose(theGradData.mGradMt[n - 1]) + theGradData.mGradMt[n - 1] * Transpose(myGradPol));
			theHessData.mCurrentHessVar += 2 * mvPolAr[n] * theGradData.mGradMt[n - 1] * Transpose(theGradData.mGradMt[n - 1]);
			theHessData.mCurrentHessVar -= 2 * mvPolAr[n] * theData.mUt[theDate - n] * theHessData.mHessMt[n-1];
		}

*/

		for (uint i = 0; i < myNParam; i++)
			for (uint j = i; j < myNParam; j++)
				myRes[i + theBegIndex][j + theBegIndex] += mvHessPolMa[i][j].BackwardPolOp(theData.mUt, theDate);
		for (uint k = 0; k < myNParam; k++)
			for (uint l = k + 1; l < myNParam; l++)
				myRes[theBegIndex + l][theBegIndex + k] = myRes[theBegIndex + k][theBegIndex + l];

		cDVector myGradPol(theGradData.GetNParam());
		for (uint n = 1; n < MIN(theDate, mvNTruncLag + 1); n++)
		{
			for (uint i = 0; i < myNParam; i++)
				myGradPol[i + theBegIndex] = mvGradPolMa[i].mCoeff[n];
			myRes -= myGradPol * Transpose(theGradData.mGradMt[n - 1]) + theGradData.mGradMt[n - 1] * Transpose(myGradPol);
			myRes -= mvPolMa[n] * theHessData.mHessMt[n - 1];
		}
		theHessData.mCurrentHessMu += myRes;
	}
	
	void cArfima::ComputeGradAndHess(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex)
	{

		uint myNAr = GetNAr();
		uint myNMa = GetNMa();
		uint myNParam = myNAr + myNMa + 1;

		cDVector myRes(theGradData.GetNParam());
		for (uint i = 0; i < myNParam; i++)
			myRes[theBegIndex + i] = mvGradPolMa[i].BackwardPolOp(theValue.mUt, theDate);
		for (int t = 1; t < MIN(theDate, mvNTruncLag + 1); t++)
			myRes -= mvPolMa[t] * theGradData.mGradMt[t - 1];
		theGradData.mCurrentGradMu += myRes;

		cDMatrix myMat = cDMatrix(theGradData.mCurrentGradMu.GetSize(), theGradData.mCurrentGradMu.GetSize());

		for (uint i = 0; i < myNParam; i++)
			for (uint j = i; j < myNParam; j++)
				myMat[i + theBegIndex][j + theBegIndex] += mvHessPolMa[i][j].BackwardPolOp(theValue.mUt, theDate);
		for (uint k = 0; k < myNParam; k++)
			for (uint l = k + 1; l < myNParam; l++)
				myMat[theBegIndex + l][theBegIndex + k] = myMat[theBegIndex + k][theBegIndex + l];

		cDVector myGradPol(theGradData.GetNParam());
		for (uint n = 1; n < MIN(theDate, mvNTruncLag + 1); n++)
		{
			for (uint i = 0; i < myNParam; i++)
				myGradPol[i + theBegIndex] = mvGradPolMa[i].mCoeff[n];
			myMat -= myGradPol * Transpose(theGradData.mGradMt[n - 1]) + theGradData.mGradMt[n - 1] * Transpose(myGradPol);
			myMat -= mvPolMa[n] * theHessData.mHessMt[n - 1];
		}
		theHessData.mCurrentHessMu += myMat;
	}
	
	void cArfima::GetParamName(uint theIndex, char** theName)
	{
	uint myIndex = theIndex;
		for (uint i = 0; i < mvAr.GetSize(); i++)
		{
			sprintf(theName[myIndex++], "AR[%d]", i + 1);

		}
		for (uint i = 0; i < mvMa.GetSize(); i++)
		{
			sprintf(theName[myIndex++], "MA[%d]", i + 1);

		}
		sprintf(theName[myIndex], "Frac. d");

	}

	void cArfima::GetParamName(uint theIndex, string theName[])
	{
	uint myIndex = theIndex;
	char myChar[100] ;
		for (uint i = 0; i < mvAr.GetSize(); i++)
		{
			sprintf(myChar, "AR[%d]", i + 1);
			theName[myIndex++] = myChar;

		}
		for (uint i = 0; i < mvMa.GetSize(); i++)
		{
			sprintf(myChar, "MA[%d]", i + 1);
			theName[myIndex++] = myChar;

		}
		sprintf(myChar, "Frac. d");
		theName[myIndex] = myChar;

	}

} // namespace
