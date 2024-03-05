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
		mvPolAr.Resize(mvNTruncLag);
	uint myNParam = theNAr + theNMa + 1;
		
		mvGradPolAr = new cPolynome[myNParam];
		mvHessPolAr = new cPolynome*[myNParam];
		for (uint i = 0; i < myNParam; i++)
		{
			mvGradPolAr[i].Resize(mvNTruncLag);
			mvHessPolAr[i] = new cPolynome[myNParam];
			for (uint j = 0; j < myNParam; j++)
				mvHessPolAr[i][j].Resize(mvNTruncLag);
		}
		MESS_CREAT("cArfima")
	}

	cArfima::cArfima(const cDVector& theAr, const cDVector& theMa, double theFracD, uint theNTruncLag) :cAbstCondMean(eArfima)
	{
		mvAr = theAr;
		mvMa = theMa;
		mvFracD = theFracD;
		mvNTruncLag = theNTruncLag;
		mvPolAr.Resize(mvNTruncLag);
		uint myNParam = theAr.GetSize() + theMa.GetSize() + 1;

		mvGradPolAr = new cPolynome[myNParam];
		mvHessPolAr = new cPolynome*[myNParam];
		for (uint i = 0; i < myNParam; i++)
		{
			mvGradPolAr[i].Resize(mvNTruncLag);
			mvHessPolAr[i] = new cPolynome[myNParam];
			for (uint j = 0; j < myNParam; j++)
				mvHessPolAr[i][j].Resize(mvNTruncLag);
		}
		MESS_CREAT("cArfima")
	}

	cArfima::cArfima(const cArfima& theArfima) :cAbstCondMean(eUnknown)
	{	
#ifdef _CLI_
		uint myNAr = theArfima.mvAr.GetSize();
		uint myNMa = theArfima.mvMa.GetSize();

//		cout << "Debut de cArfima(const cArfima& theArfima)" << endl;

		mvAr.ReAlloc(myNAr);
		mvMa.ReAlloc(myNMa);

//		cout << "ly 1" << endl;


		mvNTruncLag = theArfima.mvNTruncLag;;
		uint myNParam = myNAr + myNMa + 1;
		mvPolAr.Resize(mvNTruncLag);

//		cout << "ly 2" << endl;

		mvGradPolAr = new cPolynome[myNParam];
	 	mvHessPolAr = new cPolynome*[myNParam];

//		cout << "ly 3" << endl;

		for (uint i = 0; i < myNParam; i++)
		{
			mvGradPolAr[i].Resize(mvNTruncLag);
			mvHessPolAr[i] = new cPolynome[myNParam];
			for (uint j = 0; j < myNParam; j++)
				mvHessPolAr[i][j].Resize(mvNTruncLag);
		}
#endif _CLI_
		*this = theArfima;
//		cout << "ly 4" << endl;
		MESS_CREAT("cArfima")
	}

	cArfima:: ~cArfima()
	{
	uint myNParam = mvAr.GetSize() + mvMa.GetSize() + 1;
		mvAr.Delete();
		mvMa.Delete();
		mvPolAr.Delete();
		for (uint i = 0; i < myNParam; i++)
		{
			mvGradPolAr[i].Delete();
			for (uint j = 0; j < myNParam; j++)
				mvHessPolAr[i][j].Delete();
			delete[] mvHessPolAr[i];
		}
		delete[] mvGradPolAr;
		delete[] mvHessPolAr;

			MESS_DESTR("cArfima")
	}

	void cArfima::DeletePoly(void)
	{
		uint myNParam = mvAr.GetSize() + mvMa.GetSize() + 1;
		mvPolAr.Delete();
		for (uint i = 0; i < myNParam; i++)
		{
			mvGradPolAr[i].Delete();
			for (uint j = 0; j < myNParam; j++)
				mvHessPolAr[i][j].Delete();
			delete[] mvHessPolAr[i];
		}
		delete[] mvGradPolAr;
		delete[] mvHessPolAr;
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
		mvPolAr.Delete();
		for (uint i = 0; i < myNParam; i++)
		{
			mvGradPolAr[i].Delete();
			for (uint j = 0; j < myNParam; j++)
				mvHessPolAr[i][j].Delete();
			delete[] mvHessPolAr[i];
		}
		delete[] mvGradPolAr;
		delete[] mvHessPolAr;
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
		Rprintf("ARFIMA(%d, d, %d) model with:\n", myNAr, myNMa);
		for (uint i = 0; i < myNAr; i++)
			Rprintf("\tAR[%d]=%f\n", i + 1, mvAr[i]);
		for (uint i = 0; i < myNMa; i++)
			Rprintf("\tMA[%d]=%f\n", i + 1, mvMa[i]);
		Rprintf("\td=%f\n", mvFracD);
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

	void cArfima::ReSizePoly(const uint theNParam)
	{
		DeletePoly();
		mvPolAr.Resize(mvNTruncLag);
		mvGradPolAr = new cPolynome[theNParam];
		mvHessPolAr = new cPolynome* [theNParam];
		for (uint i = 0; i < theNParam; i++)
		{
			mvGradPolAr[i].Resize(mvNTruncLag);
			mvHessPolAr[i] = new cPolynome[theNParam];
			for (uint j = 0; j < theNParam; j++)
				mvHessPolAr[i][j].Resize(mvNTruncLag);
		}
	}

	void  cArfima::ReAlloc(const uint theSize, const uint theNumParam)
	{
		cDVector myAr, myMa;
		uint myNParam = 0;
		switch (theNumParam)
		{
			case 0:
				myMa = mvMa;
				myAr.ReAlloc(theSize);
				myNParam = myMa.GetSize() + myAr.GetSize() + 1;
				ReSizePoly(myNParam);
				mvAr = myAr;
				break;
			case 1:
				myAr = mvAr;
				myMa.ReAlloc(theSize);
				myNParam = myMa.GetSize() + myAr.GetSize() + 1;
				ReSizePoly(myNParam);
				mvMa = myMa;
				break;
			default:
				throw cError("cArfima::ReAlloc - theNumParam must be in 0 or 1");
				break;
		}
	}
	
	void  cArfima::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{
		uint myNParam = 0;
		switch (theNumParam)
		{
		case 0: // mvAr
			myNParam = theVectParam.GetSize() + mvMa.GetSize() + 1;
			ReSizePoly(myNParam);
			mvAr = theVectParam;
			break;
		case 1: // mvMa
			myNParam = theVectParam.GetSize() + mvAr.GetSize() + 1;
			ReSizePoly(myNParam);
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
		case 1: // MA
			if (theIndex < mvMa.GetSize())
				mvMa[theIndex] = theValue;
			else
				throw cError("cArfima::Set - wrong index");
			break;
		case 2: // FracD
			mvFracD = theValue;
			break;
		
		case 3: // mvNTruncLag
			mvNTruncLag = (uint)theValue;
			ReSizePoly(mvAr.GetSize() + mvMa.GetSize() + 1);
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
			ReSizePoly(theVectParam.GetSize() + mvMa.GetSize() + 1);
			mvAr = theVectParam;
			break;
		case 1:
			ReSizePoly(theVectParam.GetSize() + mvAr.GetSize() + 1);
			mvMa = theVectParam;
			break;
		case 2:
			if (theVectParam.GetSize() > 0)
				mvFracD = theVectParam[0];
			else
				throw cError("cArfima::Set - Size of theVectParam must be > 0");
			break;
		case 3:
			if (theVectParam.GetSize() > 0)
			{
				mvNTruncLag = theVectParam[0];
				ReSizePoly(mvAr.GetSize() + mvMa.GetSize() + 1);
			}
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
		case 3:
			return mvNTruncLag;
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
		case 3:
			myAux = new cDVector(1, mvNTruncLag);
			return *myAux;
			break;

		}

	}

	cArfima& cArfima::operator =(const cArfima& theSrc)
	{
		mvAr = theSrc.mvAr;
		mvMa = theSrc.mvMa ;
		mvPolAr = theSrc.mvPolAr;
		mvFracD = theSrc.mvFracD;
		mvNTruncLag = theSrc.mvNTruncLag;
		DeletePoly();
		SetCondMeanType(eArfima);
	uint myNParam = mvAr.GetSize() + mvMa.GetSize() + 1;
		mvPolAr = theSrc.mvPolAr;
		mvGradPolAr = new cPolynome[myNParam];
		mvHessPolAr = new cPolynome*[myNParam];
		for (uint i = 0; i < myNParam; i++)
		{
			mvHessPolAr[i] = new cPolynome[myNParam];
			mvGradPolAr[i] = theSrc.mvGradPolAr[i];
			for (uint j = 0; j < myNParam; j++)
			{
				mvHessPolAr[i][j] = theSrc.mvHessPolAr[i][j];
			}
		}
		return *this;
	}

	void  cArfima::ReAllocProxyMeanParameters(uint theOldNParam)
	{
		for (uint i = 0; i < theOldNParam; i++)
		{
			mvGradPolAr[i].Delete();
			for (uint j = 0; j < theOldNParam; j++)
				mvHessPolAr[i][j].Delete();
			delete[] mvHessPolAr[i];
		}
		delete[] mvGradPolAr;
		delete[] mvHessPolAr;
	uint myNParam = GetNParam();

		mvPolAr.Resize(mvNTruncLag);
		mvGradPolAr = new cPolynome[myNParam];
		mvHessPolAr = new cPolynome*[myNParam];
		for (uint i = 0; i < myNParam; i++)
		{
			mvGradPolAr[i].Resize(mvNTruncLag);
			mvHessPolAr[i] = new cPolynome[myNParam];
			for (uint j = 0; j < myNParam; j++)
				mvHessPolAr[i][j].Resize(mvNTruncLag);
		}

	}

	void  cArfima::UpdateProxyMeanParameters(void)
	{
		uint myNAr = GetNAr();
		uint myNMa = GetNMa();
		uint myNParam = myNAr + myNMa + 1;
		cPolynome myP(myNAr);
		cPolynome myQ(myNMa);
		myP[0] = myQ[0] = 1.0;
		for (uint i = 1; i <= myNAr; i++)
			myP[i] = -mvAr[i - 1];
		for (uint i = 1; i <= myNMa; i++)
			myQ[i] = mvMa[i - 1];
		cPolynome myDeltaPD(0);
		ComputeDeltaPowD(mvFracD, mvNTruncLag, myDeltaPD);
		cPolynome myLogDelta(0);
		ComputeLogDelta(mvNTruncLag, myLogDelta);
		cPolynome myUn(0);
		myUn[0] = 1.0;
		cPolynome myReste(0);
		cPolynome myQm1(0);
		IncrPowDiv(myUn, myQ, mvNTruncLag, myQm1, myReste);
		cPolynome myQm1P(0);
		myQm1P = TrunkMult(myQm1, myP, mvNTruncLag);
		cPolynome myQm1DeltaPD(0);
		myQm1DeltaPD = TrunkMult(myQm1, myDeltaPD, mvNTruncLag);
		cPolynome myQm1PDeltaPD(0);
		myQm1PDeltaPD = TrunkMult(myQm1P, myDeltaPD, mvNTruncLag);
		mvPolAr = myUn - myQm1PDeltaPD;
//		mvPolAr[0] = 0.0;

		cPolynome myQm2PDeltaPD(0);
		myQm2PDeltaPD = TrunkMult(myQm1, myQm1PDeltaPD, mvNTruncLag);
		cPolynome myQm1PDeltaPDLogDelta(0);
		myQm1PDeltaPDLogDelta = TrunkMult(myLogDelta, myQm1PDeltaPD, mvNTruncLag);
		for (uint i = 0; i < myNParam; i++)
			mvGradPolAr[i] = cPolynome(0);				
		// Derivée AR
		if (myNAr > 0)
		{
			for (uint i = 1; i <= myNAr; i++)
			{
				cPolynome myXpi(i);
				myXpi[i] = 1.0;
				mvGradPolAr[i-1] = TrunkMult(myXpi, myQm1DeltaPD, mvNTruncLag);
			}
		}
		// Derivée MA
		if (myNMa > 0)
		{
			for (uint j = 1; j <= myNMa; j++)
			{
				cPolynome myXpj(j);
				myXpj[j] = 1.0;
				mvGradPolAr[j - 1 + myNAr] = TrunkMult(myXpj, myQm2PDeltaPD, mvNTruncLag);
			}
		}
		// Derivee FracD
		mvGradPolAr[myNParam - 1] = -1.0 * myQm1PDeltaPDLogDelta;

//		for (uint i = 0; i < myNParam; i++)
//			(mvGradPolAr[i])[0] = 0.0;

		for (uint i = 0; i < myNParam; i++)
			for (uint j = 0; j < myNParam; j++)
				mvHessPolAr[i][j] = cPolynome(0);

		// Hessien AR x AR -> 0		

		// Hessien AR x MA
		cPolynome myQm2DeltaPD(0);
		myQm2DeltaPD = TrunkMult(myQm1, myQm1DeltaPD, mvNTruncLag);
		for (uint i = 1; i <= myNAr; i++)
		{
			for (uint j = i ; j <= myNMa; j++)
			{
				cPolynome myXpipj = cPolynome(i+j);
				myXpipj[i+j] = -1.0;
				mvHessPolAr[i - 1][j - 1 + myNAr] = TrunkMult(myQm2DeltaPD, myXpipj, mvNTruncLag);
				mvHessPolAr[j - 1 + myNAr][i - 1] = mvHessPolAr[i - 1][j - 1 + myNAr];
			}
		}
		// Hessien AR x FracD
		cPolynome myQPm1DeltaPDLogDelta(0);
		myQPm1DeltaPDLogDelta = TrunkMult(myLogDelta, myQm1DeltaPD, mvNTruncLag);
		for (uint i = 1; i <= myNAr; i++)
		{
			cPolynome myXpi = cPolynome(i);
			myXpi[i] = 1.0;
			mvHessPolAr[i - 1][myNAr + myNMa] = TrunkMult(myQPm1DeltaPDLogDelta, myXpi, mvNTruncLag);
			mvHessPolAr[myNAr + myNMa][i - 1] = mvHessPolAr[i - 1][myNAr + myNMa];
		}
		// Hessien MA x MA 
		cPolynome myQm3PDeltaPD(0);
		myQm3PDeltaPD = TrunkMult(myQm1, myQm2PDeltaPD, mvNTruncLag);
		for (uint i = 1; i <= myNMa; i++)
		{
			for (uint j = i; j <= myNMa; j++)
			{
				cPolynome myXpipj = cPolynome(i+j);
				myXpipj[i+j] = -2.0;
				mvHessPolAr[i - 1 + myNAr][j - 1 + myNAr] = TrunkMult(myQm3PDeltaPD, myXpipj, mvNTruncLag);
				if (i != j)
					mvHessPolAr[j - 1 + myNAr][i - 1 + myNAr] = mvHessPolAr[i - 1 + myNAr][j - 1 + myNAr];
			}
		}

		// Hessien MA x FracD
		cPolynome myQm2PDeltaPDLogDelta(0);
		myQm2PDeltaPDLogDelta = TrunkMult(myLogDelta, myQm2PDeltaPD, mvNTruncLag);
		for (uint i = 1; i <= myNMa; i++)
		{
			cPolynome myXpi = cPolynome(i);
			myXpi[i] = 1.0;
			mvHessPolAr[i + myNAr - 1][myNAr + myNMa] =  TrunkMult(myQm2PDeltaPDLogDelta, myXpi, mvNTruncLag);
			mvHessPolAr[myNAr + myNMa][i + myNAr - 1] = mvHessPolAr[i + myNAr - 1][myNAr + myNMa];
		}
		//Hessien FracD x FracD
		cPolynome myQm1PDeltaPDLogDeltaP2(0);
		myQm1PDeltaPDLogDeltaP2 = TrunkMult(myLogDelta, myQm1PDeltaPDLogDelta, mvNTruncLag);

		mvHessPolAr[myNParam - 1][myNParam - 1] = -1.0* myQm1PDeltaPDLogDeltaP2;

//		for (uint i = 0; i < myNParam; i++)
//			for (uint j = 0; j < myNParam; j++)
//				(mvHessPolAr[i][j])[0] = 0.0;
	}

	double cArfima::ComputeMean(uint theDate, const cRegArchValue& theData) const
	{
	
		return mvPolAr.BackwardPolOp(theData.mYt, theDate);

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
			myRes[theBegIndex + i] = mvGradPolAr[i].BackwardPolOp(theValue.mYt, theDate);
//		for (uint t = 1; t < MIN(theDate, mvNTruncLag + 1); t++)
//			myRes -= mvPolAr[t] * theGradData.mGradMt[t - 1];

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
				myRes[i + theBegIndex][j + theBegIndex] += mvHessPolAr[i][j].BackwardPolOp(theData.mYt, theDate);
		for (uint k = 0; k < myNParam; k++)
			for (uint l = k + 1; l < myNParam; l++)
				myRes[theBegIndex + l][theBegIndex + k] = myRes[theBegIndex + k][theBegIndex + l];

//		cDVector myGradPol(theGradData.GetNParam());
//		for (uint n = 1; n < MIN(theDate, mvNTruncLag + 1); n++)
//		{
//			for (uint i = 0; i < myNParam; i++)
//				myGradPol[i + theBegIndex] = mvGradPolAr[i].mCoeff[n];
//			myRes -= myGradPol * Transpose(theGradData.mGradMt[n - 1]) + theGradData.mGradMt[n - 1] * Transpose(myGradPol);
//			myRes -= mvPolAr[n] * theHessData.mHessMt[n - 1];
//		}
		theHessData.mCurrentHessMu += myRes;
	}
	
	void cArfima::ComputeGradAndHess(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex)
	{

		uint myNAr = GetNAr();
		uint myNMa = GetNMa();
		uint myNParam = myNAr + myNMa + 1;

		cDVector myRes(theGradData.GetNParam());
		for (uint i = 0; i < myNParam; i++)
			myRes[theBegIndex + i] = mvGradPolAr[i].BackwardPolOp(theValue.mUt, theDate);
//		for (int t = 1; t < MIN(theDate, mvNTruncLag + 1); t++)
//			myRes -= mvPolAr[t] * theGradData.mGradMt[t - 1];
		theGradData.mCurrentGradMu += myRes;

		cDMatrix myMat = cDMatrix(theGradData.mCurrentGradMu.GetSize(), theGradData.mCurrentGradMu.GetSize());

		for (uint i = 0; i < myNParam; i++)
			for (uint j = i; j < myNParam; j++)
				myMat[i + theBegIndex][j + theBegIndex] += mvHessPolAr[i][j].BackwardPolOp(theValue.mYt, theDate);
		for (uint k = 0; k < myNParam; k++)
			for (uint l = k + 1; l < myNParam; l++)
				myMat[theBegIndex + l][theBegIndex + k] = myMat[theBegIndex + k][theBegIndex + l];

//		cDVector myGradPol(theGradData.GetNParam());
//		for (uint n = 1; n < MIN(theDate, mvNTruncLag + 1); n++)
//		{
//			for (uint i = 0; i < myNParam; i++)
//				myGradPol[i + theBegIndex] = mvGradPolAr[i].mCoeff[n];
//			myMat -= myGradPol * Transpose(theGradData.mGradMt[n - 1]) + theGradData.mGradMt[n - 1] * Transpose(myGradPol);
//			myMat -= mvPolAr[n] * theHessData.mHessMt[n - 1];
//		}
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
