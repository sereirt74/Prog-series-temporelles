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
		*this = theArfima;
		MESS_CREAT("cAr")
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
			mvPolAr = cPolynome(theSrc.mvPolAr);
			SetCondMeanType(eArfima);
		}
		else
		{
			mvAr = theSrc.mvAr;
			mvMa = theSrc.mvMa ;
			mvPolAr = theSrc.mvPolAr;
		}
		mvFracD = theSrc.mvFracD;
		mvNTruncLag = theSrc.mvNTruncLag;
		
	uint myNParam = mvAr.GetSize() + mvMa.GetSize() + 1;
		mvGradPolAr = new cPolynome[myNParam];
		mvHessPolAr = new cPolynome*[myNParam];
		for (uint i = 0; i < myNParam; i++)
		{
			mvGradPolAr[i].Resize(mvNTruncLag);
			mvHessPolAr[i] = new cPolynome[myNParam];
			for (uint j = 0; j < myNParam; j++)
				mvHessPolAr[i][j].Resize(mvNTruncLag);
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
		mvPolAr = TrunkMult(myDelta, myPol1, mvNTruncLag);
		mvPolAr[0] = 0.0;
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
			mvGradPolAr[i] = myPol3;
		}
		for (uint i = myNAr; i < myNAr + myNMa; i++)
		{
			myPol4 = TrunkMult(myPol4, myPolX, mvNTruncLag);
			mvGradPolAr[i] = myPol4;
		}
		mvGradPolAr[myNAr + myNMa] = myPol7;
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
		mvPolAr = -1*myQm1PDeltaPD;
		mvPolAr[0] = 0.0;
		myQm1PDeltaPD[0] = 0;
		cPolynome myQm2PDeltaPD(0);
		myQm2PDeltaPD = TrunkMult(myQm1, myQm1PDeltaPD, mvNTruncLag);
		cPolynome myQm1PDeltaPDLogDelta(0);
		myQm1PDeltaPDLogDelta = TrunkMult(myLogDelta, myQm1PDeltaPD, mvNTruncLag);
		cPolynome myX(1);
		myX[0] = 0.0;
		myX[1] = 1.0;
		for (uint i = 0; i < myNParam; i++)
			mvGradPolAr[i] = cPolynome(0);				
		// Derivée MA
		if (myNMa > 0)
		{
			mvGradPolAr[myNAr] = TrunkMult(myX, myQm2PDeltaPD, mvNTruncLag);
			for (uint i = 1; i < myNMa; i++)
				mvGradPolAr[i+myNAr] = TrunkMult(myX, mvGradPolAr[i + myNAr - 1], mvNTruncLag);
		}
		// Derivée AR
		if (myNAr > 0)
		{
			mvGradPolAr[0] = TrunkMult(myX, myQm1DeltaPD, mvNTruncLag);
			for (uint j = 1; j < myNAr; j++)
				mvGradPolAr[j] = TrunkMult(myX, mvGradPolAr[j - 1], mvNTruncLag);
		}
		// Derivee FracD
		mvGradPolAr[myNParam - 1] = myQm1PDeltaPDLogDelta;
		mvGradPolAr[myNParam - 1] *= -1.0;

		for (uint i = 0; i < myNParam; i++)
			(mvGradPolAr[i])[0] = 0.0;

		cPolynome myQm2DeltaPD(0);
		myQm2DeltaPD = TrunkMult(myQm1, myQm1DeltaPD, mvNTruncLag);
		cPolynome myQm1DeltaPDLogDelta(0);
		myQm1DeltaPDLogDelta = TrunkMult(myLogDelta, myQm1DeltaPD, mvNTruncLag);
		cPolynome myQm2PDeltaPDLogDelta(0);
		myQm2PDeltaPDLogDelta = TrunkMult(myLogDelta, myQm2PDeltaPD, mvNTruncLag);
		cPolynome myQm3PDeltaPD(0);
		myQm3PDeltaPD = TrunkMult(myQm1, myQm2PDeltaPD, mvNTruncLag);
		cPolynome myQm1PDeltaPDLogDeltaP2(0);
		myQm1PDeltaPDLogDeltaP2 = TrunkMult(myLogDelta, myQm1PDeltaPDLogDelta, mvNTruncLag);

		for (uint i = 0; i < myNParam; i++)
			for (uint j = 0; j < myNParam; j++)
				mvHessPolAr[i][j] = cPolynome(0);

		// Hessien MA x MA		
		if (myNMa > 0)
		{
			for (uint i = 1; i <= myNMa; i++)
			{
				for (uint j = i ; j <= myNMa; j++)
				{
					cPolynome myXpipj = cPolynome(i+j);
					myXpipj[i+j] = -2.0;

					mvHessPolAr[i + myNAr - 1][j + myNAr - 1] = TrunkMult(myQm3PDeltaPD, myXpipj, mvNTruncLag);
					mvHessPolAr[j + myNAr - 1][i + myNAr - 1] = mvHessPolAr[i + myNAr - 1][j + myNAr - 1];
				}
			}
			// Hessien MA x AR
			if (myNAr > 0)
			{
				for (uint i = 1; i <= myNMa; i++)
				{
					for (uint j = i; j <= myNAr; j++)
					{
						cPolynome myXpipj = cPolynome(i + j);
						myXpipj[i + j] = -1.0;
						mvHessPolAr[i + myNAr - 1][j - 1] = TrunkMult(myQm2DeltaPD, myXpipj, mvNTruncLag);
						mvHessPolAr[j - 1][i + myNAr - 1] = mvHessPolAr[i + myNAr - 1][j - 1];
					}
				}
			}
			// Hessien MA x FracD
			for (uint i = 1; i <= myNMa; i++)
			{
				cPolynome myXpi = cPolynome(i);
				myXpi[i] = -1.0;
				mvHessPolAr[myNAr + i - 1][myNAr + myNMa] = TrunkMult(myQm1PDeltaPDLogDelta, myXpi, mvNTruncLag);
				mvHessPolAr[myNAr + myNMa][myNAr + i - 1] = mvHessPolAr[myNAr + i - 1][myNAr + myNMa];
			}
		}
		// Hessien AR x AR - Y'en n'a pas

		// Hessien AR x FracD
		if (myNAr > 0)
		{
			for (uint i = 1; i <= myNAr; i++)
			{
				cPolynome myXpi = cPolynome(i);
				myXpi[i] = 1.0;
				mvHessPolAr[i - 1][myNAr + myNMa] = TrunkMult(myQm1DeltaPDLogDelta, myXpi, mvNTruncLag);
				mvHessPolAr[myNAr + myNMa][i - 1] = mvHessPolAr[i - 1][myNAr + myNMa];
			}
		}
		//Hessien FracD x FracD
		mvHessPolAr[myNParam - 1][myNParam - 1] = myQm1PDeltaPDLogDeltaP2*(- 1.0);

		for (uint i = 0; i < myNParam; i++)
			for (uint j = 0; j < myNParam; j++)
				(mvHessPolAr[i][j])[0] = 0.0;
//		delete[] myMult;

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

	uint cArfima::GetNu(void) const
	{
		return 0;
	}

	uint cArfima::GetNh(void) const
	{
		return 0;
	}

	void cArfima::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, uint theBegIndex, cAbstResiduals* theResiduals)
	{
#ifndef _DEBUG1
	uint myNAr = GetNAr();
	uint myNMa = GetNMa();
	uint myNParam = myNAr + myNMa + 1;
		for (uint i = 0; i < myNParam; i++)
			theGradData.mCurrentGradMu[theBegIndex + i] += mvGradPolAr[i].BackwardPolOp(theValue.mYt, theDate);
#else
		NumericComputeGrad(theDate, theValue, theGradData, theBegIndex, theResiduals);
#endif //_DEBUG1

	}

	void cArfima::GetNParamF(uint theNParam[3]) const
	{
	uint myp = mvAr.GetSize();
	uint myq = mvMa.GetSize();

		theNParam[0] = myp + myq + 1;
		theNParam[1] = myq;
		theNParam[2] = 0;
	}

	void cArfima::ComputeGradF(uint theDate, const cRegArchValue& theValue, cDVector& theGradF)
	{
	}

	void cArfima::NumericComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, uint theBegIndex, cAbstResiduals* theResiduals, double theh)
	{
	uint myNAr = GetNAr();
	uint myNMa = GetNMa();
		UpdateProxyMeanParameters();
	double myF0 = ComputeMean(theDate, theValue);
		for (uint i = 0; i < myNAr; i++)
		{
		double myh = fabs(mvAr[i] * theh);
			myh = MAX(myh, 1e-8);
			mvAr[i] += myh;
			UpdateProxyMeanParameters();
		double myF1 = ComputeMean(theDate, theValue);
			theGradData.mCurrentGradMu[theBegIndex + i] = (myF1 - myF0) / myh;
			mvAr[i] -= myh;
		}
		for (uint i = 0; i < myNMa; i++)
		{
		double myh = fabs(mvMa[i] * theh);
			myh = MAX(myh, 1e-8);
			mvMa[i] += myh;
			UpdateProxyMeanParameters();
		double myF1 = ComputeMean(theDate, theValue);
			theGradData.mCurrentGradMu[theBegIndex + i + myNAr] = (myF1 - myF0) / myh;
			mvMa[i] -= myh;
		}
	double myh = fabs(mvFracD * theh);
		myh = MAX(myh, 1e-8);
		mvFracD += myh;
		UpdateProxyMeanParameters();
	double myF1 = ComputeMean(theDate, theValue);
		theGradData.mCurrentGradMu[theBegIndex + myNAr + myNMa] = (myF1 - myF0) / myh;
		mvFracD -= myh;
		UpdateProxyMeanParameters();
		myF0 = ComputeMean(theDate, theValue);
	}

	void cArfima::ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cDVector& theGradTheta, cDVector& theGradU, double& theGradH, uint& theNu, uint& theNh)
	{
		theNh = 0;
		theNu = 0;
		theGradTheta = 0.0;
		for (uint i = 0; i < mvAr.GetSize() + mvMa.GetSize() + 1; i++)
			theGradTheta[theBegIndex + i] += mvGradPolAr[i].BackwardPolOp(theValue.mYt, theDate);
	}

	void cArfima::ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cFuncMeanAndVar& theDerivM)
	{
		theDerivM.mGradTheta = 0.0;
		for (uint i = 0; i < mvAr.GetSize() + mvMa.GetSize() + 1; i++)
			theDerivM.mGradTheta[theBegIndex + i] = mvGradPolAr[i].BackwardPolOp(theValue.mYt, theDate);
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

	void cArfima::ComputeHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex, cAbstResiduals* theResiduals)
	{
	uint myNAr = GetNAr();
	uint myNMa = GetNMa();
	uint myNParam = myNAr + myNMa + 1;
		for (uint i = 0; i < myNParam; i++)
			for (uint j = i; j < myNParam; j++)
			theHessData.mCurrentHessMu[i + theBegIndex][j + theBegIndex] = 0.0;

		for (uint i = 0; i < myNParam; i++)
			for (uint j = i; j < myNParam; j++)
			{
				theHessData.mCurrentHessMu[i + theBegIndex][j + theBegIndex] += mvHessPolAr[i][j].BackwardPolOp(theData.mUt, theDate);
				for (uint t = 0; t < mvNTruncLag; t++)
				{
					theHessData.mCurrentHessMu[i+theBegIndex][j+theBegIndex] -= mvGradPolAr[i][t + 1] * theGradData.mGradMt[t][j] + mvGradPolAr[j][t + 1] * theGradData.mGradMt[t][i];
					theHessData.mCurrentHessMu[i + theBegIndex][j + theBegIndex] -= mvPolAr[t + 1] * theHessData.mHessMt[t][i][j];
				}
				theHessData.mCurrentHessMu[theBegIndex + j][theBegIndex + i] = theHessData.mCurrentHessMu[theBegIndex + i][theBegIndex + j];
			}
	}
	
	void cArfima::ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, const cDVector& theGradTheta, const cDVector& theGradU, double theGradH, uint theNu, uint theNh, cDMatrix& theHessTheta2, cDVector* theHessThetaU, cDVector& theHessThetaH, cDVector* theHessU2, cDVector& theHessUH, double& theHessH2)
	{
		for (uint i = 0; i < mvAr.GetSize() + mvMa.GetSize() + 1; i++)
			for (uint j = i ; j < mvAr.GetSize() + mvMa.GetSize() + 1; j++)
				theHessTheta2[theBegIndex + i][theBegIndex + j] = theHessTheta2[theBegIndex + j][theBegIndex + i] = mvHessPolAr[i][j].BackwardPolOp(theValue.mYt, theDate);
	}

	void cArfima::ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, cFuncMeanAndVar& theDerivM)
	{
		for (uint i = 0; i < mvAr.GetSize() + mvMa.GetSize() + 1; i++)
			for (uint j = i; j < mvAr.GetSize() + mvMa.GetSize() + 1; j++)
			{
				theDerivM.mHessTheta2[theBegIndex + i][theBegIndex + j] += mvHessPolAr[i][j].BackwardPolOp(theValue.mYt, theDate);
				theDerivM.mHessTheta2[theBegIndex + j][theBegIndex + i] = theDerivM.mHessTheta2[theBegIndex + i][theBegIndex + j];
			}

	}

	void cArfima::ComputeHessF(uint theDate, const cRegArchValue& theData, cDMatrix& theHessF)
	{
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
