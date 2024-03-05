#include "StdAfxRegArchLib.h"
/*!
\file cFigarch.cpp
\brief sources for class cFigarch methods.

\author Jean-Baptiste DURAND, Ollivier TARAMASCO
\date dec-18-2006 - Last change feb-18-2011
*/
namespace RegArchLib {
	/*!
	* \fn cFigarch::cFigarch(uint theNArch, uint theNFigarch, double theFracD, uint theNTruncLag):cAbstCondVar(eFigarch)
	* \param uint theNArch: number of ARCH lags
	* \param uint theNFigarch: number of Figarch lags
	*/
	cFigarch::cFigarch(uint theNArch, uint theNGarch, double theFracD, uint theNTruncLag)
		:cAbstCondVar(eFigarch)  // call constructor of cAbstCondVar with type eFigarch
	{
		mvArch.ReAlloc(theNArch);
		mvGarch.ReAlloc(theNGarch);
		mvConst = 0.0;
		mvFracD = theFracD;
		mvNTruncLag = theNTruncLag;
	uint myNParam = theNArch + theNGarch + 2;
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

		MESS_CREAT("cFigarch");
	}

	/*!
	* \fn cFigarch::cFigarch(double theConst, cDVector& theArch, cDVector& theFigarch):cAbstCondVar(eFigarch)
	* \param double theConst: constant part of the Figarch(p, q) model
	* \param cDVector& theFigarch theArch: ARCH parameters
	* \param cDVector& theFigarch theFigarch: Figarch parameters
	*/
	cFigarch::cFigarch(double theConst, const cDVector& theArch, const cDVector& theGarch, double theFracD, uint theNTruncLag)
		:cAbstCondVar(eFigarch)
	{
		mvConst = theConst;
		mvArch = theArch;
		mvGarch = theGarch;
		mvFracD = theFracD;
		mvNTruncLag = theNTruncLag;
	uint myNArch = theArch.GetSize();
	uint myNGarch = theGarch.GetSize();
	uint myNParam = myNArch + myNGarch + 2;
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

		MESS_CREAT("cFigarch");
	}

	/*!
	* \fn cFigarch::cFigarch(const cAbsCondVar& theFigarch):cAbstCondVar(eFigarch)
	* \param cAbsCondVar& theEFigarch: theFigarch class
	*/
	cFigarch::cFigarch(const cFigarch& theFigarch)
		:cAbstCondVar(eFigarch)
	{
		*this = theFigarch;
		MESS_CREAT("cFigarch");
	}

	/*!
	* \fn cFigarch::~cFigarch()
	*/
	cFigarch::~cFigarch()
	{
	uint myNParam = mvArch.GetSize() + mvGarch.GetSize() + 2;
		mvArch.Delete();
		mvGarch.Delete();
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

		MESS_DESTR("cFigarch");
	}

	/*!
	* \fn cAbstCondVar* cFigarch::PtrCopy()
	*/
/*	cAbstCondVar* cFigarch::PtrCopy() const
	{
		//		 cConstCondVar *myConstCondVar = new cConstCondVar(*this);
		//		 return myConstCondVar;
		return cAbstCondVarPtrCopy<cFigarch>();
	}
*/
	/*!
	* \fn void cFigarch::Delete(void)
	* \param void
	* \details Free memory
	*/
	void cFigarch::Delete(void)
	{
	uint myNParam = mvArch.GetSize() + mvGarch.GetSize() + 2;
		mvArch.Delete();
		mvGarch.Delete();
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
	/*!
	* \fn void cFigarch::Print(ostream& theOut) const
	* \param ostream& theOut: the output stream, default cout.
	*/
	void cFigarch::Print(ostream& theOut) const
	{
	uint myNArch = mvArch.GetSize();
	uint myNGarch = mvGarch.GetSize();
		theOut << "FIGARCH(" << myNArch << ", d, " << myNGarch << ") model with:" << endl;
		theOut << "\tCste=" << mvConst << endl;
		for (uint i = 0; i < myNArch; i++)
			theOut << "\tARCH[" << i + 1 << "]=" << mvArch[i] << endl;
		for (uint j = 0; j < myNGarch; j++)
			theOut << "\tGARCH[" << j + 1 << "]=" << mvGarch[j] << endl;
		theOut << "\td=" << mvFracD << endl;
	}
#else
	/*!
	* \fn void cFigarch::Print()
	*/
	void cFigarch::Print()
	{
		uint myNArch = mvArch.GetSize();
		uint myNGarch = mvGarch.GetSize();
		Rprintf("FIGARCH(%d, d,%d) model with:", myNArch, myNGarch);
		Rprintf("\tCste=%f\n", mvConst);
		for (uint i = 0; i < myNArch; i++)
			Rprintf("\tARCH[%d]=%f\n", i + 1, mvArch[i]);
		for (uint j = 0; j < myNGarch; j++)
			Rprintf("\tGARCH[%d]=%f\n", j + 1, mvGarch[j]);
		Rprintf("\td=%f\n", mvFracD);
	}
#endif //_RDLL_

	void cFigarch::SetDefaultInitPoint(double theMean, double theVar)
	{
		mvConst = theVar*0.1;
	uint myNArch = mvArch.GetSize();
	uint myNFigarch = mvGarch.GetSize();
	uint i;
		for (i = 0; i < myNArch; i++)
			mvArch[i] = 0.1 / (double)myNArch;
		for (i = 0; i < myNFigarch; i++)
			mvGarch[i] = 0.8 / (double)myNFigarch;
		mvFracD = 0.2;
	}

	void cFigarch::SetDefaultInitPoint(cRegArchValue& theValue)
	{
	double myVar;
		theValue.ComputeVar(myVar);
		mvConst = myVar*0.1;
	uint myNArch = mvArch.GetSize();
	uint myNFigarch = mvGarch.GetSize();
	uint i;
		for (i = 0; i < myNArch; i++)
			mvArch[i] = 0.1 / (double)myNArch;
		for (i = 0; i < myNFigarch; i++)
			mvGarch[i] = 0.8 / (double)myNFigarch;
		mvFracD = 0.2;
	}

	/*!
	* \fn void cFigarch::ReAlloc(const uint theSize, const uint theNumParam)
	* \param const uint theSize: new size of mvArch or mvGarch
	* \param const uint theNumParam: 0 for mvArch, 1 for mvGarch.
	* \details new allocation of mvArch or mvGarch
	*/
	void cFigarch::ReAlloc(const uint theSize, const uint theNumParam)
	{
		switch (theNumParam)
		{
		case 1:
			mvArch.ReAlloc(theSize);
			break;
		case 2:
			mvGarch.ReAlloc(theSize);
			break;
		default:
//			throw cError("cFigarch::ReAlloc - theNumParam must be in 1, 2.");
			break;
		}
	}

	/*!
	* \fn void cFigarch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	* \param const cDVector& theVectParam: the vector of Const, ARCH or Figarch coefficients
	* \param const uint theNumParam: =0, the constant part; =1 the ARCH coefficients; =2 theFigarch Coefficients
	* \details new allocation of mvArch or mvConst
	*/
	void cFigarch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{
		switch (theNumParam)
		{
		case 0: // mvConst
			if (theVectParam.GetSize() > 0)
				mvConst = theVectParam[0];
			else
				throw cError("cFigarch::ReAlloc - Size of theVectParam must be > 0");
			break;
		case 1: // mvArch
			mvArch = theVectParam;
			break;
		case 2: // mvGarch
			mvGarch = theVectParam;
			break;
		case 3:
			mvFracD = theVectParam[0];
			break;
		default:
			throw cError("cFigarch::ReAlloc - theNumParam must be in 0 .. 3");
			break;
		}
	}

	/*!
	* \fn void cFigarch::Set(const double theValue, const uint theIndex, const uint theNumParam)
	* \brief fill the parameters vector
	* \param const double theValue: the value of the "theIndex" th lag. Default 0.
	* \param const uint theIndex: the index.
	* \param const uint theNumParam: =0, mvConst, =1, ARCH parameters; =2, Figarch parameters
	* \details mvArch[theIndex] = theValue or mvGarch[theIndex]= theValue or mvConst = theValue
	*/
	void cFigarch::Set(const double theValue, const uint theIndex, const uint theNumParam)
	{
		switch (theNumParam)
		{
		case 0:
			mvConst = theValue;
			break;
		case 1:
			if (theIndex < mvArch.GetSize())
				mvArch[theIndex] = theValue;
			else
				throw cError("cFigarch::Set - wrong index");
			break;
		case 2:
			if (theIndex < mvGarch.GetSize())
				mvGarch[theIndex] = theValue;
			else
				throw cError("cFigarch::Set - wrong index");
			break;
		case 3:
			mvFracD = theValue;
			break;
		default:
			throw cError("cFigarch::Set - theNumParam must be in 0 .. 3");
			break;
		}
	}

	/*!
	* \fn void cFigarch::Set(const cDVector& theVectParam, const uint theNumParam)
	* \brief fill the parameters vector
	* \param const cDVector& theVectParam: the vector of values
	* \param const uint theNumParam: =0, mvConst, =1, ARCH parameters; =2, Figarch parameters
	* \details mvAr = theValue
	*/
	void cFigarch::Set(const cDVector& theVectParam, const uint theNumParam)
	{
		switch (theNumParam)
		{
		case 0:
			if (theVectParam.GetSize() > 0)
				mvConst = theVectParam[0];
			else
				throw cError("cFigarch::Set - Size of theVectParam must be > 0");
			break;
		case 1:
			mvArch = theVectParam;
			break;
		case 2:
			mvGarch = theVectParam;
			break;
		case 3:
			if (theVectParam.GetSize() > 0)
				mvFracD = theVectParam[0];
			else
				throw cError("cFigarch::Set - Size of theVectParam must be > 0");
			break;
		default:
			throw cError("cFigarch::Set - theNumParam must be in 0, 1, 2");
			break;
		}
	}

	double  cFigarch::Get(const uint theIndex, const uint theNumParam)
	{
		switch (theNumParam)
		{
		case 0:
			return mvConst;
			break;
		case 1:
			return mvArch[theIndex];
			break;
		case 2:
			return mvGarch[theIndex];
			break;
		case 3:
			return mvFracD;
			break;
		}
	}

	cDVector& cFigarch::Get(const uint theNumParam)
	{
	cDVector *myAux;	
		switch (theNumParam)
		{
		case 0:
			myAux = new cDVector(1, mvConst);
			return *myAux;
			break;
		case 1:
			return mvArch;
			break;
		case 2:
			return mvGarch;
			break;
		case 3:
			myAux = new cDVector(1, mvFracD);
			return *myAux;
			break;
		}
	}

	cFigarch& cFigarch::operator =(const cFigarch& theSrc)
	{
		mvArch = theSrc.mvArch;
		mvGarch = theSrc.mvGarch;
		mvPolMa = theSrc.mvPolMa;
		SetCondVarType(eFigarch);
		mvConst = theSrc.mvConst;
		mvFracD = theSrc.mvFracD;
		mvNTruncLag = theSrc.mvNTruncLag;
	uint myNParam = mvArch.GetSize() + mvGarch.GetSize() + 2;
	uint myOldNParam = theSrc.GetNParam();

		mvPolMa = theSrc.mvPolMa;
		mvGradPolMa = new cPolynome[myNParam];
		mvHessPolMa = new cPolynome*[myNParam];
		for (uint i = 0; i < myNParam; i++)
			mvHessPolMa[i] = new cPolynome[myNParam];

		for (uint i = 0; i < MIN(myNParam, myOldNParam); i++)
		{
			mvGradPolMa[i] = theSrc.mvGradPolMa[i];
			for (uint j = 0; j < MIN(myNParam, myOldNParam); j++)
				mvHessPolMa[i][j] = theSrc.mvHessPolMa[i][j];
		}

		return *this;
	}

	void cFigarch::ReAllocProxyVarParameters(uint theNOldParam)
	{
		for (uint i = 0; i < theNOldParam; i++)
		{
			mvGradPolMa[i].Delete();
			for (uint j = 0; j < theNOldParam; j++)
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
	
	static void ComputePhiAndTeta(cFigarch& theFigarch, cPolynome& thePhi, cPolynome& theTeta)
	{
	uint myNArch = theFigarch.GetNArch();
	uint myNGarch = theFigarch.GetNGarch();
		thePhi.Resize(MAX(myNArch, myNGarch));
		theTeta.Resize(myNGarch);
		thePhi[0] = theTeta[0] = 1.0;
		for (uint i = 0; i < myNArch; i++)
			thePhi[i + 1] -= theFigarch.Get(i, 1);
		for (uint j = 0; j < myNGarch; j++)
		{
			thePhi[j + 1] -= theFigarch.Get(j, 2);
			theTeta[j + 1] = -theFigarch.Get(j, 2);
		}
	}

	static void ComputePolForGrad(cPolynome& thePhi, cPolynome& theTeta, double theD, uint theNMaxLag, cPolynome& theARMAPol, cPolynome& thePolGradAlpha, cPolynome& thePolGradBeta, cPolynome& thePolGradD)
	{
		cPolynome myDeltaD(theNMaxLag + 2);
		ComputeDeltaPowD(theD, theNMaxLag + 2, myDeltaD);
		cPolynome myPolAux1(theNMaxLag + 2), myRest(0);
		IncrPowDiv(myDeltaD, theTeta, theNMaxLag + 2, myPolAux1, myRest);
		thePolGradAlpha = -1 * myPolAux1;
		cPolynome myPolAux2(theNMaxLag + 2);
		IncrPowDiv(myPolAux1, theTeta, theNMaxLag + 2, myPolAux2, myRest);
		thePolGradBeta = myPolAux2 - myPolAux1;
		cPolynome myLogDeltaD(theNMaxLag + 2);
		for (uint i = 1; i <= theNMaxLag + 2; i++)
			myLogDeltaD[i] = -1 / (double)i;
		cPolynome myPolAux3(theNMaxLag + 2), myPolAux4(theNMaxLag + 2);
		myPolAux3 = TrunkMult(myPolAux1, thePhi, theNMaxLag + 2);
		myPolAux4 = TrunkMult(myPolAux3, myLogDeltaD, theNMaxLag + 2);
		theARMAPol = -1 * myPolAux3;
		thePolGradD = -1 * myPolAux4;
	}

	static void UpdateAllArchParam(cArch* thecArch, double theConst, cDVector& theARCH, cDVector& theGARCH, uint theNLag, double theFracD)
	{
	cPolynome myDeltaD(theNLag + 2);
		ComputeDeltaPowD(theFracD, theNLag + 2, myDeltaD);
	uint myNArch = theARCH.GetSize();
	uint myNGarch = theGARCH.GetSize();
	cPolynome myPhi(MAX(myNArch, myNGarch));
		for (uint i = 1; i <= myNArch; i++)
			myPhi[i] -= theARCH[i - 1];
		for (uint i = 1; i <= myNGarch; i++)
			myPhi[i] -= theGARCH[i - 1];
		myPhi[0] = 1;
	cPolynome myAux(theNLag + 2);
		myAux = TrunkMult(myPhi, myDeltaD, theNLag + 2);
	cPolynome myTeta(myNGarch);
		for (uint i = 1; i <= myNGarch; i++)
			myTeta[i] = -theGARCH[i - 1];
		myTeta[0] = 1;
	cPolynome myQuot(theNLag), myRest(theNLag);
		IncrPowDiv(myAux, myTeta, theNLag, myQuot, myRest);
		for (uint i = 1; i <= theNLag; i++)
			thecArch->Set(-myQuot[i], i - 1, 1);
	double mySum = Sum(theGARCH);
	double myOmega0 = theConst / (1 - mySum);
		thecArch->Set(myOmega0, 0, 0);
	}

	void cFigarch::UpdateProxyVarParameters(void)
	{
	uint myNArch = mvArch.GetSize();
	uint myNGarch = mvGarch.GetSize();
	uint myNParam = myNArch + myNGarch + 2;
	cPolynome myPhi(0), myTeta(0);
		ComputePhiAndTeta(*this, myPhi, myTeta);
	cPolynome myDeltaPD(0);
		ComputeDeltaPowD(mvFracD, mvNTruncLag, myDeltaPD);
	cPolynome myPhiDeltaPD(0);
		myPhiDeltaPD = TrunkMult(myPhi, myDeltaPD, mvNTruncLag);
	cPolynome myReste(0);
	cPolynome myUn(0);
		myUn[0] = 1.0;
	cPolynome myTetam1(0);
		IncrPowDiv(myUn, myTeta, mvNTruncLag, myTetam1, myReste);
	cPolynome myTetam1PhiDeltaPD(0);
		myTetam1PhiDeltaPD = TrunkMult(myPhiDeltaPD, myTetam1, mvNTruncLag);
		mvPolMa = myUn - myTetam1PhiDeltaPD;
	double myS = (1.0 - Sum(mvGarch));
	double myOmega0 = mvConst / myS;
		mvPolMa[0] = myOmega0;
	cPolynome myTetam2PhiDeltaPD(0);
		myTetam2PhiDeltaPD = TrunkMult(myTetam1, myTetam1PhiDeltaPD, mvNTruncLag);
	cPolynome myTetam1DeltaPD(0);
		myTetam1DeltaPD = TrunkMult(myTetam1, myDeltaPD, mvNTruncLag);
	cPolynome myLogDelta(0);
		ComputeLogDelta(mvNTruncLag, myLogDelta);
	cPolynome myTetam1PhiDeltaPDLogDelta(0);
		myTetam1PhiDeltaPDLogDelta = TrunkMult(myLogDelta, myTetam1PhiDeltaPD, mvNTruncLag);
	cPolynome myX(1);
		myX[0] = 0.0;
		myX[1] = 1.0;
	cPolynome* myMult = new cPolynome[myNParam + 1];
		myMult[0] = myUn;
		for (uint i = 0; i < myNParam; i++)
		{
			mvGradPolMa[i] = cPolynome(0);
			myMult[i + 1] = myMult[i] * myX;
		}
		for (uint k = 0; k < myNParam; k++)
		{
			for (uint l = 0; l < myNParam; l++)
				mvHessPolMa[k][l].mCoeff = 0.0;
		}
		// Pol Grad ARCH
		for (uint k = 1; k <= myNArch; k++)
		{
			mvGradPolMa[k] = myMult[k] * myTetam1DeltaPD;
		}
	// Pol Grad Garch
		for (uint k = 1; k <= myNGarch; k++)
		{
			mvGradPolMa[k+myNArch] += myMult[k] * myTetam1DeltaPD;
		}
		for (uint k = 1; k <= myNGarch; k++)
		{
			mvGradPolMa[k+myNArch] -= myMult[k] * myTetam2PhiDeltaPD;

		}
	// Pol Grad FracD
		mvGradPolMa[myNParam - 1] -=  myTetam1PhiDeltaPDLogDelta;

		// Final Grad
	double myS2 = myS * myS;
	double myAuxDouble = mvConst / myS2;
	cDVector myDOmega0(myNParam);
		myDOmega0[0] = 1/myS;
		for (uint k = myNArch + 1; k < myNArch + myNGarch + 1; k++)
			myDOmega0[k] = myAuxDouble;
		for (uint k = 0; k < myNParam; k++)
		{
		cPolynome myAuxPol(0);
			myAuxPol = TrunkPoly(mvGradPolMa[k], mvNTruncLag);
			myAuxPol[0] = myDOmega0[k];
			mvGradPolMa[k] = myAuxPol;
		}
	// Hessien
	cDMatrix myD2Omega0(myNParam, myNParam,0.0);
	double myS3 = myS * myS2;
	double myAux1 = 1 / myS2;
	double myAux2 = 2 * mvConst / myS3;
		for (uint j = myNArch+1 ; j < myNArch + myNGarch + 1 ; j++)
		{
			myD2Omega0[0][j] = myD2Omega0[j][0] = myAux1;
			for (uint i = myNArch + 1; i < myNArch + myNGarch + 1; i++)
				myD2Omega0[i][j] = myD2Omega0[j][i] = myAux2;
		}
	cPolynome myTetam2DeltaPD(0);
		myTetam2DeltaPD = TrunkMult(myTetam1DeltaPD, myTetam1, mvNTruncLag);
	cPolynome myTetam2PhiDeltaPDLogDelta(0);
		myTetam2PhiDeltaPDLogDelta = TrunkMult(myTetam2PhiDeltaPD, myLogDelta, mvNTruncLag);
	cPolynome myTetam1DeltaPDLogDelta(0);
		myTetam1DeltaPDLogDelta = TrunkMult(myTetam1DeltaPD, myLogDelta, mvNTruncLag);
	cPolynome myTetam3PhiDeltaPD(0);
		myTetam3PhiDeltaPD = TrunkMult(myTetam1, myTetam2PhiDeltaPD, mvNTruncLag);
	cPolynome myTetam1PhiDeltaPDLogDeltaP2(0);
		myTetam1PhiDeltaPDLogDeltaP2 = TrunkMult(myTetam1PhiDeltaPDLogDelta, myLogDelta, mvNTruncLag);
	// dPhi*dTeta + dTeta*dPhi
		for (uint i = myNArch + 1; i < myNArch + myNGarch + 1; i++)
		{
		int myInd = i - myNArch + 1;
			for (uint j = 1; j < myNArch + myNGarch + 1; j++)
			{
			cPolynome myAuxPol = myTetam2DeltaPD * myMult[myInd];
				mvHessPolMa[i][j] += myAuxPol;
				mvHessPolMa[j][i] = mvHessPolMa[i][j];
				if (j == myNArch)
					myInd = i - myNArch + 1;
				else
					myInd++;
			}
		}

		//dTeta*dFracD + dFracD*dTeta
		for (uint i = myNArch + 1; i < myNArch + myNGarch + 1; i++)
		{
		cPolynome myAuxPol = myTetam2PhiDeltaPDLogDelta * myMult[i-myNArch];
			mvHessPolMa[i][myNParam - 1] -= myAuxPol;
			mvHessPolMa[myNParam - 1][i] -= myAuxPol;
		}
	//dPhi*dFracD + dFracd*dPhi
		for (uint i = 1 ; i < myNArch + myNGarch + 1; i++)
		{
		int myInd = 1;
		cPolynome myAuxPol = myTetam1DeltaPDLogDelta * myMult[myInd];
			mvHessPolMa[i][myNParam - 1] += myAuxPol;
			mvHessPolMa[myNParam - 1][i] += myAuxPol;
			if (i == myNArch)
				myInd = 1;
			else
				myInd++;
		}
	//dTeta*dTeta
		for (uint i = myNArch + 1; i < myNArch + myNGarch + 1; i++)
		{
		int myInd = i - myNArch + 1;
			for (uint j = 0; j < myNArch + myNGarch + 1; j++)
			{
			cPolynome myAuxPol = myTetam3PhiDeltaPD * myMult[myInd];
				myAuxPol *= 2.0;
				mvHessPolMa[i][j] -= myAuxPol;
//				if (i != j)
//					mvHessPolMa[j][i] -= myAuxPol;
				myInd++;
			}
		}

	//dFracD * dFracD
		mvHessPolMa[myNParam - 1][myNParam - 1] -= myTetam1PhiDeltaPDLogDeltaP2;
	// Final Hessien	
		for (uint i = 0; i < myNParam; i++)
			for (uint j = i; j < myNParam; j++)
			{
				mvHessPolMa[i][j] = TrunkPoly(mvHessPolMa[i][j], mvNTruncLag);
				(mvHessPolMa[i][j])[0] = myD2Omega0[i][j];
				(mvHessPolMa[j][j])[0] = (mvHessPolMa[i][j])[0];
			}
		delete[] myMult;

    }

	/*!
	* \fn double cFigarch::ComputeVar(uint theDate, const cReFigarchValue& theData) const
	* \param int theDate: date of computation
	* \param const cReFigarchValue& theData: past datas
	* \details theData is not updated here.
	*/
	double cFigarch::ComputeVar(uint theDate, const cRegArchValue& theData) const
	{
		double myVar =  mvPolMa.BackwardPolOp(theData.mUt, theDate, 2);
		return myVar;

	}

	uint cFigarch::GetNParam(void) const
	{
		return 2 + mvArch.GetSize() + mvGarch.GetSize();
	}

	uint cFigarch::GetNArch(void) const
	{
		return mvArch.GetSize();
	}

	uint cFigarch::GetNGarch(void) const
	{
		return mvGarch.GetSize();
	}

	uint cFigarch::GetNLags(void) const
	{
	uint myMax = MAX(mvArch.GetSize(), mvGarch.GetSize());
		myMax = MAX(myMax, mvNTruncLag);
		return  myMax;
	}

	uint cFigarch::GetNu(void) const
	{
		return mvArch.GetSize();
	}

	uint cFigarch::GetNh(void) const
	{
		return MAX(mvGarch.GetSize(), mvArch.GetSize());
	}

	void cFigarch::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, cAbstResiduals* theResiduals)
	{
	uint myNArch = mvArch.GetSize();
	uint myNGarch = mvGarch.GetSize();
	uint myNParam = myNArch + myNGarch + 2;
	uint myBegIndex = theGradData.GetNMeanParam();
		theGradData.mCurrentGradVar = 0.0;
		for (uint i = 0; i < myNParam; i++)
			theGradData.mCurrentGradVar[myBegIndex + i] = mvGradPolMa[i].BackwardPolOp(theValue.mUt, theDate, 2);			
		for (uint t = 1; t < MIN(theDate, mvNTruncLag+1); t++)
			theGradData.mCurrentGradVar -= 2*mvPolMa[t] * theValue.mUt[theDate-t]*theGradData.mGradMt[t-1];
	}

	void cFigarch::ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cFuncMeanAndVar& theDerivM, cAbstResiduals* theResids)
	{
	}

	void  cFigarch::NumericComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, uint theBegIndex, cAbstResiduals* theResiduals, double theh)
	{
	uint myNArch = mvArch.GetSize();
	uint myNGarch = mvGarch.GetSize();
		UpdateProxyVarParameters();
		double myF0 = ComputeVar(theDate, theValue);
		double myh = fabs(mvConst * theh);
		myh = MAX(myh, 1e-8);
		mvConst += myh;
		UpdateProxyVarParameters();
		double myF1 = ComputeVar(theDate, theValue);
		theGradData.mCurrentGradVar[theBegIndex] += (myF1 - myF0) / myh;
		mvConst -= myh;

		for (uint i = 0; i < myNArch; i++)
		{
			myh = fabs(mvArch[i] * theh);
			myh = MAX(myh, 1e-8);
			mvArch[i] += myh;
			UpdateProxyVarParameters();
			myF1 = ComputeVar(theDate, theValue);
			theGradData.mCurrentGradVar[theBegIndex + i + 1] += (myF1 - myF0) / myh;
			mvArch[i] -= myh;
		}
		for (uint i = 0; i < myNGarch; i++)
		{
			myh = fabs(mvGarch[i] * theh);
			myh = MAX(myh, 1e-8);
			mvGarch[i] += myh;
			UpdateProxyVarParameters();
			myF1 = ComputeVar(theDate, theValue);
			theGradData.mCurrentGradVar[theBegIndex + i + myNArch + 1] += (myF1 - myF0) / myh;
			mvGarch[i] -= myh;
		}
		myh = fabs(mvFracD * theh);
		myh = MAX(myh, 1e-8);
		mvFracD += myh;
		UpdateProxyVarParameters();
		myF1 = ComputeVar(theDate, theValue);
		theGradData.mCurrentGradVar[theBegIndex + myNArch + myNGarch + 1] += (myF1 - myF0) / myh;
		mvFracD -= myh;
	}

	void cFigarch::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{
	uint mySize = GetNParam();
		if (theDestVect.GetSize() < mySize + theIndex)
			throw cError("Wrong size");
		theDestVect[theIndex] = mvConst;
		mvArch.SetSubVectorWithThis(theDestVect, theIndex + 1);
		mvGarch.SetSubVectorWithThis(theDestVect, theIndex + 1 + mvArch.GetSize());
		theDestVect[theIndex + 1 + mvArch.GetSize() + mvGarch.GetSize()] = mvFracD;
	}

	void cFigarch::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
		uint mySize = theSrcVect.GetSize();
		if (GetNParam() + theIndex > mySize)
			throw cError("Wrong size");
		mvConst = theSrcVect[theIndex];
		mvArch.SetThisWithSubVector(theSrcVect, theIndex + 1);
		mvGarch.SetThisWithSubVector(theSrcVect, theIndex + 1 + mvArch.GetSize());
		mvFracD = theSrcVect[theIndex + 1 + mvArch.GetSize() + mvGarch.GetSize()];
	}

	void cFigarch::ComputeHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
		uint myNArch = mvArch.GetSize();
		uint myNGarch = mvGarch.GetSize();
		uint myNParam = myNArch + myNGarch + 2;
		uint myBegIndex = theGradData.GetNMeanParam();
 		for (uint i = 0; i < myNParam; i++)
			for (uint j = i; j < myNParam; j++) 
				theHessData.mCurrentHessVar[myBegIndex + i][myBegIndex + j] = mvHessPolMa[i][j].BackwardPolOp(theData.mUt, theDate, 2);
		for (uint k = 0; k < myNParam; k++)
			for (uint l = k; l < myNParam; l++)
				theHessData.mCurrentHessVar[myBegIndex + l][myBegIndex + k] = theHessData.mCurrentHessVar[myBegIndex + k][myBegIndex + l];
		cDVector myGradPol(theGradData.GetNParam());
		for (uint t = 1; t < MIN(theDate, mvNTruncLag + 1); t++)
		{
			for (uint i = 1; i < myNParam; i++)
				myGradPol[i+myBegIndex] = mvGradPolMa[i][t];
			theHessData.mCurrentHessVar -= 2 * theData.mUt[theDate - t] * (myGradPol*Transpose(theGradData.mGradMt[t - 1]) + theGradData.mGradMt[t - 1] * Transpose(myGradPol));
			theHessData.mCurrentHessVar += 2 * mvPolMa[t] * theGradData.mGradMt[t - 1] * Transpose(theGradData.mGradMt[t - 1]);
			theHessData.mCurrentHessVar -= 2 * mvPolMa[t] * theData.mUt[theDate - t] * theHessData.mHessMt[t-1];
		}
	}

	void cFigarch::ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, cFuncMeanAndVar& theDerivM, cAbstResiduals* theResids)
	{

	}

	void cFigarch::GetParamName(uint theIndex, char** theName)
	{
		uint myIndex = theIndex;
		sprintf(theName[myIndex++], "CST VAR");
		for (uint i = 0; i < mvArch.GetSize(); i++)
		{
			sprintf(theName[myIndex++], "ARCH[%d]", i + 1);

		}
		for (uint i = 0; i < mvGarch.GetSize(); i++)
		{
			sprintf(theName[myIndex++], "GARCH[%d]", i + 1);

		}
		sprintf(theName[myIndex++], "Frac. d");
	}
	 
	void cFigarch::GetParamName(uint theIndex, string theName[])
	{
	uint myIndex = theIndex;
	char myChar[100];
		sprintf(myChar, "CST VAR");
		theName[myIndex++] = myChar;
		for (uint i = 0; i < mvArch.GetSize(); i++)
		{
			sprintf(myChar, "ARCH[%d]", i + 1);
			theName[myIndex++] = myChar;

		}
		for (uint i = 0; i < mvGarch.GetSize(); i++)
		{
			sprintf(myChar, "GARCH[%d]", i + 1);
			theName[myIndex++] = myChar;

		}
		sprintf(myChar, "Frac. d");
		theName[myIndex++] = myChar;
	}

}//namespace
