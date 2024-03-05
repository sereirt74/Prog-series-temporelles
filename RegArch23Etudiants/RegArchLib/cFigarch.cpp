#include "StdAfxRegArchLib.h"
#ifdef _DEBUG
#include <Windows.h>
#endif // _DEBUG

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

		MESS_CREAT("cFigarch");
	}

	/*!
	* \fn cFigarch::cFigarch(const cAbsCondVar& theFigarch):cAbstCondVar(eFigarch)
	* \param cAbsCondVar& theEFigarch: theFigarch class
	*/
	cFigarch::cFigarch(const cFigarch& theFigarch)
		:cAbstCondVar(eFigarch)
	{
#ifdef _CLI_
		uint myNArch = theFigarch.mvArch.GetSize();
		uint myNGarch = theFigarch.mvGarch.GetSize();

//		cout << "ly 1" << endl;

		mvArch.ReAlloc(myNArch);
		mvGarch.ReAlloc(myNGarch);


//		cout << "ly 2" << endl;

		mvNTruncLag = theFigarch.mvNTruncLag;;
		uint myNParam = myNArch + myNGarch + 2;
		mvPolAr.Resize(mvNTruncLag);
		mvGradPolAr = new cPolynome[myNParam];
		mvHessPolAr = new cPolynome * [myNParam];
		for (uint i = 0; i < myNParam; i++)
		{
			mvGradPolAr[i].Resize(mvNTruncLag);
			mvHessPolAr[i] = new cPolynome[myNParam];
			for (uint j = 0; j < myNParam; j++)
				mvHessPolAr[i][j].Resize(mvNTruncLag);
		}
#endif _CLI_

#ifdef _CLI_
//		cout << "ly 3" << endl;
#endif // _CLI_
		*this = theFigarch;
#ifdef _CLI_
//		cout << "ly 4" << endl;
#endif // _CLI_
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
	void cFigarch::DeletePoly(void)
	{
		uint myNParam = mvArch.GetSize() + mvGarch.GetSize() + 2;
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

	void cFigarch::Delete(void)
	{
		mvArch.Delete();
		mvGarch.Delete();
		DeletePoly();
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

	void cFigarch::ReSizePoly(const uint theNParam)
	{
		DeletePoly();
		mvPolAr.Resize(mvNTruncLag);
		mvGradPolAr = new cPolynome[theNParam];
		mvHessPolAr = new cPolynome * [theNParam];
		for (uint i = 0; i < theNParam; i++)
		{
			mvGradPolAr[i].Resize(mvNTruncLag);
			mvHessPolAr[i] = new cPolynome[theNParam];
			for (uint j = 0; j < theNParam; j++)
				mvHessPolAr[i][j].Resize(mvNTruncLag);
		}
	}

	/*!
	* \fn void cFigarch::ReAlloc(const uint theSize, const uint theNumParam)
	* \param const uint theSize: new size of mvArch or mvGarch
	* \param const uint theNumParam: 0 for mvArch, 1 for mvGarch.
	* \details new allocation of mvArch or mvGarch
	*/
	void cFigarch::ReAlloc(const uint theSize, const uint theNumParam)
	{
		cDVector myArch, myGarch;
		uint myNParam = 0;
		switch (theNumParam)
		{
		case 1:
			myNParam = mvGarch.GetSize() + theSize + 2;
			ReSizePoly(myNParam);
			myArch.ReAlloc(theSize);
			mvArch = myArch;
			break;
		case 2:
			myNParam = mvArch.GetSize() + theSize + 2;
			ReSizePoly(myNParam);
			myGarch.ReAlloc(theSize);
			mvGarch = myGarch;
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
		uint myNParam = 0;
		switch (theNumParam)
		{
		case 0: // mvConst
			if (theVectParam.GetSize() > 0)
				mvConst = theVectParam[0];
			else
				throw cError("cFigarch::ReAlloc - Size of theVectParam must be > 0");
			break;
		case 1: // mvArch
			myNParam = mvGarch.GetSize() + theVectParam.GetSize() + 2;
			ReSizePoly(myNParam);
			mvArch = theVectParam;
			break;
		case 2: // mvGarch
			myNParam = theVectParam.GetSize() + mvArch.GetSize() + 2;
			ReSizePoly(myNParam);
			mvGarch = theVectParam;
			break;
		case 3:
			if (theVectParam.GetSize() > 0)
				mvFracD = theVectParam[0];
			else
				throw cError("cFigarch::ReAlloc - Size of theVectParam must be > 0");
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
		case 4:
			mvNTruncLag = theValue;
			ReSizePoly(mvArch.GetSize() + mvGarch.GetSize() + 2);
		break;
		default:
			throw cError("cFigarch::Set - theNumParam must be in 0 .. 4");
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
		case 4:
			if (theVectParam.GetSize() > 0)
			{
				mvNTruncLag = theVectParam[0];
				ReSizePoly(mvArch.GetSize() + mvGarch.GetSize() + 2);
			}
			else
				throw cError("cFigarch::Set - Size of theVectParam must be > 0");
			break;
		default:
			throw cError("cFigarch::Set - theNumParam must be in 0, ... 4");
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
		case 4:
			return mvNTruncLag;
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
		case 4:
			myAux = new cDVector(1, (double)mvNTruncLag);
			return *myAux;
			break;
		}
	}

	cFigarch& cFigarch::operator =(const cFigarch& theSrc)
	{

		mvArch = theSrc.mvArch;
		mvGarch = theSrc.mvGarch;
		DeletePoly();
		mvPolAr = theSrc.mvPolAr;
		SetCondVarType(eFigarch);
		mvConst = theSrc.mvConst;
		mvFracD = theSrc.mvFracD;
		mvNTruncLag = theSrc.mvNTruncLag;
	uint myNParam = theSrc.mvArch.GetSize() + theSrc.mvGarch.GetSize() + 2;

		mvGradPolAr = new cPolynome[myNParam];
		mvHessPolAr = new cPolynome* [myNParam];
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

	void cFigarch::ReAllocProxyVarParameters(uint theNOldParam)
	{
		for (uint i = 0; i < theNOldParam; i++)
		{
			mvGradPolAr[i].Delete();
			for (uint j = 0; j < theNOldParam; j++)
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
		cPolynome myAlpha(myNArch), myBeta(myNGarch);
		myAlpha[0] = myBeta[0] = 0;
		for (uint i = 1; i <= myNArch; i++)
			myAlpha[i] = mvArch[i - 1];
		for (uint j = 1; j <= myNGarch; j++)
			myBeta[j] = mvGarch[j - 1];
		cPolynome myOne(0);
		myOne[0] = 1;
		cPolynome myAlphapBeta(MAX(myNArch, myNGarch));
		myAlphapBeta = myAlpha + myBeta;
		cPolynome myPsi(MAX(myNArch, myNGarch));
		myPsi = myOne - myAlphapBeta;
		cPolynome myDeltaPD(0);
		ComputeDeltaPowD(mvFracD, mvNTruncLag, myDeltaPD);
		cPolynome myPsiDeltaPD(mvNTruncLag);
		myPsiDeltaPD = TrunkMult(myPsi, myDeltaPD, mvNTruncLag);
		cPolynome myOnemBeta(myNGarch);
		myOnemBeta = myOne - myBeta;
		cPolynome myOnemBetam1(0), myReste(0);
		IncrPowDiv(myOne, myOnemBeta, mvNTruncLag, myOnemBetam1, myReste);
		cPolynome myOnemBetam1PsiDeltaPD(mvNTruncLag);
		myOnemBetam1PsiDeltaPD = TrunkMult(myOnemBetam1, myPsiDeltaPD, mvNTruncLag);
		cPolynome myLambda(mvNTruncLag);
		myLambda = myOne - myOnemBetam1PsiDeltaPD;
		double myS = (1 - Sum(mvGarch));
		double myOmega0 = mvConst / myS;
		mvPolAr = myLambda;
		mvPolAr.Resize(mvNTruncLag);
		mvPolAr[0] = myOmega0;
		cPolynome* myGradPsi = new cPolynome[myNParam];
		for (int n = 0; n < myNParam; n++)
			myGradPsi[n] = cPolynome(mvNTruncLag);

		for (uint i = 1; i <= myNArch; i++)
		{
			cPolynome myXpi(i);
			myXpi[i] = -1;
			myGradPsi[i] = myXpi;
		}
		for (uint j = 1; j <= myNGarch; j++)
		{
			cPolynome myXpj(j);
			myXpj[j] = -1;
			myGradPsi[myNArch + j] = myXpj;
		}
		cPolynome* myGradLambda = new cPolynome[myNParam];
		for (int n = 0; n < myNParam; n++)
			myGradLambda[n] = cPolynome(mvNTruncLag);
		cPolynome myOnemBetam2(mvNTruncLag);
		IncrPowDiv(myOnemBetam1, myOnemBeta, mvNTruncLag, myOnemBetam2, myReste);
		cPolynome myOnemBetam2PsiDeltaPD(mvNTruncLag);
		myOnemBetam2PsiDeltaPD = TrunkMult(myOnemBetam2, myPsiDeltaPD, mvNTruncLag);
		for (uint j = 1; j <= myNGarch; j++)
		{
			cPolynome myXpj(j);
			myXpj[j] = 1;
			myGradLambda[myNArch + j] -= TrunkMult(myXpj, myOnemBetam2PsiDeltaPD, mvNTruncLag);
		}
		cPolynome myOnemBetam1DeltaPD(mvNTruncLag);
		myOnemBetam1DeltaPD = TrunkMult(myOnemBetam1, myDeltaPD, mvNTruncLag);
		for (uint i = 1; i <= myNArch + myNGarch; i++)
			myGradLambda[i] -= TrunkMult(myOnemBetam1DeltaPD, myGradPsi[i], mvNTruncLag);
		cPolynome myLogDelta(0);
		ComputeLogDelta(mvNTruncLag, myLogDelta);
		cPolynome myOnemBetam1PsiDeltaPDLogDelta(mvNTruncLag);
		myOnemBetam1PsiDeltaPDLogDelta = TrunkMult(myOnemBetam1PsiDeltaPD, myLogDelta, mvNTruncLag);
		myGradLambda[myNArch + myNGarch + 1] -= myOnemBetam1PsiDeltaPDLogDelta;

		for (uint i = 0; i <= myNArch + myNGarch + 1; i++)
			mvGradPolAr[i] = myGradLambda[i];
		mvGradPolAr[0][0] = 1.0 / myS;

		for (uint j = 1; j <= myNGarch; j++)
			mvGradPolAr[myNArch + j][0] = mvConst / pow(myS, 2.0);

		cPolynome** myHessLambda = new cPolynome*[myNParam];
		for (uint i = 0; i < myNParam; i++)
			myHessLambda[i] = new cPolynome[myNParam];
		for (uint i = 0; i < myNParam; i++)
			for (uint j = 0; j < myNParam; j++)
				myHessLambda[i][j] = cPolynome(mvNTruncLag);

		cPolynome myOnemBetam3(mvNTruncLag);
		IncrPowDiv(myOnemBetam2, myOnemBeta, mvNTruncLag, myOnemBetam3, myReste);
		cPolynome myOnemBetam3PsiDeltaPD(mvNTruncLag);
		myOnemBetam3PsiDeltaPD = TrunkMult(myOnemBetam3, myPsiDeltaPD, mvNTruncLag);
		for (uint i = 1; i <= myNGarch; i++)
		{
			for (uint j = 1; j <= myNGarch; j++)
			{
				cPolynome myXipj(i + j);
				myXipj[i + j] = -2.0;
				myHessLambda[i + myNArch][j + myNArch] += TrunkMult(myXipj, myOnemBetam3PsiDeltaPD, mvNTruncLag);
			}
		}

		cPolynome myOnemBetam2DeltaPD(mvNTruncLag);
		myOnemBetam2DeltaPD = TrunkMult(myOnemBetam2, myDeltaPD, mvNTruncLag);
		for (uint i = 1; i <= myNGarch; i++)
		{
			for (uint j = 1; j <= myNArch; j++)
			{
				cPolynome myXpipj(i + j);
				myXpipj[i + j] = 1.0;
				cPolynome myAuxPol = TrunkMult(myXpipj, myOnemBetam2DeltaPD, mvNTruncLag);
				myHessLambda[i+myNArch][j] += myAuxPol;
				myHessLambda[j][i + myNArch] += myAuxPol;
			}
			for (uint j = 1; j <= myNGarch; j++)
			{
				cPolynome myXpipj(i + j);
				myXpipj[i + j] = 1.0;
				myHessLambda[i+myNArch][j + myNArch] += TrunkMult(myXpipj, myOnemBetam2DeltaPD, mvNTruncLag);
				myHessLambda[j + myNArch][i + myNArch]+= TrunkMult(myXpipj, myOnemBetam2DeltaPD, mvNTruncLag);
			}
		}

/*
		for (uint i = 1; i <= myNGarch; i++)
		{
			for (uint j = 1; j <= myNArch; j++)
				myHessLambda[j][i + myNArch] = myHessLambda[i + myNArch][j];
			for (uint j = 1; j <= myNGarch; j++)
				myHessLambda[j + myNArch][i + myNArch] = myHessLambda[i + myNArch][j + myNArch];
		}
*/
		cPolynome myPsiDeltaPDLogDelta(mvNTruncLag);

		myPsiDeltaPDLogDelta = TrunkMult(myPsiDeltaPD, myLogDelta, mvNTruncLag);
		cPolynome myOnemBetam2PsiLogDeltaDeltaPD = TrunkMult(myOnemBetam2, myPsiDeltaPDLogDelta, mvNTruncLag);
		for (uint i = 1; i <= myNGarch; i++)
		{
			cPolynome myXpi(i);
			myXpi[i] = -1.0;

			myHessLambda[i + myNArch][1 + myNArch + myNGarch] += TrunkMult(myXpi, myOnemBetam2PsiLogDeltaDeltaPD, mvNTruncLag);
			myHessLambda[1 + myNArch + myNGarch][i + myNArch] = myHessLambda[i + myNArch][1 + myNArch + myNGarch];
		}
		cPolynome myOnemBetam1LogDelta(mvNTruncLag);
		myOnemBetam1LogDelta = TrunkMult(myOnemBetam1, myLogDelta, mvNTruncLag);
		cPolynome myOnemBetam1LogDeltaDeltaPD(mvNTruncLag);
		myOnemBetam1LogDeltaDeltaPD = TrunkMult(myDeltaPD, myOnemBetam1LogDelta, mvNTruncLag);
		for (uint i = 1; i <= myNArch; i++)
		{
			cPolynome myXpi(i);
			myXpi[i] = 1.0;
			myHessLambda[i][1 + myNArch + myNGarch] += TrunkMult(myXpi, myOnemBetam1LogDeltaDeltaPD, mvNTruncLag);
		}
		for (uint i = 1; i <= myNGarch; i++)
		{
			cPolynome myXpi(i);
			myXpi[i] = 1.0;
			myHessLambda[i + myNArch][1 + myNArch + myNGarch] += TrunkMult(myXpi, myOnemBetam1LogDeltaDeltaPD, mvNTruncLag);
		}
		for (uint i = 1; i <= myNArch + myNGarch; i++)
			myHessLambda[1 + myNArch + myNGarch][i] = myHessLambda[i][1 + myNArch + myNGarch];

		cPolynome myOnemBetam1LogDeltap2DeltaPD(mvNTruncLag);
		myOnemBetam1LogDeltap2DeltaPD = TrunkMult(myLogDelta, myOnemBetam1LogDeltaDeltaPD, mvNTruncLag);
		cPolynome  myOnemBetam1LogDeltap2PsiDeltaPD(mvNTruncLag);
		myOnemBetam1LogDeltap2PsiDeltaPD = TrunkMult(myOnemBetam1LogDeltap2DeltaPD, myPsi, mvNTruncLag);

		myHessLambda[myNArch + myNGarch + 1][myNArch + myNGarch + 1] -= myOnemBetam1LogDeltap2PsiDeltaPD;
		double myS2 = myS * myS;
		double myS3 = myS2 * myS;
		cPolynome myCstPol1(0);
		myCstPol1[0] = 1.0 / myS2;
		cPolynome myCstPol2(0);
			myCstPol2[0] = 2 * mvConst / myS3;
		for (uint i = 1; i <= myNGarch; i++)
		{
			myHessLambda[0][i + myNArch] += myCstPol1;
			myHessLambda[i + myNArch][0] += myCstPol1;
			for (uint j = i; j <= myNGarch; j++)
			{
				myHessLambda[j + myNArch][i + myNArch] += myCstPol2;
				myHessLambda[i + myNArch][j + myNArch] = myHessLambda[j + myNArch][i + myNArch];
			}
		}
		for (uint i = 0; i < myNParam; i++)
			for (uint j = 0; j < myNParam; j++)
				mvHessPolAr[i][j] = myHessLambda[i][j];

		for (uint i = 0; i < myNParam; i++)
		{
			myGradLambda[i].Delete();
			myGradPsi[i].Delete();
			for (int j = 0; j < myNParam; j++)
				myHessLambda[i][j].Delete();
			delete[] myHessLambda[i];
		}
		delete[] myHessLambda;
		delete[] myGradPsi;
		delete[] myGradLambda;
	}

	/*!
	* \fn double cFigarch::ComputeVar(uint theDate, const cReFigarchValue& theData) const
	* \param int theDate: date of computation
	* \param const cReFigarchValue& theData: past datas
	* \details theData is not updated here.
	*/
	double cFigarch::ComputeVar(uint theDate, const cRegArchValue& theData) const
	{
		double myVar =  mvPolAr.BackwardPolOp(theData.mUt, theDate, 2);
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

	void cFigarch::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, cAbstResiduals* theResiduals)
	{
	uint myNArch = mvArch.GetSize();
	uint myNGarch = mvGarch.GetSize();
	uint myNParam = myNArch + myNGarch + 2;
	uint myBegIndex = theGradData.GetNMeanParam();
		theGradData.mCurrentGradVar = 0.0;
		for (uint i = 0; i < myNParam; i++)
			theGradData.mCurrentGradVar[myBegIndex + i] = mvGradPolAr[i].BackwardPolOp(theValue.mUt, theDate, 2);			
		for (uint n = 1; n < MIN(theDate, mvNTruncLag+1); n++)
			theGradData.mCurrentGradVar -= 2*mvPolAr[n] * theValue.mUt[theDate-n] * theGradData.mGradMt[n-1];
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

	void cFigarch::ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
		uint myNArch = mvArch.GetSize();
		uint myNGarch = mvGarch.GetSize();
		uint myNParam = myNArch + myNGarch + 2;
		uint myBegIndex = theGradData.GetNMeanParam();
		theHessData.mCurrentHessVar = 0.0;
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
	}

	void cFigarch::ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
		uint myNArch = mvArch.GetSize();
		uint myNGarch = mvGarch.GetSize();
		uint myNParam = myNArch + myNGarch + 2;
		uint myBegIndex = theGradData.GetNMeanParam();
		for (uint i = 0; i < myNParam; i++)
			for (uint j = i; j < myNParam; j++)
				theHessData.mCurrentHessVar[myBegIndex + i][myBegIndex + j] = mvHessPolAr[i][j].BackwardPolOp(theData.mUt, theDate, 2);
		for (uint k = 0; k < myNParam; k++)
			for (uint l = k; l < myNParam; l++)
				theHessData.mCurrentHessVar[myBegIndex + l][myBegIndex + k] = theHessData.mCurrentHessVar[myBegIndex + k][myBegIndex + l];
		cDVector myGradPol(theGradData.GetNParam());
		for (uint t = 1; t < MIN(theDate, mvNTruncLag + 1); t++)
		{
			for (uint i = 1; i < myNParam; i++)
				myGradPol[i + myBegIndex] = mvGradPolAr[i][t];
			theHessData.mCurrentHessVar -= 2 * theData.mUt[theDate - t] * (myGradPol * Transpose(theGradData.mGradMt[t - 1]) + theGradData.mGradMt[t - 1] * Transpose(myGradPol));
			theHessData.mCurrentHessVar += 2 * mvPolAr[t] * theGradData.mGradMt[t - 1] * Transpose(theGradData.mGradMt[t - 1]);
			theHessData.mCurrentHessVar -= 2 * mvPolAr[t] * theData.mUt[theDate - t] * theHessData.mHessMt[t - 1];
		}
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
