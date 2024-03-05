#include "StdAfxRegArchLib.h"

/*!
	\file cAparch.cpp
	\brief sources for abstract class cAparch methods.

	\author Jean-Baptiste DURAND, Ollivier TARAMASCO 
	\date dec-18-2006 - Last change feb-18-2011
*/

namespace RegArchLib {
	cAparch::cAparch(int theNArch, int theNGarch)
	:cAbstCondVar(eAparch)  // call constructor of cAbstCondVar with type eAparch
	{
		mvCste = 0.0L ;
		mvDelta = 0.0L ;
		mvArch.ReAlloc(theNArch) ;
		mvGamma.ReAlloc(theNArch) ;
		mvGarch.ReAlloc(theNGarch) ; 
		MESS_CREAT("cAparch") ;
	}

	cAparch::cAparch(const cAparch& theAparch)
	:cAbstCondVar(eAparch)  // call constructor of cAbstCondVar with type eAparch
	{
		*this = theAparch;
		MESS_CREAT("cAparch")
	}

	cAparch::~cAparch()
	{
		mvArch.Delete() ;
		mvGamma.Delete() ;
		mvGarch.Delete() ;
		MESS_DESTR("cAparch") ;
	}

	/*!
	 * \fn cAbstCondVar* cAparch::PtrCopy()
	 */
/*
cAbstCondVar* cAparch::PtrCopy() const
	{

cAparch* myAparch = new cAparch(*this);

		 return myAparch;

		return cAbstCondVarPtrCopy<cAparch>();
	}
*/
	void cAparch::Delete(void)
	{
		mvArch.Delete() ;
		mvGamma.Delete() ;
		mvGarch.Delete() ;
	}

#ifdef _RDLL_
	void cAparch::Print(void)
	{
		Rprintf("APARCH(%d, %d) model with:\n", mvArch.GetSize(), mvGarch.GetSize());
		Rprintf("Const=%f\n", mvCste);
		Rprintf("Delta=%f\n", mvDelta);
		uint i;
		for (i = 0; i < mvArch.GetSize(); i++)
			Rprintf("Arch[%d]=%f\n", i + 1, mvArch[i]);
		for (i = 0; i < mvGamma.GetSize(); i++)
			Rprintf("Gamma[%d]=%f\n", i + 1, mvGamma[i]);
		for (i = 0; i < mvGarch.GetSize(); i++)
			Rprintf("Garch[d]=%f\n", i + 1, mvGarch[i]);
	}
#else
	void cAparch::Print(ostream& theOut) const
		{
			theOut << "APARCH(" << mvArch.GetSize() << ", " << mvGarch.GetSize() << ") model:" << endl ;
			theOut << "Const=" << mvCste << endl ;
			theOut << "Delta=" << mvDelta << endl ;
		uint i ;
			for (i = 0 ; i < mvArch.GetSize() ; i++)
				theOut << "Arch[" << i+1 << "]=" << mvArch[i] << endl ;
			for (i = 0 ; i < mvGamma.GetSize() ; i++)
				theOut << "Gamma[" << i+1 << "]=" << mvGamma[i] << endl ;
			for (i = 0 ; i < mvGarch.GetSize() ; i++)
				theOut << "Garch[" << i+1 << "]=" << mvGarch[i] << endl ;
		}
#endif // _RDLL_	

	void cAparch::SetDefaultInitPoint(double theMean, double theVar)
	{
		mvCste = theVar ;
		mvDelta = 2.0 ;
	uint i ;
	for (i = 0 ; i < mvArch.GetSize() ; i++)
			mvArch[i] = 0.0 ;
		for (i = 0 ; i < mvGamma.GetSize() ; i++)
			mvGamma[i] = 0.0 ;
		for (i = 0 ; i < mvGarch.GetSize() ; i++)
			mvGarch[i] = 0.0 ;
	}

	void cAparch::SetDefaultInitPoint(cRegArchValue& theValue)
	{
		double myVar;
		theValue.ComputeVar(myVar);
		mvCste = myVar;
		mvDelta = 2.0;
		uint i;
		for (i = 0; i < mvArch.GetSize(); i++)
			mvArch[i] = 0.0;
		for (i = 0; i < mvGamma.GetSize(); i++)
			mvGamma[i] = 0.0;
		for (i = 0; i < mvGarch.GetSize(); i++)
			mvGarch[i] = 0.0;
	}

	void cAparch::ReAlloc(const uint theSize, const uint theNumParam)
	{
		switch (theNumParam)
		{	case 2 :
				if (mvArch.GetSize() != theSize)
					mvArch.ReAlloc(theSize) ;
				if (mvGamma.GetSize() != theSize)
					mvGamma.ReAlloc(theSize);
				break ;
			case 3:
				if (mvGamma.GetSize() != theSize)
					mvGamma.ReAlloc(theSize);
				if (mvArch.GetSize() != theSize)
					mvArch.ReAlloc(theSize);
				break;
			case 4 :
				if (mvGarch.GetSize() != theSize)
					mvGarch.ReAlloc(theSize);
				break ;
			default :
				break ;
		}
	}

	void cAparch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{
		switch (theNumParam)
		{	case 2 :
				mvArch = theVectParam;
			break;
			case 3:
				mvGamma = theVectParam;
			break;
			case 4 :
				mvGarch = theVectParam;
			break ;
			default:
			break;
		}
	}

	void cAparch::Set(const cDVector& theDVector, const uint thePlace) 
	{
		uint mySize = 0;
		switch (thePlace)
		{	case 0 :
				mvCste = theDVector[0] ;
			break ;
			case 1 :
				mvDelta = theDVector[0] ;
			break ;
			case 2:
				mvArch = theDVector ;
				mySize = mvArch.GetSize();
				if (mvGamma.GetSize() != mySize)
					mvGamma.ReAlloc(mySize);
				break;
			case 3:
				mvGamma = theDVector;
				mySize = mvGamma.GetSize();
				if (mvArch.GetSize() != mySize)
					mvArch.ReAlloc(mySize);
				break;
			case 4 :
				mvGarch = theDVector;
			break ;
			default:
				break;
		}
	}

	void cAparch::Set(const double theValue, const uint theIndex, const uint theNumParam) 
	{
	uint mySize ;
		switch (theNumParam)
		{	
			case 0 :
				mvCste = theValue ;
			break ;
			case 1 :
				mvDelta = theValue ;
			break ;
			case 2 :
			{	mySize = mvArch.GetSize() ;
				if ( mySize < theIndex)
				{
					mvArch.ReAlloc(theIndex+1) ;
					mvGamma.ReAlloc(theIndex+1) ;
				}
				mvArch[theIndex] = theValue ;
			}
			break;
			case 3 :
			{	mySize = mvGamma.GetSize() ;
				if ( mySize < theIndex)
				{
					mvArch.ReAlloc(theIndex+1) ;
					mvGamma.ReAlloc(theIndex+1) ;
				}
				mvGamma[theIndex] = theValue ;
			}
			break;
			case 4 :
			{	mySize = mvGarch.GetSize() ;
				if ( mySize < theIndex)
					mvGarch.ReAlloc(theIndex+1) ;
				mvGarch[theIndex] = theValue ;	
			}
			break ;
			default :
			break ;
		}
	}

	double cAparch::Get(const uint theIndex, const uint theNumParam)
	{
		switch (theNumParam)
		{	case 0 :
				return mvCste ;
			break ;
			case 1 :
				return mvDelta ;
			break ;
			case 2 :
				return mvArch[theIndex] ;
			break ;
			case 3 :
				return mvGamma[theIndex] ;
			break ;
			case 4 :
				return mvGarch[theIndex] ;
			break ;
			default:
				break;
		}
	}

	cDVector& cAparch::Get(const uint theNumParam)
	{
	cDVector* myAux;
		switch (theNumParam)
		{
			case 0:
				myAux = new cDVector(1, mvCste);
				return *myAux;
				break;
			case 1:
				myAux = new cDVector(1, mvDelta);
				return *myAux;
				break;
			case 2:
				return mvArch;
				break;
			case 3:
				return mvGamma;
				break;
			case 4:
				return mvGarch;
				break;
			default:
				break;
		}
	}

	cAparch& cAparch::operator =(const cAparch& theSrc)
	{
		mvArch = theSrc.mvArch;
		mvGarch = theSrc.mvGarch;
		mvCste = theSrc.mvCste;
		mvGamma = theSrc.mvGamma;
		mvDelta = theSrc.mvDelta;
		return *this ;
	}

	double cAparch::ComputeVar(uint theDate, const cRegArchValue& theValue) const
	{
	uint myp = mvArch.GetSize() ;
	uint myq = mvGarch.GetSize() ;

	double	myRes = mvCste,
			myTemp = 0.0,
			myDeltaDiv2 = mvDelta/2.0 ;
	
		for (uint i = 1 ; i <= MIN(myp, theDate) ; i++)
		{	myTemp = abs(theValue.mUt[theDate-i]) - mvGamma[i-1] * theValue.mUt[theDate-i];
			myRes += mvArch[i-1] * pow(myTemp, mvDelta) ;
		}

		for (uint i = 1 ; i <= MIN(myq, theDate) ; i++)
			myRes += mvGarch[i-1] * pow(theValue.mHt[theDate-i], myDeltaDiv2) ;

		return pow(myRes, 1.0/myDeltaDiv2) ;
	}

	uint cAparch::GetNParam(void) const
	{
		return (2*(int)mvArch.GetSize() + (int)mvGarch.GetSize() + 2) ;
	}

	uint cAparch::GetNLags(void) const
	{
		return MAX(mvArch.GetSize(), mvGarch.GetSize()) ;
	}

	void cAparch::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{
		uint mySize = GetNParam();
		uint myIndex = theIndex;
	
		if (theDestVect.GetSize() < mySize + theIndex)
			throw cError("Wrong size") ;
		theDestVect[myIndex++] = mvCste ;
		theDestVect[myIndex++] = mvDelta ;
		mvArch.SetSubVectorWithThis(theDestVect, myIndex) ;
		myIndex += mvArch.GetSize() ;
		mvGamma.SetSubVectorWithThis(theDestVect, myIndex) ;
		myIndex += mvGamma.GetSize() ;
		mvGarch.SetSubVectorWithThis(theDestVect, myIndex) ;
	}

	void cAparch::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	uint mySize = theSrcVect.GetSize(),
		myIndex = theIndex ;
		if (GetNParam() + theIndex > mySize)
			throw cError("Wrong size") ;
		mvCste = theSrcVect[myIndex++] ;
		mvDelta = theSrcVect[myIndex++] ;
		mvArch.SetThisWithSubVector(theSrcVect,myIndex) ;
		myIndex += mvArch.GetSize() ;
		mvGamma.SetThisWithSubVector(theSrcVect,myIndex) ;
		myIndex += mvGamma.GetSize() ;
		mvGarch.SetThisWithSubVector(theSrcVect,myIndex) ;
	}

	void cAparch::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, cAbstResiduals* theResiduals)
	{
	uint myp = mvArch.GetSize();
	uint myq = mvGarch.GetSize();
	uint myNParam = theGradData.GetNParam();
	uint myNMeanParam = theGradData.GetNMeanParam();
	uint myNVarParam = theGradData.GetNMeanParam();
	uint myBegIndex = myNMeanParam;

	double myDeltaDiv2 = mvDelta / 2.0;
	cDVector* myDZ = new cDVector[myp];
	cDVector myZ(myNParam);
	cDVector myVt(myNParam);
	cDVector* myDVt = new cDVector[myq];
	for (uint i = 0; i < MIN(theDate, myp); i++)
	{
		myZ[i] = fabs(theValue.mUt[theDate - i - 1]) - mvGamma[i] * theValue.mUt[theDate - i - 1];
		if (theValue.mUt[theDate - i - 1] > 0)
			myDZ[i] = (mvGamma[i] - 1.0) * theGradData.mGradMt[i];
		else
			myDZ[i] = (mvGamma[i] + 1) * theGradData.mGradMt[i];
		myDZ[i][myBegIndex + i + myp + 2] -= theValue.mUt[theDate - i - 1];
	}
	for (uint j = 0; j < MIN(myq, theDate); j++)
	{
		myVt[j] = pow(theValue.mHt[theDate - j - 1], myDeltaDiv2);
		myDVt[j] = myDeltaDiv2 * pow(theValue.mHt[theDate - j - 1], myDeltaDiv2 - 1) * theGradData.mGradHt[j];
		myDVt[j][myBegIndex + 1] += 0.5 * log(theValue.mHt[theDate - j - 1]) * myVt[j];
	}

	cDVector myRes = theGradData.mCurrentGradVar = 0.0;
	myRes[myBegIndex] = 1.0;
	for (uint i = 0; i < MIN(myp, theDate); i++)
	{
		myRes += mvArch[i] * mvDelta * pow(myZ[i], mvDelta - 1.0) * myDZ[i];
		myRes[myBegIndex + i + 2] += pow(myZ[i], mvDelta);
		myRes[myBegIndex + 1] += mvArch[i] * log(myZ[i]) * pow(myZ[i], mvDelta);
	}
	for (uint j = 0; j < MIN(myq, theDate); j++)
	{
		myRes += (mvGarch[j]) * myDVt[j];
		myRes[myBegIndex + 2 * myp + 2 + j] += myVt[j];
	}
	double myAux = 2.0 / (mvDelta * pow(theValue.mHt[theDate], myDeltaDiv2 - 1.0));
	theGradData.mCurrentGradVar = myAux * myRes;
	theGradData.mCurrentGradVar[myBegIndex + 1] -= 0.5 * log(theValue.mHt[theDate]) * pow(theValue.mHt[theDate], myDeltaDiv2) * myAux;
	for (uint i = 0; i < myp; i++)
	{
		myDZ[i].Delete();
		myDVt[i].Delete();
	}
	delete[] myDZ;
	delete[] myDVt;
/*
	cDVector* myDArch = new cDVector[myp];
	cDVector* myDGarch = new cDVector[myq];
	cDVector* myDGamma = new cDVector[myp];
		for (uint i = 0; i < myp; i++)
		{
			myDArch[i].ReAlloc(myNParam, 0.0);
			myDArch[i][myBegIndex + i + 2] = 1.0;
			myDGamma[i].ReAlloc(myNParam, 0.0);
			myDGamma[i][myBegIndex + myp + i + 2] = 1.0;
		}
		for (uint j = 0; j < myq; j++)
		{
			myDGarch[j].ReAlloc(myNParam, 0.0);
			myDGarch[j][myBegIndex + 2 * myp + 2 + j] = 1.0;
		}
	cDVector myDDelta(myNParam, 0.0);
		myDDelta[myBegIndex + 1] = 1.0;
	cDVector myDOmega(myNParam, 0.0);
		myDOmega[myBegIndex] = 1.0;
	cDVector myGradNu = myDOmega;
		for (uint i = 0; i < MIN(myp, theDate); i++)
		{
		double myU = theData.mUt[theDate - i - 1];
		double mySign = SIGN(myU);
		double myCrochet = myU * (mySign - mvGamma[i]);
		double myCrochetpDelta = pow(myCrochet, mvDelta);
		double myCrochetpDeltam1 = pow(myCrochet, mvDelta - 1);
		double myCrochetpDeltam2 = pow(myCrochet, mvDelta - 2);
		double myLogCrochet = log(myCrochet);
		cDVector myDCrochet = -myU * myDGamma[i] - (mySign - mvGamma[i])*theGradData.mGradMt[i];
			myGradNu += myDArch[i] * myCrochetpDelta;
			myGradNu += mvArch[i] * mvDelta * myDCrochet * myCrochetpDeltam1;
			myGradNu += mvArch[i] * myLogCrochet * myCrochetpDelta * myDDelta;
		}
		for (uint j = 0; j < MIN(myq, theDate); j++)
		{
		double	myNu = pow(theData.mHt[theDate - j - 1], mvDelta/2);
		double	myNupm1 = pow(theData.mHt[theDate - j - 1], mvDelta/2 - 1.0);
		double	myNupm2 = pow(theData.mHt[theDate - j - 1], mvDelta/2 - 2.0);
		double myLogh = log(theData.mHt[theDate - j - 1]);
		cDVector myDNu = mvDelta/2 * myNupm1 * theGradData.mGradHt[j] + myLogh/2 * myNu * myDDelta;
			myGradNu += myDGarch[j] * myNu + mvGarch[j] * myDNu;

		}
		// On a caluler Grad nu = Grad h(t)^(delta/2), on doit calculer Grad h(t)
	double	myNu = pow(theData.mHt[theDate], mvDelta / 2);
	double	myNupm1 = pow(theData.mHt[theDate], mvDelta / 2 - 1.0);
	double myLogh = log(theData.mHt[theDate]);
		theGradData.mCurrentGradVar = myGradNu - myLogh / 2 * myNu * myDDelta;
		theGradData.mCurrentGradVar /= mvDelta / 2 * myNupm1;

		for (uint i = 0; i < myp; i++)
		{
			myDArch[i].Delete();
			myDGamma[i].Delete();
		}
		for (uint j = 0; j < myq; j++)
		{
			myDGarch[j].Delete();
		}
		delete[] myDArch;
		delete[] myDGarch;
		delete[] myDGamma;

*/
	}

	void cAparch::ComputeHess(uint theDate, const cRegArchValue& theValue, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
		uint myp = mvArch.GetSize();
		uint myq = mvGarch.GetSize();
		uint myNParam = theGradData.GetNParam();
		uint myNMeanParam = theGradData.GetNMeanParam();
		uint myNVarParam = theGradData.GetNMeanParam();
		uint myBegIndex = myNMeanParam;

		double myDeltaDiv2 = mvDelta / 2.0;
		cDVector* myDZ = new cDVector[myp];
		cDVector myZ(myNParam);
		cDVector myVt(myNParam);
		cDVector* myDVt = new cDVector[myq];
		for (uint i = 0; i < MIN(theDate, myp); i++)
		{
			myZ[i] = fabs(theValue.mUt[theDate - i - 1]) - mvGamma[i] * theValue.mUt[theDate - i - 1];
			if (theValue.mUt[theDate - i - 1] > 0)
				myDZ[i] = (mvGamma[i] - 1.0) * theGradData.mGradMt[i];
			else
				myDZ[i] = (mvGamma[i] + 1) * theGradData.mGradMt[i];
			myDZ[i][myBegIndex + i + myp + 2] -= theValue.mUt[theDate - i - 1];
		}
		for (uint j = 0; j < MIN(myq, theDate); j++)
		{
			myVt[j] = pow(theValue.mHt[theDate - j - 1], myDeltaDiv2);
			myDVt[j] = myDeltaDiv2 * pow(theValue.mHt[theDate - j - 1], myDeltaDiv2 - 1) * theGradData.mGradHt[j];
			myDVt[j][myBegIndex + 1] += 0.5 * log(theValue.mHt[theDate - j - 1]) * myVt[j];
		}
		cDVector myX(myNParam);
		cDVector* myDX = new cDVector[myp];
		for (uint i = 0; i < MIN(theDate, myp); i++)
		{
			myX[i] = pow(myZ[i], mvDelta);
			myDX[i] = mvDelta * myX[i] / myZ[i] * myDZ[i];
			myDX[i][myBegIndex + 1] += log(myZ[i]) * myX[i];
		}

		cDMatrix* myD2Z = new cDMatrix[myp];
		cDMatrix* myD2X = new cDMatrix[myp];
		cDMatrix* myD2Vt = new cDMatrix[myp];
		for (uint i = 0; i < MIN(myp, theDate); i++)
		{
			if (theValue.mUt[theDate - i - 1] > 0)
				myD2Z[i] = (mvGamma[i] - 1) * theHessData.mHessMt[i];
			else
				myD2Z[i] = (mvGamma[i] + 1) * theHessData.mHessMt[i];
			for (uint j = 0; j < myNParam; j++)
			{
				myD2Z[i][myBegIndex + myp + 2 + i][j] += theGradData.mGradMt[i][j];
				myD2Z[i][j][myBegIndex + myp + 2 + i] += theGradData.mGradMt[i][j];
			}
			myD2X[i] = mvDelta / myZ[i] * myX[i] * myD2Z[i];
			myD2X[i] -= mvDelta / pow(myZ[i], 2) * myX[i] * myDZ[i] * Transpose(myDZ[i]);
			for (uint k = 0; k < myNParam; k++)
			{
				myD2X[i][myBegIndex + 1][k] += myX[i] / myZ[i] * myDZ[i][k];
				myD2X[i][k][myBegIndex + 1] += myX[i] / myZ[i] * myDZ[i][k];
			}
			myD2X[i] += 1.0 / myX[i] * myDX[i] * Transpose(myDX[i]);
		}
		double myV = pow(theValue.mHt[theDate], myDeltaDiv2);
		for (uint j = 0; j < MIN(myq, theDate); j++)
		{
			myD2Vt[j] = myDeltaDiv2 * myVt[j] / theValue.mHt[theDate - j - 1] * theHessData.mHessHt[j];
			myD2Vt[j] -= myDeltaDiv2 / pow(theValue.mHt[theDate - j - 1], 2.0) * myVt[j] * theGradData.mGradHt[j] * Transpose(theGradData.mGradHt[j]);
			for (uint k = 0; k < myNParam; k++)
			{
				myD2Vt[j][myBegIndex + 1][k] += 0.5 * myVt[j] / theValue.mHt[theDate - j - 1] * theGradData.mGradHt[j][k];
				myD2Vt[j][k][myBegIndex + 1] += 0.5 * myVt[j] / theValue.mHt[theDate - j - 1] * theGradData.mGradHt[j][k];
			}
			myD2Vt[j] += 1.0 / myVt[j] * myDVt[j] * Transpose(myDVt[j]);
		}
		cDMatrix myMat = theHessData.mCurrentHessVar = 0.0;
		for (uint i = 0; i < MIN(myp, theDate); i++)
		{
			myMat += mvArch[i] * myD2X[i];
			for (uint k = 0; k < myNParam; k++)
			{
				myMat[myBegIndex + i + 2][k] += myDX[i][k];
				myMat[k][myBegIndex + i + 2] += myDX[i][k];
			}
		}
		for (uint j = 0; j < MIN(theDate, myq); j++)
		{
			myMat += mvGarch[j] * myD2Vt[j];
			for (uint k = 0; k < myNParam; k++)
			{
				myMat[myBegIndex + 2 * myp + 2 + j][k] += myDVt[j][k];
				myMat[k][myBegIndex + 2 * myp + 2 + j] += myDVt[j][k];
			}
		}
		cDVector myDV = myDeltaDiv2 / theValue.mHt[theDate] * theGradData.mCurrentGradVar;
		myDV[myBegIndex + 1] += log(theValue.mHt[theDate]) / 2.0;
		myDV *= myV;
		theHessData.mCurrentHessVar = 2.0 / (mvDelta * myV) * myMat;
		theHessData.mCurrentHessVar -= 1.0 / (myDeltaDiv2 * pow(myV, 2)) * myDV * Transpose(myDV);
		for (uint k = 0; k < myNParam; k++)
		{
			theHessData.mCurrentHessVar[k][myBegIndex + 1] -= 2.0 / (pow(mvDelta, 2) * myV) * myDV[k];
			theHessData.mCurrentHessVar[myBegIndex + 1][k] -= 2.0 / (pow(mvDelta, 2) * myV) * myDV[k];
		}
		theHessData.mCurrentHessVar[myBegIndex + 1][myBegIndex + 1] += 4.0 / pow(mvDelta, 3) * log(myV);
		theHessData.mCurrentHessVar *= theValue.mHt[theDate];
		cDVector myVect = 1.0 / (myDeltaDiv2 * myV) * myDV;
		myVect[myBegIndex + 1] -= 2.0 / pow(mvDelta, 2) * log(myV);
		theHessData.mCurrentHessVar += myVect * Transpose(theGradData.mCurrentGradVar);
		for (uint i = 0; i < myp; i++)
		{
			myDZ[i].Delete();
			myDX[i].Delete();
			myD2Z[i].Delete();
			myD2X[i].Delete();
		}
		for (uint j = 0; j < myq; j++)
		{
			myDVt[j].Delete();
			myD2Vt[j].Delete();
		}
		delete[] myDZ;
		delete[] myDX;
		delete[] myDVt;
		delete[] myD2Z;
		delete[] myD2X;
		delete[] myD2Vt;

	}

	void cAparch::ComputeGradAndHess(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
		uint myp = mvArch.GetSize();
		uint myq = mvGarch.GetSize();
		uint myNParam = theGradData.GetNParam();
		uint myNMeanParam = theGradData.GetNMeanParam();
		uint myNVarParam = theGradData.GetNMeanParam();
		uint myBegIndex = myNMeanParam;

		double myDeltaDiv2 = mvDelta / 2.0;
		cDVector* myDZ = new cDVector[myp];
		cDVector myZ(myNParam);
		cDVector myVt(myNParam);
		cDVector* myDVt = new cDVector[myq];
		for (uint i = 0; i < MIN(theDate, myp); i++)
		{
			myZ[i] = fabs(theValue.mUt[theDate - i - 1]) - mvGamma[i] * theValue.mUt[theDate - i - 1];
			if (theValue.mUt[theDate - i - 1] > 0)
				myDZ[i] = (mvGamma[i] - 1.0) * theGradData.mGradMt[i];
			else
				myDZ[i] = (mvGamma[i] + 1) * theGradData.mGradMt[i];
			myDZ[i][myBegIndex + i + myp + 2] -= theValue.mUt[theDate - i - 1];
		}
		for (uint j = 0; j < MIN(myq, theDate); j++)
		{
			myVt[j] = pow(theValue.mHt[theDate - j - 1], myDeltaDiv2);
			myDVt[j] = myDeltaDiv2 * pow(theValue.mHt[theDate - j - 1], myDeltaDiv2 - 1) * theGradData.mGradHt[j];
			myDVt[j][myBegIndex + 1] += 0.5 * log(theValue.mHt[theDate - j - 1]) * myVt[j];
		}

		cDVector myRes = theGradData.mCurrentGradVar = 0.0;
		myRes[myBegIndex] = 1.0;
		for (uint i = 0; i < MIN(myp, theDate); i++)
		{
			myRes += mvArch[i] * mvDelta * pow(myZ[i], mvDelta - 1.0) * myDZ[i];
			myRes[myBegIndex + i + 2] += pow(myZ[i], mvDelta);
			myRes[myBegIndex + 1] += mvArch[i] * log(myZ[i]) * pow(myZ[i], mvDelta);
		}
		for (uint j = 0; j < MIN(myq, theDate); j++)
		{
			myRes += (mvGarch[j]) * myDVt[j];
			myRes[myBegIndex + 2 * myp + 2 + j] += myVt[j];
		}
		double myAux = 2.0 / (mvDelta * pow(theValue.mHt[theDate], myDeltaDiv2 - 1.0));
		theGradData.mCurrentGradVar = myAux * myRes;
		theGradData.mCurrentGradVar[myBegIndex + 1] -= 0.5 * log(theValue.mHt[theDate]) * pow(theValue.mHt[theDate], myDeltaDiv2) * myAux;

		cDVector myX(myNParam);
		cDVector* myDX = new cDVector[myp];
		for (uint i = 0; i < MIN(theDate, myp); i++)
		{
			myX[i] = pow(myZ[i], mvDelta);
			myDX[i] = mvDelta * myX[i] / myZ[i] * myDZ[i];
			myDX[i][myBegIndex + 1] += log(myZ[i]) * myX[i];
		}

		cDMatrix* myD2Z = new cDMatrix[myp];
		cDMatrix* myD2X = new cDMatrix[myp];
		cDMatrix* myD2Vt = new cDMatrix[myp];
		for (uint i = 0; i < MIN(myp, theDate); i++)
		{
			if (theValue.mUt[theDate - i - 1] > 0)
				myD2Z[i] = (mvGamma[i] - 1) * theHessData.mHessMt[i];
			else
				myD2Z[i] = (mvGamma[i] + 1) * theHessData.mHessMt[i];
			for (uint j = 0; j < myNParam; j++)
			{
				myD2Z[i][myBegIndex + myp + 2 + i][j] += theGradData.mGradMt[i][j];
				myD2Z[i][j][myBegIndex + myp + 2 + i] += theGradData.mGradMt[i][j];
			}
			myD2X[i] = mvDelta / myZ[i] * myX[i] * myD2Z[i];
			myD2X[i] -= mvDelta / pow(myZ[i], 2) * myX[i] * myDZ[i] * Transpose(myDZ[i]);
			for (uint k = 0; k < myNParam; k++)
			{
				myD2X[i][myBegIndex + 1][k] += myX[i] / myZ[i] * myDZ[i][k];
				myD2X[i][k][myBegIndex + 1] += myX[i] / myZ[i] * myDZ[i][k];
			}
			myD2X[i] += 1.0 / myX[i] * myDX[i] * Transpose(myDX[i]);
		}
		double myV = pow(theValue.mHt[theDate], myDeltaDiv2);
		for (uint j = 0; j < MIN(myq, theDate); j++)
		{
			myD2Vt[j] = myDeltaDiv2 * myVt[j] / theValue.mHt[theDate - j - 1] * theHessData.mHessHt[j];
			myD2Vt[j] -= myDeltaDiv2 / pow(theValue.mHt[theDate - j - 1], 2.0) * myVt[j] * theGradData.mGradHt[j] * Transpose(theGradData.mGradHt[j]);
			for (uint k = 0; k < myNParam; k++)
			{
				myD2Vt[j][myBegIndex + 1][k] += 0.5 * myVt[j] / theValue.mHt[theDate - j - 1] * theGradData.mGradHt[j][k];
				myD2Vt[j][k][myBegIndex + 1] += 0.5 * myVt[j] / theValue.mHt[theDate - j - 1] * theGradData.mGradHt[j][k];
			}
			myD2Vt[j] += 1.0 / myVt[j] * myDVt[j] * Transpose(myDVt[j]);
		}
		cDMatrix myMat = theHessData.mCurrentHessVar = 0.0;
		for (uint i = 0; i < MIN(myp, theDate); i++)
		{
			myMat += mvArch[i] * myD2X[i];
			for (uint k = 0; k < myNParam; k++)
			{
				myMat[myBegIndex + i + 2][k] += myDX[i][k];
				myMat[k][myBegIndex + i + 2] += myDX[i][k];
			}
		}
		for (uint j = 0; j < MIN(theDate, myq); j++)
		{
			myMat += mvGarch[j] * myD2Vt[j];
			for (uint k = 0; k < myNParam; k++)
			{
				myMat[myBegIndex + 2 * myp + 2 + j][k] += myDVt[j][k];
				myMat[k][myBegIndex + 2 * myp + 2 + j] += myDVt[j][k];
			}
		}
		cDVector myDV = myDeltaDiv2 / theValue.mHt[theDate] * theGradData.mCurrentGradVar;
		myDV[myBegIndex + 1] += log(theValue.mHt[theDate]) / 2.0;
		myDV *=  myV;
		theHessData.mCurrentHessVar = 2.0 / (mvDelta * myV) * myMat;
		theHessData.mCurrentHessVar -= 1.0 / (myDeltaDiv2 * pow(myV, 2)) * myDV * Transpose(myDV);
		for (uint k = 0; k < myNParam; k++)
		{
			theHessData.mCurrentHessVar[k][myBegIndex + 1] -= 2.0 / (pow(mvDelta, 2) * myV) * myDV[k];
			theHessData.mCurrentHessVar[myBegIndex + 1][k] -= 2.0 / (pow(mvDelta, 2) * myV) * myDV[k];
		}
		theHessData.mCurrentHessVar[myBegIndex + 1][myBegIndex + 1] += 4.0 / pow(mvDelta, 3) * log(myV);
		theHessData.mCurrentHessVar *= theValue.mHt[theDate];
		cDVector myVect = 1.0 / (myDeltaDiv2 * myV) * myDV;
		myVect[myBegIndex + 1] -= 2.0 / pow(mvDelta, 2) * log(myV);
		theHessData.mCurrentHessVar += myVect * Transpose(theGradData.mCurrentGradVar);
		for (uint i = 0; i < myp; i++)
		{
			myDZ[i].Delete();
			myDX[i].Delete();
			myD2Z[i].Delete();
			myD2X[i].Delete();
		}
		for (uint j = 0; j < myq; j++)
		{
			myDVt[j].Delete();
			myD2Vt[j].Delete();
		}
		delete[] myDZ;
		delete[] myDX;
		delete[] myDVt;
		delete[] myD2Z;
		delete[] myD2X;
		delete[] myD2Vt;

	}

	void cAparch::GetParamName(uint theIndex, char** theName)
	{
	uint myIndex = theIndex;
		sprintf(theName[myIndex++], "CST VAR");
		sprintf(theName[myIndex +++1], "DELTA");
		for (uint i = 0; i < mvArch.GetSize(); i++)
		{
			sprintf(theName[myIndex++], "ARCH[%d]", i + 1);

		}
		for (uint i = 0; i < mvGamma.GetSize(); i++)
		{
			sprintf(theName[myIndex++], "GAMMA[%d]", i + 1);

		}
		for (uint i = 0; i < mvGarch.GetSize(); i++)
		{
			sprintf(theName[myIndex++], "GARCH[%d]", i + 1);

		}

	}

	void cAparch::GetParamName(uint theIndex, string theName[])
	{
	uint myIndex = theIndex;
	char myChar[100];

		sprintf(myChar, "CST VAR");
		theName[myIndex++] = myChar;
		sprintf(myChar, "DELTA");
		theName[myIndex++] = myChar;
		for (uint i = 0; i < mvArch.GetSize(); i++)
		{
			sprintf(myChar, "ARCH[%d]", i + 1);
			theName[myIndex++] = myChar;

		}
		for (uint i = 0; i < mvGamma.GetSize(); i++)
		{
			sprintf(myChar, "GAMMA[%d]", i + 1);
			theName[myIndex++] = myChar;

		}
		for (uint i = 0; i < mvGarch.GetSize(); i++)
		{
			sprintf(myChar, "GARCH[%d]", i + 1);
			theName[myIndex++] = myChar;

		}

	}

}//namespace
