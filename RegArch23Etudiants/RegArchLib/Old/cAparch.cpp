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
		{	case 0 :
			case 1 :
				mvArch.ReAlloc(theSize) ;
				mvGamma.ReAlloc(theSize) ;
			break ;
			case 2 :
				mvGarch.ReAlloc(theSize) ;
			break ;
			default :
			break ;
		}
	}

	void cAparch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
	{
		switch (theNumParam)
		{	case 0 :
			case 1 :
				mvArch.ReAlloc((int)theVectParam[0]) ;
				mvGamma.ReAlloc((int)theVectParam[0]) ;
			break ;
			case 2 :
				mvGarch.ReAlloc((int)theVectParam[0]) ;
			break ;
		}
	}

	void cAparch::Set(const cDVector& theDVector, const uint thePlace) 
	{
		switch (thePlace)
		{	case 0 :
				mvCste = theDVector[0] ;
			break ;
			case 1 :
				mvDelta = theDVector[0] ;
			break ;
			case 2:
				mvArch = theDVector ;
			case 3:
				mvGamma = theDVector;
			case 4 :
				mvGarch = theDVector;
			break ;
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
				{	mvArch.ReAlloc(theIndex+1) ;
					mvGamma.ReAlloc(theIndex+1) ;
				}
				mvArch[theIndex] = theValue ;
			}
			case 3 :
			{	mySize = mvGamma.GetSize() ;
				if ( mySize < theIndex)
				{	mvArch.ReAlloc(theIndex+1) ;
					mvGamma.ReAlloc(theIndex+1) ;
				}
				mvGamma[theIndex] = theValue ;
			}
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

	double  cAparch::Get(const uint theIndex, const uint theNumParam)
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

	uint cAparch::GetNu(void) const
	{
		return mvArch.GetSize();
	}
	
	uint cAparch::GetNh(void) const
	{
		return mvGarch.GetSize();
	}

	void cAparch::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
	{
	uint	mySize = GetNParam(),
			myIndex = theIndex		;
	
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
	{uint	mySize = theSrcVect.GetSize(),
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

	/*
	void cAparch::ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cAbstResiduals* theResiduals)
	{
	uint myp = mvArch.GetSize(),
		myq = mvGarch.GetSize(),
		myBegIndex = theGradData.GetNMeanParam() 	;
		theGradData.mCurrentGradVar = 0.0 ;
	double myDeltaDiv2 = mvDelta/2.0 ;
	
	// dérivée par rapport à la constante 
		theGradData.mCurrentGradVar[myBegIndex++] = 1.0 ;
	uint i ;

	// dérivée par rapport à Delta 
		theGradData.mCurrentGradVar[myBegIndex] = 0.0 ;
		for (i = 1 ; i <= MIN(myp, theDate) ; i++)
		{	
		double	myTemp1 = abs(theData.mUt[theDate-i]) - mvGamma[i-1]*theData.mUt[theDate-i],
				myTemp2 = pow(myTemp1, mvDelta) ;
			theGradData.mCurrentGradVar[myBegIndex] += mvArch[i-1]* log(myTemp1) * myTemp2 ;
		}
		for (i = 1 ; i <= MIN(myq, theDate) ; i++)
		{	
		double	myTemp1 = pow(theData.mHt[theDate-i], myDeltaDiv2) ;
			theGradData.mCurrentGradVar[myBegIndex] += mvGarch[i-1] * myTemp1 * log(theData.mHt[theDate-i])/2.0 ;
		}

	// d�riv�e par rapport � Arch 
		for (i = 1 ; i <= MIN(myp, theDate) ; i++)
		{	
		double	myTemp1 = abs(theData.mUt[theDate-i]) - mvGamma[i-1]*theData.mUt[theDate-i],
				myTemp2 = pow(myTemp1, mvDelta) ;
			theGradData.mCurrentGradVar[myBegIndex+i] = myTemp2 ;
		}

	// d�riv�e par rapport � Gamma 
		myBegIndex += myp ;
		for (i = 1 ; i <= MIN(myp, theDate) ; i++)
		{	
		double	myTemp1 = abs(theData.mUt[theDate-i]) - mvGamma[i-1]*theData.mUt[theDate-i],
				myTemp2 = pow(myTemp1, mvDelta) ;
			theGradData.mCurrentGradVar[myBegIndex+i] = -mvArch[i-1] * mvDelta * theData.mUt[theDate-i] * myTemp2/myTemp1 ;
		}

	// d�riv�e par rapport � Garch
		myBegIndex += myp ;
		for (i = 1 ; i <= MIN(myq, theDate) ; i++)
			theGradData.mCurrentGradVar[myBegIndex+i] = pow(theData.mHt[theDate-i], myDeltaDiv2) ;
	
	// Il faut maintenant rajouter les d�riv�es par rapport aux u(t-i) et h(t-j) 
		for ( i = 1 ; i <= MIN(myp, theDate) ; i++)
		{
		double myTemp1 = abs(theData.mUt[theDate-i])/theData.mUt[theDate-i] - mvGamma[i-1] ;

			theGradData.mCurrentGradVar -= mvArch[i-1] * mvDelta * myTemp1 * theGradData.mGradMt[i-1] ; 
		}
	
		for ( i = 1 ; i <= MIN(myq, theDate) ; i++)
		{
		double myTemp1 = mvGarch[i-1] * myDeltaDiv2 * pow(theData.mHt[theDate-i], myDeltaDiv2 - 1.0) ;
			theGradData.mCurrentGradVar += myTemp1 * theGradData.mGradHt[i-1] ; 
		}

		// on a calcul� d(ht^Delta/2)/dTeta ... Il faut calculer dht/dTeta 
		myBegIndex = (int)theGradData.GetNMeanParam() + 1 ; // On recup�re l'indice de Delta
	double myAux = theGradData.mCurrentGradVar[myBegIndex] ;
		// Formule valable pour tous les param�tres sauf Delta 
	double myTemp = pow(theData.mHt[theDate], myDeltaDiv2-1) ;
		theGradData.mCurrentGradVar /= myDeltaDiv2 * myTemp  ;
		// Et pour Delta 
		theGradData.mCurrentGradVar[myBegIndex] = (myAux - 0.5*log(theData.mHt[theDate]) * myTemp * theData.mHt[theDate])/(myDeltaDiv2 * myTemp) ;
	}
	*/
	void cAparch::ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cAbstResiduals* theResiduals)
	{
	uint myp = mvArch.GetSize();
	uint myq = mvGarch.GetSize();
	uint myNParam = theGradData.GetNParam();
	uint myNMeanParam = theGradData.GetNMeanParam();
	uint myNVarParam = theGradData.GetNMeanParam();
	uint myBegIndex = myNMeanParam;

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
	}

	void cAparch::ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cFuncMeanAndVar& theDerivM, cAbstResiduals* theResids)
	{
	}
	
	void cAparch::ComputeHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	uint myp = mvArch.GetSize();
	uint myq = mvGarch.GetSize();
	uint myNParam = theGradData.GetNParam();
	uint myNDistrParam = theGradData.GetNDistrParameter();
	uint myNMeanParam = theGradData.GetNMeanParam();
	uint myNVarParam = theGradData.GetNMeanParam();
	uint myBegIndex = myNMeanParam;

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
			myDGarch[j][myBegIndex + 2*myp + 2 + j] = 1.0;
		}
	cDVector myDDelta(myNParam, 0.0);
		myDDelta[myBegIndex + 1]= 1.0;	
	cDMatrix myHessNu(myNParam, myNParam);
		for (uint i = 0; i < MIN(myp, theDate); i++)
		{
		double myU = theData.mUt[theDate - i - 1];
		double mySign = SIGN(myU);
		double myCrochet = myU * (mySign - mvGamma[i]);
		double myCrochetpDelta = pow(myCrochet, mvDelta);
		double myCrochetpDeltam1 = pow(myCrochet, mvDelta-1.0);
		double myCrochetpDeltam2 = pow(myCrochet, mvDelta - 2.0);
		double myLogCrochet = log(myCrochet);
		cDVector myDCrochet = -myU * myDGamma[i] - (mySign - mvGamma[i])*theGradData.mGradMt[i];
			myHessNu += myLogCrochet * myCrochetpDelta*(myDArch[i] * Transpose(myDDelta) + myDDelta * Transpose(myDArch[i]));
			myHessNu += mvDelta * myCrochetpDeltam1 * (myDArch[i] * Transpose(myDCrochet) + myDCrochet * Transpose(myDArch[i]));
			myHessNu += mvArch[i] * myCrochetpDeltam1 * (myDDelta * Transpose(myDCrochet) + myDCrochet * Transpose(myDDelta)); 
			myHessNu += mvArch[i] * myLogCrochet * mvDelta * myCrochetpDeltam1 * (myDDelta * Transpose(myDCrochet) + myDCrochet * Transpose(myDDelta));
			myHessNu += mvArch[i] * mvDelta * (mvDelta - 1) * myCrochetpDeltam2 * myDCrochet * Transpose(myDCrochet);
			myHessNu += mvArch[i] * pow(myLogCrochet, 2) * myCrochetpDelta * myDDelta * Transpose(myDDelta);
		cDMatrix myD2Crochet = -(mySign - mvGamma[i]) * theHessData.mHessMt[i] + myDGamma[i] * Transpose(theGradData.mGradMt[i]) + theGradData.mGradMt[i] * Transpose(myDGamma[i]);
			myHessNu += mvArch[i] * mvDelta * myCrochetpDeltam1 * myD2Crochet;
		}
		double myDeltaDiv2 = mvDelta/2.0;
		for (uint j = 0; j < MIN(myq, theDate); j++)
		{
		double	myNu = pow(theData.mHt[theDate - j - 1], myDeltaDiv2);
		double	myNupm1 = pow(theData.mHt[theDate - j - 1], myDeltaDiv2 -1.0);
		double	myNupm2 = pow(theData.mHt[theDate - j - 1], myDeltaDiv2 - 2.0);
		double myLogh = log(theData.mHt[theDate - j - 1]);
		cDVector myDNu = myDeltaDiv2 * myNupm1 * theGradData.mGradHt[j] + myLogh/2 * myNu * myDDelta;
		cDMatrix myD2Nu = myNupm1/2 * (myDDelta * Transpose(theGradData.mGradHt[j]) + theGradData.mGradHt[j] * Transpose(myDDelta));
			myD2Nu += mvDelta/4.0 * myLogh * myNupm1 *  (myDDelta * Transpose(theGradData.mGradHt[j]) + theGradData.mGradHt[j] * Transpose(myDDelta));
			myD2Nu += myDeltaDiv2 * (myDeltaDiv2 - 1.0) * myNupm2 * theGradData.mGradHt[j] * Transpose(theGradData.mGradHt[j]);
			myD2Nu += pow(myLogh/2.0, 2.0) * myNu * myDDelta * Transpose(myDDelta);
			myD2Nu += myDeltaDiv2 * myNupm1 * theHessData.mHessHt[j];

			myHessNu += myDGarch[j] * Transpose(myDNu) + myDNu * Transpose(myDGarch[j]);
			myHessNu += mvGarch[j] * myHessNu;
		}
	// On a calculé Hess h(t)^(Delta/2), il faut calculer Hess H(t)
		
	double	myNu = pow(theData.mHt[theDate], myDeltaDiv2);
	double	myNupm1 = pow(theData.mHt[theDate], myDeltaDiv2 - 1.0);
	double	myNupm2 = pow(theData.mHt[theDate], myDeltaDiv2 - 2.0);
	double myLogh = log(theData.mHt[theDate]);
	cDVector myDNu = myDeltaDiv2 * myNupm1 * theGradData.mCurrentGradVar + myLogh/2 * myNu * myDDelta;
		theHessData.mCurrentHessVar = myHessNu - myNupm1/2 * (myDDelta * Transpose(theGradData.mCurrentGradVar) + theGradData.mCurrentGradVar * Transpose(myDDelta));
		theHessData.mCurrentHessVar -= mvDelta/4.0 * myLogh * myNupm1 *  (myDDelta * Transpose(theGradData.mCurrentGradVar) + theGradData.mCurrentGradVar * Transpose(myDDelta));
		theHessData.mCurrentHessVar -= myDeltaDiv2 * (myDeltaDiv2 - 1.0) * myNupm2 * theGradData.mCurrentGradVar * Transpose(theGradData.mCurrentGradVar);
		theHessData.mCurrentHessVar -= pow(myLogh/2.0, 2.0) * myNu * myDDelta * Transpose(myDDelta);
		theHessData.mCurrentHessVar /= myDeltaDiv2 * myNupm1;
		
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
	}

	void cAparch::ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, cFuncMeanAndVar& theDerivM, cAbstResiduals* theResids)
	{

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
