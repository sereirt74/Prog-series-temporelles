#include "StdAfxRegArchLib.h"

namespace RegArchLib {

	cPolynome::cPolynome(int theDegree)
	{
		if (theDegree < 0)
		{
			mDegree = -1;
			mCoeff = NULL;
		}
		else
		{
			mCoeff = new double[theDegree + 1];
			mDegree = theDegree;
			for (uint i = 0; i <= (uint)mDegree; i++)
				mCoeff[i] = 0.0;
		}
	}

	cPolynome::cPolynome(int theDegree, double* theCoeff)
	{
		if (theDegree < 0)
		{
			mDegree = -1;
			mCoeff = NULL;
		}
		else
		{
			mCoeff = new double[theDegree + 1];
			mDegree = theDegree;
			for (int i = 0; i <= mDegree; i++)
				mCoeff[i] = theCoeff[i];
		}
	}
		
	cPolynome::cPolynome(const cPolynome &theSrc)
	{
		Delete();
		if (theSrc.mDegree >= 0)
		{
			mDegree = theSrc.mDegree;
			mCoeff = new double[mDegree + 1];
			for (int i = 0; i <= mDegree; i++)
				mCoeff[i] = theSrc.mCoeff[i];
		}
	}

	cPolynome::~cPolynome()
	{
		if (mDegree >= 0)
			delete[] mCoeff;
		mCoeff = NULL;
		mDegree = -1;
		MESS_DESTR("cArfima")
	}

	void cPolynome::Delete(void)
	{
		if (mDegree >= 0)
			delete[] mCoeff;
		mDegree = -1;
		mCoeff = NULL;
	}

	void cPolynome::Resize(int theDegree)
	{
		this->Delete();
		if (theDegree >= 0)
		{
			this->mCoeff = new double[theDegree + 1];
			this->mDegree = theDegree;
			for (uint i = 0; i <= (uint)mDegree; i++)
				mCoeff[i] = 0.0;
		}
	}

	double& cPolynome::operator[](int theIndex) const
	{
		return mCoeff[theIndex];
	}

	cPolynome& cPolynome::operator=(const cPolynome& theSrc)
	{
		int myDegree = theSrc.mDegree;
		if (myDegree < 0)
			Delete();
		else
		{
			if (myDegree != mDegree)
				Resize(myDegree);
			for (int i = 0; i <= myDegree; i++)
					mCoeff[i] = theSrc.mCoeff[i];
		}

		return *this;
	}

	cPolynome& cPolynome::operator +=(const cPolynome& theP)
	{

		int myp1 = mDegree;
		int myp2 = theP.mDegree;
		int myp = MAX(myp1, myp2);

		cPolynome myTmpPoly(0);
		myTmpPoly = *this;
		Resize(myp);
		for (int i = 0; i <= myp; i++)
		{
			if (i <= myp1)
				mCoeff[i] = myTmpPoly.mCoeff[i];
			else
				mCoeff[i] = 0.0;
			if (i <= myp2)
				mCoeff[i] += theP.mCoeff[i];

		}
		return *this;
	}

	cPolynome& cPolynome::operator -=(const cPolynome& theP)
	{
		int myp1 = mDegree;
		int myp2 = theP.mDegree;
		int myp = MAX(myp1, myp2);

		cPolynome myTmpPoly(0);
		myTmpPoly = *this;
		Resize(myp);
		for (int i = 0; i <= myp; i++)
		{
			if (i <= myp1)
				mCoeff[i] = myTmpPoly.mCoeff[i];
			else
				mCoeff[i] = 0.0;
			if (i <= myp2)
				mCoeff[i] -= theP.mCoeff[i];

		}
		return *this;
	}

	cPolynome& cPolynome::operator *=(const cPolynome& theP)
	{
		if (mDegree < 0)
			return *this;
		int myp1 = mDegree;
		int myp2 = theP.mDegree;
		int myp = myp1 + myp2;
		cPolynome myTmpPoly(0);
		myTmpPoly = *this;
		Resize(myp);
		for (int i = 0; i <= myp1; i++)
			for (uint j = 0; j <= (uint)myp2; j++)
					mCoeff[i + j] += myTmpPoly.mCoeff[i] * theP.mCoeff[j];

		return *this;
	}

	cPolynome& cPolynome::operator *=(double theVal)
	{
		if (mDegree < 0)
			return *this;
		for (int i = 0; i <= mDegree; i++)
			mCoeff[i] *= theVal;

		return *this;
	}

	cPolynome& cPolynome::operator+(const cPolynome& theP)
	{
		cPolynome* myTmpPoly = new cPolynome(0);
		*myTmpPoly = *this;
		*myTmpPoly += theP;
		return *myTmpPoly;
	}

	cPolynome& cPolynome::operator-(const cPolynome& theP)
	{

		cPolynome* myTmpPoly = new cPolynome(0);
		*myTmpPoly = *this;
		*myTmpPoly -= theP;
		return *myTmpPoly;
	}

	cPolynome& cPolynome::operator*(const cPolynome& theP)
	{
		cPolynome* myTmpPoly = new cPolynome(0);
		*myTmpPoly = *this;
		*myTmpPoly *= theP;
		return *myTmpPoly;
	}

#ifndef _RDLL_
	void cPolynome::Print(void)
	{
	int myDegree = mDegree;
		if (mCoeff[myDegree] != 0)
			cout << mCoeff[myDegree] << "*x^" << myDegree;
		for (int i = myDegree - 1; i >= 1; i--)
		{
			if (mCoeff[i] > 0)
				cout << "+" << mCoeff[i] << "*x^" << i;
			else
				if (mCoeff[i] < 0)
					cout << mCoeff[i] << "*x^" << i;
		}
		if (mCoeff[0] > 0)
			cout << "+" << mCoeff[0] << endl;
		else
			if (mCoeff[0] < 0)
				cout << mCoeff[0] << endl;
			else 
				cout << endl;
	}
#else
	void cPolynome::Print(void)
	{
		uint myDegree = mDegree;
		if ((mCoeff)[myDegree] != 0)
			Rprintf("%f*x^%d", (mCoeff)[myDegree],myDegree);
		for (uint i = myDegree - 1; i >= 1; i--)
		{
			if ((mCoeff)[i] > 0)
				Rprintf("+%f*x^%d", (mCoeff)[i], i);
			else
				if ((mCoeff)[i] < 0)
					Rprintf("%f*x^%d", (mCoeff)[i], i);
		}
		if ((mCoeff)[0] > 0)
			Rprintf("+%f\n", (mCoeff)[0]);
		else
			if ((mCoeff)[0] < 0)
				Rprintf("%f\n", (mCoeff)[0]);
			else 
				Rprintf("\n");
	}
#endif //_RDLL_

	cPolynome& operator *(cPolynome& theP, double theVal)
	{
		cPolynome* myTmpPoly = new cPolynome(0);
		*myTmpPoly = theP;
		*myTmpPoly *= theVal;
		return *myTmpPoly;
	}

	cPolynome& operator *(double theVal, cPolynome& theP)
	{
		cPolynome* myTmpPoly = new cPolynome(0);
		*myTmpPoly = theP;
		*myTmpPoly *= theVal;
		return *myTmpPoly;
	}

	void ComputeDeltaPowD(double theD, uint theDegree, cPolynome& theRes)
	{
		theRes.Resize(theDegree);
		theRes[0] = 1.0;
		for (uint i = 1; i <= theDegree; i++)
			theRes[i] = -theRes[i - 1] * (theD - (double)(i - 1)) / (double)i;
	}

	void ComputeLogDelta(uint theDegree, cPolynome& theRes)
	{
		theRes.Resize(theDegree);
		theRes[0] = 0.0;
		for (uint i = 1; i <= theDegree; i++)
			theRes[i] = -1.0 / (double)i;
	}

	void IncrPowDiv(cPolynome& theNum, cPolynome& theDen, uint theOrder, cPolynome& theQuot, cPolynome& theRest)
	{
	uint myp1 = theNum.mDegree;
	uint myp2 = theDen.mDegree;
		theQuot.Resize(theOrder);
		theRest.Resize(MAX(theOrder + myp2, myp1));
		for (uint k = 0; k <= myp1; k++)
			theRest[k] = theNum[k];
		for (uint k = 0; k <= theOrder; k++)
		{
		double myCoeff = theQuot[k] = theRest[k] / theDen[0];
			if (myCoeff != 0)
			{
				for (uint j = 0; j <= myp2; j++)
					theRest[j+k] -= myCoeff * theDen[j];

			}
		}
	}

	cPolynome& TrunkMult(cPolynome& theP, cPolynome& theQ, uint theMaxDegree)
	{
		cPolynome* myPoly = new cPolynome(theMaxDegree);
		for (uint i = 0; i <= (uint)theP.mDegree; i++)
			for (uint j = 0; j <= (uint)theQ.mDegree; j++)
				if (i + j <= theMaxDegree)
					myPoly->mCoeff[i + j] += (theP.mCoeff)[i] * (theQ.mCoeff)[j];
		return *myPoly;
	}

	double cPolynome::BackwardPolOp(const cDVector& theYt, uint theIndex0, double thePow) const
	{
		int myDegree = mDegree;
		double myRes = mCoeff[0];
		if (thePow == 1.0)
		{
			for (int i = 1; i <= MIN(myDegree, (int)theIndex0); i++)
				myRes += mCoeff[i] * theYt[(int)theIndex0 - i];
		}
		else
		{
			for (int i = 1; i <= MIN(myDegree, (int)theIndex0); i++)
				myRes += mCoeff[i] * pow(theYt[(int)theIndex0 - i], thePow);
		}

		return myRes;
	}

	cPolynome& TrunkPoly(cPolynome& theP, uint theMaxDegree)
	{
		static cPolynome myPoly((uint)(MIN(theP.mDegree, theMaxDegree)));
		for (uint i = 0; i <= (uint)(MIN(theP.mDegree, theMaxDegree)); i++)
			myPoly[i] = theP[i];
		return myPoly;
	}


} //namespace
