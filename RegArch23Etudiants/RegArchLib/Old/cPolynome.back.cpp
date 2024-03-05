#include "StdAfxRegArchLib.h"

namespace RegArchLib {
	using namespace std;

	cPolynome::cPolynome(int theDegree)
	{
		if (theDegree < 0)
		{
			mCoeff = vector<double>();
		}
		else
		{
			mCoeff = vector<double>(theDegree + 1);
		}
	}

	cPolynome::cPolynome(int theDegree, double* theCoeff)
	{
		if (theDegree < 0)
		{
			mCoeff = vector<double>();
		}
		else
		{
			mCoeff = vector<double>(theDegree + 1);
			for (uint i = 0; i <= theDegree; i++)
				mCoeff[i] = theCoeff[i];
		}
	}

	cPolynome::cPolynome(cDVector& theCoeff)
	{
		int mySize = theCoeff.GetSize();
		if (mySize > 0)
		{
			mCoeff = vector<double>(mySize);
			for (uint i = 0; i < mySize; i++)
				mCoeff[i] = theCoeff[i];
		}
		else
		{
			mCoeff = vector<double>();
		}
	}
		
	cPolynome::cPolynome(const cPolynome &theSrc)
	{
		*this = theSrc;
	}

	cPolynome::~cPolynome()
	{
		mCoeff.~vector();
	}

	void cPolynome::Delete(void)
	{
		mCoeff.~vector();
	}

	void cPolynome::Resize(uint theDegree)
	{
		if (theDegree >= 0)
		{
			mCoeff.resize(theDegree + 1);
		}
	}

	int cPolynome::GetNDegree(void) const
	{
		return mCoeff.size() - 1;
	}

	double& cPolynome::operator[](int theIndex) 
	{
		return mCoeff[theIndex];
	}

	cPolynome& cPolynome::operator=(const cPolynome& theSrc)
	{
		mCoeff = theSrc.mCoeff;
		return *this;
	}

	cPolynome& cPolynome::operator+(const cPolynome& theP)
	{
	int myp1 = GetNDegree();
	int myp2 = theP.GetNDegree();
	int myp = MAX(myp1, myp2);
	static cPolynome myTmpPoly;
		myTmpPoly.Resize(myp);
	
		for (int i = 0; i <= myp; i++)
		{
			if (i <= myp1)
				myTmpPoly.mCoeff[i] =  mCoeff[i];
			if (i <= myp2)
				myTmpPoly.mCoeff[i] += theP.mCoeff[i];
		}
		return myTmpPoly;
	}

	cPolynome& cPolynome::operator-(const cPolynome& theP)
	{
	int myp1 = GetNDegree();
	int myp2 = theP.GetNDegree();
	int myp = MAX(myp1, myp2);

	static cPolynome myTmpPoly;
		myTmpPoly.Resize(myp);
		for (int i = 0; i <= myp; i++)
		{
			if (i <= myp1)
				myTmpPoly.mCoeff[i] = mCoeff[i];
			if (i <= myp2)
				myTmpPoly.mCoeff[i] -= theP.mCoeff[i];
		}
		return myTmpPoly;
	}

	cPolynome& cPolynome::operator*(const cPolynome& theP)
	{
	int myp1 = GetNDegree();
	int myp2 = theP.GetNDegree();
	int myp = myp1 + myp2;
	static cPolynome myTmpPoly;
		myTmpPoly.Resize(myp);
		for (int i = 0; i <= myp1; i++)
			for (int j = 0; j <= myp2; j++)
				myTmpPoly.mCoeff[i + j] += mCoeff[i] * theP.mCoeff[j];
		return myTmpPoly;
	}

#ifndef _RDLL_
	void cPolynome::Print(void)
	{
	uint myDegree = GetNDegree();
		if (mCoeff[myDegree] != 0)
			cout << mCoeff[myDegree] << "*x^" << myDegree;
		for (uint i = myDegree - 1; i >= 1; i--)
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
		uint myDegree = GetNDegree();
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
	int myp = theP.GetNDegree();
	static cPolynome myTmpPoly(myp);

		for (int i = 0; i <= myp; i++)
			myTmpPoly.mCoeff[i] = theP.mCoeff[i] * theVal;
		return myTmpPoly;
	}

	cPolynome& operator *(double theVal, cPolynome& theP)
	{
	int myp = theP.GetNDegree();
		static cPolynome myTmpPoly(myp);

		for (int i = 0; i <= myp ; i++)
			myTmpPoly.mCoeff[i] = theP.mCoeff[i] * theVal;
		return myTmpPoly;
	}

	cPolynome& cPolynome::operator +=(const cPolynome& theP)
	{
	int myp1 = GetNDegree();
	int myp2 = theP.GetNDegree();
	int myp = MAX(myp1, myp2);
		
		cPolynome myTmpPoly(myp);
		for (int i = 0; i <= myp; i++)
		{
			if (i <= myp1)
				myTmpPoly[i] = mCoeff[i];
			if (i <= myp2)
				myTmpPoly[i] += theP.mCoeff[i];
		}

		*this = myTmpPoly;
		myTmpPoly.Delete();
		return *this;
	}

	cPolynome& cPolynome::operator -=(const cPolynome& theP)
	{
	int myp1 = GetNDegree();
	int myp2 = theP.GetNDegree();
	int myp = MAX(myp1, myp2);
		
		cPolynome myTmpPoly(myp);
		for (int i = 0; i <= myp; i++)
		{
			if (i <= myp1)
				myTmpPoly.mCoeff[i] = mCoeff[i];
			if (i <= myp2)
				myTmpPoly.mCoeff[i] -= theP.mCoeff[i];
		}

		*this = myTmpPoly;
		myTmpPoly.Delete();
		return *this;
	}

	cPolynome& cPolynome::operator *=(const cPolynome& theP)
	{
	cPolynome myAuxPol(0);
		myAuxPol = (*this) * theP;
		*this = myAuxPol;
		myAuxPol.Delete();
		return *this;
	}

	cPolynome& cPolynome::operator *=(double theVal)
	{
		for (int i = 0 ; i < mCoeff.size() ; i++)
			mCoeff[i] *= theVal;
		return *this;
	}

	void ComputeDeltaPowD(double theD, uint theDegree, cPolynome& theRes)
	{
		theRes.Resize(theDegree);
		theRes[0] = 1.0;
		for (int i = 1; i <= theDegree; i++)
			theRes[i] = -theRes[i - 1] * (theD - (double)(i - 1)) / (double)i;
	}

	void ComputeLogDelta(uint theDegree, cPolynome& theRes)
	{
		theRes.Resize(theDegree);
		theRes[0] = 0.0;
		for (int i = 1; i <= theDegree; i++)
			theRes[i] = -1.0 / (double)i;
	}

	/*
	void IncrPowDiv(cPolynome& theNum, cPolynome& theDen, uint theOrder, cPolynome& theQuot, cPolynome& theRest)
	{
		//Division of two polynomials by increasing powers

		uint p1 = theNum.GetNDegree();

		theQuot.Resize(theOrder + p1);
		theRest.Resize(theOrder + p1);
		cPolynome myAux(theOrder + p1);
		cPolynome myMonome(theOrder + p1);

		for (uint i = 0; i <= p1; i++)
			theRest[i] = theNum[i];
		uint j = 0;
		while (j <= theOrder)
		{
			uint k = j;
		bool myFini = (k > theOrder + p1);
			while (! myFini)
			{
				if (theRest[k] == 0)
				{
					k++;
					myFini = (k > theOrder + p1);
				}
				else
					myFini = true;
			}
			if (k <= theOrder + p1)
			{
				(myMonome.mCoeff)[k] = (theQuot.mCoeff)[k] = theRest[k] / theDen[0];
				myAux = myMonome * theDen;
				(myMonome.mCoeff)[k] = 0;
				theRest -= myAux;
				j++;
			}
			else
				j = theOrder + 1;

		}
		myAux.Delete();
		myMonome.Delete();
	}
*/
	void IncrPowDiv(cPolynome& theNum, cPolynome& theDen, uint theOrder, cPolynome& theQuot, cPolynome& theRest)
	{
	int myp1 = theNum.GetNDegree();
	int myp2 = theDen.GetNDegree();
		theQuot.Resize(theOrder);
		theRest.Resize(MAX(theOrder + myp2, myp1));
		for (int k = 0; k <= myp1; k++)
			theRest[k] = theNum[k];
		for (int k = 0; k <= theOrder; k++)
		{
		double myCoeff = theQuot[k] = theRest[k] / theDen[0];
			if (myCoeff != 0)
			{
				for (int j = 0; j <= myp2; j++)
					theRest[j+k] -= myCoeff * theDen[j];

			}
		}
	}


	cPolynome& TrunkMult(cPolynome& theP, cPolynome& theQ, uint theMaxDegree)
	{
	static cPolynome myTmpPoly(theMaxDegree);

		for (int i = 0; i <= theP.GetNDegree(); i++)
			for (int j = 0; j <= theQ.GetNDegree(); j++)
				if (i + j <= theMaxDegree)
					myTmpPoly.mCoeff[i + j] += (theP.mCoeff)[i] * (theQ.mCoeff)[j];
		return myTmpPoly;
	}

	double cPolynome::BackwardPolOp(const cDVector& theYt, uint theDate, double thePow) const
	{
	int myDegree = GetNDegree();
	double myRes = mCoeff[0];
		if (thePow == 1.0)
		{
			for (int i = 1; i <= MIN(myDegree, theDate); i++)
				myRes += mCoeff[i] * theYt[theDate - i];
		}
		else
		{
			for (int i = 1; i <= MIN(myDegree, theDate); i++)
				myRes += mCoeff[i] * pow(theYt[theDate - i], thePow);
		}

		return myRes;
	}

	cPolynome& TrunkPoly(cPolynome& theP, uint theMaxDegree)
	{
	static cPolynome myTmpPoly(MIN(theP.GetNDegree(), theMaxDegree));

		for (int i = 0; i <= MIN(theP.GetNDegree(), theMaxDegree); i++)
			myTmpPoly[i] = theP[i];
		return myTmpPoly;
	}

	cPolynome& InsertCoef(vector<double> theCoef)
	{
		static cPolynome myTmpPoly(theCoef.size() - 1);
		myTmpPoly.mCoeff = theCoef;
		return myTmpPoly;
	}

} //namespace
