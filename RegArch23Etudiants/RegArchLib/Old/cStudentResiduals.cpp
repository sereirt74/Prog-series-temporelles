#include "StdAfxRegArchLib.h"
/*!
 \file cStudentResiduals.cpp
 \brief implementation of the class for Student conditional distribution.

 \author Jean-Baptiste DURAND, Ollivier TARAMASCO
 \date dec-18-2006 - Last change apr-7-2019
*/
namespace RegArchLib {
	/*!
	 * \fn cStudentResiduals::cStudentResiduals(double theDof, bool theSimulFlag)
	 * \param double theDof: number of degrees of freedom 
	 * \param bool theSimulFlag: true if created for simulation
	 * \details: mvBool is initialised by ce cAbstResiduals constructor
	 */
	cStudentResiduals::cStudentResiduals(double theDof, bool theSimulFlag): cAbstResiduals(eStudent, NULL, theSimulFlag)
	{
		mDistrParameter = new cDVector(1, theDof);
		MESS_CREAT("cStudentResiduals")
	}

	/*!
	 * \fn cStudentResiduals::cStudentResiduals(const cDVector* theDistrParameter, bool theSimulFlag): cAbstResiduals(eStudent, theDistrParameter, theSimulFlag)
	 * \param const cDVector* theDistrParameter: theDistrParameter[0] = d.o.f.
	 * \param bool theSimulFlag: true if created for simulation
	 * \details: mvBool is initialised by ce cAbstResiduals constructor
	 */
	cStudentResiduals::cStudentResiduals(cDVector* theDistrParameter, bool theSimulFlag): cAbstResiduals(eStudent, theDistrParameter, theSimulFlag)
	{
		MESS_CREAT("cStudentResiduals")
	}

	cStudentResiduals::cStudentResiduals(const cStudentResiduals& theSrc) : cAbstResiduals(eStudent, theSrc.mDistrParameter, (theSrc.mtR != NULL))
	{
		MESS_CREAT("cStudentResiduals")
	}


	/*!
	 * \fn cStudentResiduals::~cStudentResiduals()
	 */
	cStudentResiduals::~cStudentResiduals()
	{
		MESS_DESTR("cStudentResiduals")
	}

	/*!
	 * \fn cAbstCondVar* cStudentResiduals::::PtrCopy()
	 */

/*
cAbstResiduals* cStudentResiduals::PtrCopy() const
	{
		return cAbstResidualsPtrCopy<cStudentResiduals>();
	}
*/
	/*!
	 * \fn void cStudentResiduals::Generate(const uint theNSample, cDVector& theYt) const
	 * \param const uint theNSample: the sample size
	 * \param cDVector& theYt: the output vector
	 */
	void cStudentResiduals::Generate(const uint theNSample, cDVector& theYt) const
	{
		theYt.ReAlloc(theNSample) ;
		if ((*mDistrParameter)[0] <= 2.0)
			throw cError("wrong d.o.f.") ;

	double myStd = sqrt((*mDistrParameter)[0]/((*mDistrParameter)[0]-2.0)) ;
		for (uint t = 0 ; t < theNSample ; t++)
			theYt[t] = gsl_ran_tdist(mtR, (*mDistrParameter)[0])/myStd ;
	}

	/*!
	 * \fn void cStudentResiduals::Print(ostream& theOut) const
	 * \param ostream& theOut: the output stream, default cout.
	 */
#ifndef _RDLL_
	void cStudentResiduals::Print(ostream& theOut) const
	{
		theOut << "Conditional Student Distribution with " << (*mDistrParameter)[0] << " d. o. f." << endl ;
	}
#else
	void cStudentResiduals::Print(void)
	{
		Rprintf("Conditional Student Distribution with %f d.o.f.\n", (*mDistrParameter)[0]);
	}
#endif // _RDLL_

	void cStudentResiduals::SetDefaultInitPoint(void) 
	{
		(*mDistrParameter)[0] = 10.0 ;
	}

	double cStudentResiduals::LogDensity(double theX) const
	{
	double myStd = sqrt((*mDistrParameter)[0]/((*mDistrParameter)[0]-2.0)) ;
		return StudentLogDensity(theX*myStd, (*mDistrParameter)[0]) + log(myStd) ;

	}

	/*!
	 * \fn double cStudentResiduals::GetNParam(void) const
	 * \param void.
	 * \brief return 1: One parameter for St(n) residuals.
	 */
	uint cStudentResiduals::GetNParam(void) const
	{
		return 1 ;
	}

	/*!
	 * \fn static void StudentGradLogDensity(double theX, double theDof, cDVector& theGrad)
	 * \brief Compute the derivative of log density of a Student distribution (unstandardized)
	 * with respect to the random variable (theGrad[0]) \e and the gradient
	 * of log density with respect to the model parameters (other components in theGrad)
	 * \param theX double: value of the random variable
	 * \param theDof double: value of the distribution parameter
	 * \param theGrad cDVector&: concatenation of derivatives with respect to the random variable and the model parameters
	 */
	double cStudentResiduals::DiffLogDensity(double theX) const
	{
		double myDof = (*mDistrParameter)[0];
		return -(myDof + 1)*theX / (theX*theX + myDof - 2.0);
	}

	/*!
	 * \fn static void GradLogDensity(double theX, cDVector& theGrad, cDVector& theDistrParam)
	 * \brief Compute the derivative of log density with respect to the random variable (theGrad[0]) \e and the gradient
	 * of log density with respect to the model parameters (other components in theGrad)
	 * \param theX double: value of the random variable
	 * \param theGrad cDVector&: concatenation of derivatives with respect to the random variable and the model parameters
	 * \param theDistrParam cDVector&: value of the distribution parameters
	 */
	static void GradLogDensity(double theX, cDVector& theGrad, const cDVector& theDistrParam, uint theBegIndex)
	{
	double	myDof = theDistrParam[0];
	double myt2 = theX*theX;
	double myAux1 = myt2 + myDof - 2;
	double myRes = -log(myAux1 / (myDof - 2)) / 2;
		myRes += (gsl_sf_psi((myDof + 1) / 2) - gsl_sf_psi(myDof / 2)) / 2;
		myRes += (myDof*myt2 - myDof + 2) / ((myDof - 2)*(myt2 + myDof - 2)) / 2;
		theGrad[theBegIndex] = myRes;
	}

	/*!
	 * \fn void cStudentResiduals::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData)
	 * \brief Compute the derivative of log density with respect to the random variable (theGradData[0]) \e and the gradient
	 * of log density with respect to the model parameters (other components in theGradData)
	 * \param theDate uint: time at which gradient is computed
	 * \param theValue const cRegArchValue&: value of the random variable
	 * \param theGradData cRegArchGradient&: concatenation of derivatives with respect to the random variable and the model parameters
	 */
	void cStudentResiduals::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData) const
	{
 	uint myBegIndex = theGradData.GetNMeanParam() + theGradData.GetNVarParam();
		GradLogDensity(theValue.mEpst[theDate], theGradData.mCurrentGradLogDens, *mDistrParameter, myBegIndex) ;
		theGradData.mCurrentDiffLogDensity = DiffLogDensity(theValue.mEpst[theDate]);
	}

	void cStudentResiduals::RegArchParamToVector(cDVector& theDestVect, uint theIndex) const
	{
		if (theDestVect.GetSize() < theIndex )
			throw cError("Wrong size") ;
		theDestVect[theIndex] = (*mDistrParameter)[0] ;
	}

	void cStudentResiduals::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
		if (theIndex >= theSrcVect.GetSize())
			throw cError("Wrong size") ;
		(*mDistrParameter)[0] = theSrcVect[theIndex] ;
	}

	/*
	 * (2*sqrt(n-2)*gamma((n+1)/2))/(sqrt(PI)*(n-1)*gamma(n/2))
	 */
	double cStudentResiduals::ComputeEspAbsEps(void)
	{
	double myDof = (*mDistrParameter)[0] ;
		return gsl_sf_gamma((myDof+1.0)/2.0)/gsl_sf_gamma(myDof/2.0)
			*2.0*sqrt(myDof-2.0)/(myDof-1.0)/SQRT_PI ;
	}

	/*
	 *  [(n^2 - 3*n + 2)*(psi0((n+1)/2) - psi0(n/2)) - n + 3]/[sqrt(PI)*sqrt(n-2)*(n-1)^2]*gamma((n+1)/2)/gamma(n/2)
	 */
	void cStudentResiduals::ComputeGradBetaEspAbsEps(cDVector& theGrad)
	{
	double myDof = (*mDistrParameter)[0],
			myAux = myDof*myDof-3*myDof+2 ;
	
		theGrad[0] = (myAux*gsl_sf_psi((myDof+1)/2.0)-myAux*gsl_sf_psi(myDof/2.0)- myDof + 3.0)
			*gsl_sf_gamma((myDof+1)/2.0)
			/(SQRT_PI*sqrt(myDof-2)*(myDof*myDof-2*myDof+1)*gsl_sf_gamma(myDof/2.0)) ;
	}

	void cStudentResiduals::ComputeHessBetaEspAbsEps(cDMatrix &theHess)
	{
		double myDof = (*mDistrParameter)[0];
		double myAux1 = pow(myDof, 4) - 6 * pow(myDof, 3) + 13 * pow(myDof, 2) - 12 * myDof + 4;
		double myPsiks2 = gsl_sf_psi(myDof/2.0);
		double myPsi1ks2 = gsl_sf_psi_1(myDof / 2);
		double myGammaks2 = gsl_sf_gamma(myDof / 2);
		double myGammakp1s2 = gsl_sf_gamma((myDof+1)/2);

		theHess[0][0]=(((pow(myDof, 4) - 6 * pow(myDof, 3) + 13 * pow(myDof, 2) - 12 * myDof + 4)*gsl_sf_psi_1((myDof + 1) / 2)
			+ (pow(myDof, 4) - 6 * pow(myDof, 3) + 13 * pow(myDof, 2) - 12 * myDof + 4)*pow(gsl_sf_psi((myDof + 1) / 2),2)
			+ (-2 * myPsiks2*pow(myDof, 4) + (12 * myPsiks2 - 2)*pow(myDof, 3)
				+ (12 - 26 * myPsiks2)*pow(myDof, 2) + (24 * myPsiks2 - 22)*myDof - 8 * myPsiks2
				+ 12)*gsl_sf_psi((myDof + 1) / 2) + (pow(myPsiks2, 2) - myPsi1ks2)*pow(myDof, 4)
			+ (6 * myPsi1ks2 - 6 * pow(myPsiks2, 2) + 2 * myPsiks2)*pow(myDof, 3)
			+ (-13 * myPsi1ks2 + 13 * pow(myPsiks2, 2) - 12 * myPsiks2 + 3)*pow(myDof, 2)
			+ (12 * myPsi1ks2 - 12 * pow(myPsiks2, 2) + 22 * myPsiks2 - 18)*myDof - 4 * myPsi1ks2 + 4 * pow(myPsiks2, 2)
			- 12 * myPsiks2 + 23)*myGammakp1s2 / (sqrt(PI)*sqrt(myDof - 2)*(2 * myGammaks2*pow(myDof, 4)
				- 10 * myGammaks2*pow(myDof, 3) + 18 * myGammaks2*pow(myDof, 2) - 14 * myGammaks2*myDof + 4 * myGammaks2)));
	}

	double cStudentResiduals::Diff2LogDensity(double theX) const
	{
		double myDof = (*mDistrParameter)[0];
		double myt2 = theX*theX;
		return (myDof + 1)*(myt2 - myDof + 2) / pow(myt2 + myDof - 2, 2);
	}

	void cStudentResiduals::GradDiffLogDensity(double theX, const cDVector& theDistrParam, cDVector& theGrad)
	{
	double myDof = theDistrParam[0];
	double myt2 = theX*theX;
		theGrad[0] = -(theX*(myt2- 3)) / pow(myt2 + myDof - 2, 2);
	}

	static void HessLogDensity(double theX, cDMatrix& theHess, const cDVector& theDistrParam, uint theBegIndex)
	{
	double myDof = theDistrParam[0];
	double myt2 = theX*theX;
	double myAux1 = myt2 + myDof - 2;
	double myRes = ((myDof - 4)*myt2*myt2 - (myDof - 2)*(4*myt2-myDof+2))/(2*pow((myDof-2)*myAux1, 2));
		myRes += (gsl_sf_psi_1((myDof+1)/2)-gsl_sf_psi_1(myDof/2)) / 4;
		theHess[theBegIndex][theBegIndex] = myRes;
	}

	void cStudentResiduals::ComputeHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
	{
	double myX = theData.mEpst[theDate];
	cDVector myGradDiffLogDens = cDVector(1);
		GradDiffLogDensity(myX, *mDistrParameter, myGradDiffLogDens);
	uint myBegIndex = theHessData.GetNMeanParam() + theHessData.GetNVarParam();
		theHessData.mCurrentGradDiffLogDensity[myBegIndex] = myGradDiffLogDens[0];
		HessLogDensity(myX, theHessData.mCurrentHessDens, *mDistrParameter, myBegIndex);
	}

	void cStudentResiduals::GetParamName(uint theIndex, char** theName)
	{
		uint myIndex = theIndex;
		sprintf(theName[myIndex++], "Student d.o.f.");
	}

	void  cStudentResiduals::GetParamName(uint theIndex, string theName[])
	{
	uint myIndex = theIndex;
	char myChar[100];
		sprintf(myChar, "Student d.o.f.");
		theName[myIndex++] = myChar;
	}


}//namespace
