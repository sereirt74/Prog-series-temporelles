#include "StdAfxRegArchLib.h"

namespace RegArchLib {

	double TrueLogDensity(double theX, const cDVector theTeta, LogDensFunction theLogDens, StdevFunction theStdev)
	{
	double myStdev = theStdev(theTeta);

		return theLogDens(theX*myStdev, theTeta) + log(myStdev);
	}

	double DiffTrueDensity(double theX, const cDVector& theTeta, LogDensFunction theDiffLogDens, StdevFunction theStdev)
	{
	double myStdev = theStdev(theTeta);
		return myStdev * theDiffLogDens(theX*myStdev, theTeta);
	}

	double Diff2TrueDensity(double theX, const cDVector& theTeta, LogDensFunction theDiff2LogDens, StdevFunction theStdev)
	{
	double myStdev = theStdev(theTeta);
		return myStdev * myStdev * theDiff2LogDens(theX*myStdev, theTeta);
	}

	void GradTrueLogDensity(double theX, const cDVector& theTeta, LogDensFunction theDiffLogDens, StdevFunction theStdev, GradStdevFunction theGradStdev, GradLogDensFunction theGradLogDensity, cDVector& theGrad)
	{
	double myStdev = theStdev(theTeta);
	double myX = myStdev * theX;
	cDVector myGradStd = cDVector(theGrad.GetSize());
		theGradStdev(theTeta, myGradStd);
	cDVector myGradLogDens = cDVector(theGrad.GetSize());
		theGradLogDensity(myX, theTeta, myGradLogDens);
		theGrad = (theX*theDiffLogDens(myX, theTeta) + 1 / myStdev)*myGradStd + myGradLogDens;
	}

	void GradTrueDiffLogDensity(double theX, const cDVector& theTeta, LogDensFunction theDiffLogDens, LogDensFunction theDiff2LogDens, StdevFunction theStdev, GradStdevFunction theGradStdev, GradLogDensFunction theGradDiffLogDensity, cDVector& theGrad)
	{
	double myStdev = theStdev(theTeta);
	double myX = myStdev * theX;
	cDVector myGradStd = cDVector(theGrad.GetSize());
		theGradStdev(theTeta, myGradStd);
	cDVector myGradDifLogDens = cDVector(theGrad.GetSize());
		theGradDiffLogDensity(myX, theTeta, myGradDifLogDens);
	double myAux = theDiffLogDens(myX, theTeta) + myX *theDiff2LogDens(myX, theTeta);
		theGrad = myAux * myGradStd + myStdev*myGradDifLogDens;
	}

	void HessTrueLogDensity(double theX, const cDVector& theTeta, LogDensFunction theDiffLogDens, LogDensFunction theDiff2LogDens, StdevFunction theStdev, GradStdevFunction theGradStdev, GradLogDensFunction theGradDiffLogDensity, HessStdevFunction theHessStdev, HessLogDensFunction theHessLogDens, cDMatrix& theHess)
	{
	int myNParam = theHess.GetNRow();
	double myStdev = theStdev(theTeta);
	double myX = myStdev * theX;
	cDVector myGradStd = cDVector(myNParam);
		theGradStdev(theTeta, myGradStd);
	cDVector myGradDifLogDens = cDVector(myNParam);
		theGradDiffLogDensity(myX, theTeta, myGradDifLogDens);
	double myAux = theX * theX * theDiff2LogDens(myX, theTeta) - 1.0 / (myStdev*myStdev);
		theHess = myAux * myGradStd * Transpose(myGradStd);
		theHess += theX * (myGradDifLogDens * Transpose(myGradStd) + myGradStd * Transpose(myGradDifLogDens));
	cDMatrix myHessStd = cDMatrix(myNParam, myNParam);
		theHessStdev(theTeta, myHessStd);
		myAux = 1 / myStdev + theX * theDiffLogDens(myX, theTeta);
		theHess += myAux * myHessStd;
	cDMatrix myHessLogDens = cDMatrix(myNParam, myNParam);
		theHessLogDens(myX, theTeta, myHessLogDens);
		theHess += myHessLogDens;
	}


} // namespace

