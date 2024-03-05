// Test.cpp : définit le point d'entrée pour l'application console.
//

#include "StdAfxTestCPlusPlus.h"
 
using namespace std;
using namespace ErrorNameSpace;
using namespace VectorAndMatrixNameSpace;
using namespace RegArchLib ;

#ifdef _MSC_VER
int _tmain(int argc, _TCHAR* argv[])
#else
int main(int argc, char* argv[])
#endif // _MSC_VER
{

	cout.precision(4);
	fixed;

	cConst myConst(10);

//	cAr	myAr(2);
	cAr	myAr(1);
	myAr.Set(.8, 0);
//	myAr.Set(-.2, 1);

	cMa myMa(2);
//	cMa myMa(1);
	myMa.Set(0.4, 0);
	myMa.Set(0.6, 1) ;

//	cArfima myArfima(2, 2, 0.3, 10);
//	cArfima myArfima(1, 0, 0.0, 10);
//	cArfima myArfima(0, 1, 0.0, 10);
	cArfima myArfima(1, 1, 0.3, 5);	
	myArfima.Set(-.2, 0, 0);
	myArfima.Set(0.4, 0, 1);
//	myArfima.Set(-.4, 0, 1);
//	myArfima.Set(.5, 1, 1);

	cStdDevInMean myStdDevInMean(.6);
	cVarInMean myVarInMean(.6);
	cLinReg myBeta(2);
	myBeta.Set(3, 0, 0);
	myBeta.Set(-0.5, 0, 1);
	bool myLinRegBool = false;

	cConstCondVar myConstVar(1.0);

	cArch myArch(1);
	myArch.Set(0.4, 0, 0);
	myArch.Set(0.6, 0, 1);

	cGarch myGarch(1, 1);
	myGarch.Set(0.1, 0, 0);
	myGarch.Set(0.1, 0, 1);
	myGarch.Set(0.8, 0, 2);
	
	cNgarch myNgarch(1, 1);
	myNgarch.Set(0.1, 0, 0);
	myNgarch.Set(0.5, 0, 1);
	myNgarch.Set(0.1, 0, 2);
	myNgarch.Set(0.8, 0, 3);

	cNormResiduals myNormResid;

	cAparch myAparch(1, 1);
	myAparch.Set(0.1, 0, 0); // cste ;
	myAparch.Set(1.3, 0, 1); // Delta ;
	myAparch.Set(0.1, 0, 2); // Arch ;
	myAparch.Set(0.3, 0, 3); // Gamma ;
	myAparch.Set(0.8, 0, 4); // Garch ;


	cTarch myTarch(1);
	myTarch.Set(0.1, 0, 0);
	myTarch.Set(0.4, 0, 1);
	myTarch.Set(0.8, 0, 2);

	cUgarch myUgarch(false, 2, 1, 1);
	myUgarch.Set(.1, 0, 1);
	myUgarch.Set(.2, 1, 1);
	myUgarch.Set(.2, 0, 2);
	myUgarch.Set(.7, 0, 3);
	bool myUgarchBool = false;

	cFigarch myFigarch(1, 1, 0.3, 5);
	myFigarch.Set(.1, 0, 0);
	myFigarch.Set(.1, 0, 1);
	//		myFigarch.Set(.2, 1, 1);
	myFigarch.Set(.5, 0, 2);
	//		myFigarch.Set(.1, 1, 2);

	/*
		eLoggarch

	*/
	cGtarch myGtarch(1, 1);
	myGtarch.Set(0.1, 0, 0); // CST
	myGtarch.Set(0.1, 0, 1); // ARCH+
	myGtarch.Set(0.2, 0, 2); // ARCH-
	myGtarch.Set(0.7, 0, 3); // GARCH

	cTsgarch myTsgarch(1, 1);
	myTsgarch.Set(0.1, 0, 0); // CST
	myTsgarch.Set(0.1, 0, 1); // ARCH
	myTsgarch.Set(0.8, 0, 2); // GARCH

	cNagarch myNagarch(1, 1);
	myNagarch.Set(0.1, 0, 1); // CST
	myNagarch.Set(0.2, 0, 2); // ARCH
	myNagarch.Set(0.2, 0, 3); // GAMMA
	myNagarch.Set(0.7, 0, 4); // GARCH

	cSqrgarch mySqrgarch(1, 1);
	mySqrgarch.Set(0.1, 0, 1); // CST
	mySqrgarch.Set(0.2, 0, 2); // ARCH
	mySqrgarch.Set(0.2, 0, 3); // GAMMA
	mySqrgarch.Set(0.7, 0, 4); // GARCH

	cStgarch myStgarch(1, 1);
	mySqrgarch.Set(0.1, 0, 1); // CST
	mySqrgarch.Set(0.2, 0, 2); // ARCH
	mySqrgarch.Set(0.2, 0, 3); // GAMMA
	mySqrgarch.Set(0.7, 0, 4); // GARCH

	cLoggarch myLoggarch(1, 1);
	myLoggarch.Set(-2.3, 0, 0); // CST
	myLoggarch.Set(-2.3, 0, 1); // ARCH
	myLoggarch.Set(-1.5, 0, 2); // GARCH

	cStudentResiduals myStudent(5, true);
	cGedResiduals myGedResiduals(1.5);
	cMixNormResiduals myMixNorm(.5, 1, 2);

	cEgarch myEgarch(&myNormResid, 1, 1);
	//cEgarch myEgarch(&myStudent, 1, 1);	
	//cEgarch myEgarch(&myMixNorm, 1, 1);	
//	cEgarch myEgarch(&myGedResiduals, 1, 1);
	myEgarch.Set(0.0001, 0, 1); // cste
	myEgarch.Set(1.0, 0, 2); // ARCH
	myEgarch.Set(0.75, 0, 3); // GARCH
	myEgarch.Set(-0.3, 0, 4); // TETA
	myEgarch.Set(0.1, 0, 5); //GAMMA

	cCondMean myCondMean;

//		myCondMean.AddOneMean(myConst);
//	    myCondMean.AddOneMean(myAr);
//		myCondMean.AddOneMean(myMa);
//		myCondMean.AddOneMean(myArfima);
//		myCondMean.AddOneMean(myStdDevInMean);
//		myCondMean.AddOneMean(myVarInMean);
//		myCondMean.AddOneMean(myBeta);
//		myLinRegBool = true;

	cRegArchModel myModel;

	myModel.SetMean(myCondMean);


	myModel.SetVar(myConstVar);	
//	myModel.SetVar(myArch);
//	myModel.SetVar(myGarch);
//	myModel.SetVar(myNgarch) ;
//	myModel.SetVar(myEgarch);
//	myModel.SetVar(myAparch);
//	myModel.SetVar(myTarch);
//	myModel.SetVar(myFigarch) ;
//	myModel.SetVar(myUgarch) ;
//	myUgarchBool = true;
//	myModel.SetVar(myGtarch);
//	myModel.SetVar(myTsgarch);
//	myModel.SetVar(myNagarch);
//	myModel.SetVar(mySqrgarch);
//	myModel.SetVar(myStgarch);
//	myModel.SetVar(myLoggarch);

	myModel.SetResid(myNormResid);
//		myModel.SetResid(myStudent) ;
//		myModel.SetResid(myMixNorm);
//		myModel.SetResid(myGedResiduals);

	cout << "Modele : " << endl;
	myModel.Print();

	uint myNSimul = 10;

	cDMatrix myXt = cDMatrix(myNSimul, 2, 1);
	if (myLinRegBool)
	{
		for (uint i = 0; i < myNSimul; i++)
		{
			myXt[i][0] = (double)(i + 1);
			myXt[i][1] = (double)(i + 1)*(i + 1);
		}
	}

	cRegArchValue myValue(myNSimul);
	if (myLinRegBool)
		myValue.mXt = myXt;
	else
		myValue.mXt = NULL;

	if (myUgarchBool)
		myValue.mXvt = myXt;
	else
		myValue.mXvt = NULL;

	RegArchSimul(myNSimul, myModel, myValue);

	//	cout << myValue << endl;

	cRegArchGradient myGradient = cRegArchGradient(&myModel);
	cRegArchGradient myGradNum = cRegArchGradient(&myModel);
	cRegArchHessien myHessien = cRegArchHessien(&myModel);
	cRegArchHessien myHessNum = cRegArchHessien(&myModel);
	uint myNParam = myModel.GetNParam();
	uint myNDistrParam = myModel.mResids->GetNParam();
	uint myNMeanVarParam = myNParam - myNDistrParam;
	cDVector myGradlt(myNParam);
	cDMatrix myHesslt(myNParam, myNParam);
	double myh = 1e-3;
	cDVector myVectParam(myNParam);
	cNumericDerivative myNumDeriv = cNumericDerivative(myNParam, myNDistrParam, myh, myValue.mYt, &myValue.mXt, &myValue.mXvt);
	uint myDate = 7;
	cRegArchValue myAuxValue(myValue);
	myModel.mVar->UpdateProxyVarParameters();
	if (myModel.mMean != NULL)
		myModel.mMean->UpdateProxyMeanParameters();

	for (uint t = 0; t < myNSimul; t++)
	{
		myModel.ComputeGradAndHess(t, myValue, myGradient, myHessien, myModel.mResids);
		myModel.NumericComputeGradAndHess(t, myAuxValue, myGradNum, myHessNum, myModel.mResids, myNumDeriv);
		if (t >= myNSimul - 10) //  (t <= 10)
		{
			cout << "t = " << t << endl;
			cout << "------" << endl;
			cout << "Grad Mu" << endl;
			cout << myGradient.mCurrentGradMu << endl;
			cout << "Grad Mu Numerique" << endl;
			cout << myGradNum.mCurrentGradMu << endl;
			cout << "Grad Var" << endl;
			cout << myGradient.mCurrentGradVar << endl;
			cout << "Grad Var Numerique" << endl;
			cout << myGradNum.mCurrentGradVar << endl << endl;
			cout << "Grad Sigma" << endl;
			cout << myGradient.mCurrentGradSigma << endl;
			cout << "Grad Sigma Numerique" << endl;
			cout << myGradNum.mCurrentGradSigma << endl << endl;
			cout << "Grad Eps" << endl;
			cout << myGradient.mCurrentGradEps << endl;
			cout << "Grad Eps Numerique" << endl;
			cout << myGradNum.mCurrentGradEps << endl << endl;
			cout << "Grad Dens" << endl;
			cout << myGradient.mCurrentGradLogDens << endl;
			cout << "Grad Dens Numerique" << endl;
			cout << myGradNum.mCurrentGradLogDens << endl << endl;
			cout << "Diff LogDens" << endl;
			cout << myGradient.mCurrentDiffLogDensity << endl;
			cout << "Grad diff LogDens Numerique" << endl;
			cout << myGradNum.mCurrentDiffLogDensity << endl << endl;
			cout << "Hessien Mu" << endl;
			cout << myHessien.mCurrentHessMu << endl;
			cout << "Hessien numérique Mu" << endl;
			cout << myHessNum.mCurrentHessMu << endl;
			cout << "Hessien Var" << endl;
			cout << myHessien.mCurrentHessVar << endl;
			cout << "Hessien numérique Var" << endl;
			cout << myHessNum.mCurrentHessVar << endl;
			cout << "Hessien Sigma" << endl;
			cout << myHessien.mCurrentHessSigma << endl;
			cout << "Hessien numérique Sigma" << endl;
			cout << myHessNum.mCurrentHessSigma << endl;
			cout << "Hessien Eps" << endl;
			cout << myHessien.mCurrentHessEps << endl;
			cout << "Hessien numérique Eps" << endl;
			cout << myHessNum.mCurrentHessEps << endl;
			cout << "Grad Diff LogDens" << endl;
			cout << myHessien.mCurrentGradDiffLogDensity << endl;
			cout << "Grad Diff LogDens numérique" << endl;
			cout << myHessNum.mCurrentGradDiffLogDensity << endl;
			cout << "Hessien Dens" << endl;
			cout << myHessien.mCurrentHessDens << endl;
			cout << "Hessien numérique Dens" << endl;
			cout << myHessNum.mCurrentHessDens << endl;
		}
		myGradient.Update();
		myHessien.Update();
		myGradNum.Update();
		myHessNum.Update();
	}
	

	return 0;
}
