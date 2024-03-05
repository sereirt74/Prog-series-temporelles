#pragma once 
#ifndef _CFIGARCH_H_
#define _CFIGARCH_H_

#include "cAbstCondVar.h"
#include "cRegArchValue.h"
#include "cRegArchGradient.h"
/*!
\file cFigarch.h
\brief Definition of the FIGARCH(p, d, q) class

\author Jean-Baptiste DURAND, Ollivier TARAMASCO
\date jan-26-2016 - Last change jan-26-2016
*/
namespace RegArchLib {

	/*!
	* \class cFiarch
	* \brief Class to implement a 'pure' FIGARCH(p, d, q) model
	*/
	class cFigarch : public cAbstCondVar
	{
	private:
		double mvConst; ///< Constant part of FIGARCH(p, q) variance model.
		cDVector mvArch; ///< Vector of ARCH coefficients. 
		cDVector mvGarch; ///< Vector of GARCH coefficients.
		double mvFracD; ///< Value of the fractal d parameter
	private :
		uint mvNTruncLag; ///< Number of lags for MA representation
		cPolynome mvPolMa;
		cPolynome* mvGradPolMa;
		cPolynome** mvHessPolMa;
	public:
		cFigarch(uint theNArch = 0, uint theNGarch = 0, double theFracD = 0, uint theNTruncLag=50); ///< A simple constructor
		cFigarch(double theConst, const cDVector& theArch, const cDVector& theGarch, double theFracD, uint theNTruncLag=20); ///< Another constructor
		cFigarch(const cFigarch& teFigarch);
		virtual ~cFigarch(); ///< A simple destructor
		//cAbstCondVar* PtrCopy() const; /// < Return a copy of *this				
		void Delete(void); /// Delete
#ifndef _RDLL_
		void Print(ostream& theOut = cout) const; ///< print the parameters
#else
		void Print(void);
#endif //_RDLL_
		void SetDefaultInitPoint(double theMean, double theVar);
		void SetDefaultInitPoint(cRegArchValue& theValue);
		void ReAlloc(const uint theSize, const uint theNumParam = 0); ///< Allocation of the model parameters
		void ReAlloc(const cDVector& theVectParam, const uint theNumParam = 0); ///< Allocation of the model parameters
		void Set(const double theValue, const uint theIndex = 0, const uint theNumParam = 0); ///< Set model parameters.
		void Set(const cDVector& theVectParam, const uint theNumParam = 0); ///< Set model parameters.
		double Get(const uint theIndex, const uint theNumParam);
		cDVector& Get(const uint theNumParam);
		cFigarch& operator =(const cFigarch& theSrc); ///< Standard affectation
		void ReAllocProxyVarParameters(uint theOldNParam=0) ;
		void UpdateProxyVarParameters(void);
		double ComputeVar(uint theDate, const cRegArchValue& theData) const;	///< Return conditional variance.
		uint GetNParam(void) const; ///< Number of parameters in that model part
		uint GetNArch(void) const; ///< Number of ARCH parameters
		uint GetNGarch(void) const; ///< Number of GARCH parameters
		uint GetNLags(void) const; ///< Number of past gradients required to compute gradient at current time t.
		uint GetNu(void) const;
		uint GetNh(void) const;
		void ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cAbstResiduals* theResiduals);
		void ComputeGradForM(uint theDate, const cRegArchValue& theValue, int theBegIndex, cFuncMeanAndVar& theDerivM, cAbstResiduals* theResids = NULL);
		void NumericComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, uint theBegIndex, cAbstResiduals* theResiduals, double theh=1e-4);
		void RegArchParamToVector(cDVector& theDestVect, uint theIndex = 0);
		void VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex = 0);
		void ComputeHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals);
		void ComputeHessForM(uint theDate, const cRegArchValue& theValue, uint theBegIndex, cFuncMeanAndVar& theDerivM, cAbstResiduals* theResids = NULL);
		void GetParamName(uint theIndex, char** theName);
		void GetParamName(uint theIndex, string theName[]);
#ifdef _DEBUG
		void PolMaPrint(void) { mvPolMa.Print(); };
		void PolGradPrint(void) {
									uint myNParam = mvArch.GetSize() + mvGarch.GetSize() + 2;
									for (uint i = 0; i < myNParam; i++) 
										mvGradPolMa[i].Print(); 
								};		
		void PolHessPrint(void) {
									uint myNParam = mvArch.GetSize() + mvGarch.GetSize() + 2;
									for (uint i = 0; i < myNParam; i++)
										for (uint j = i; j < myNParam; j++)
										{
										#ifndef _RDLL_	
											cout << "[" << i << "][" << j << "] : ";
										#else
											Rprintf("[%d][%d] : ", i, j);
										#endif // _RDLL_
											mvHessPolMa[i][j].Print();
										}
								};


#endif //_DEBUG
	};

}
#endif // _CFIGARCH_H_
