#pragma once
#ifndef _CARFIMA_H_
#define _CARFIMA_H_

#include "cAbstCondMean.h"
#include "cRegArchValue.h"
#include "cRegArchGradient.h"
#include "cRegArchHessien.h"

//#include <vector>
//#include <array>

/*!
\file cARFIMA.h
\brief Definition of the ARFIMA(p, d, q) class

\author Jean-Baptiste DURAND, Ollivier TARAMASCO
\date jan-26-2016 - Last change jan-26-2016
*/
namespace RegArchLib {

	/*!
	* \class cArfima
	* \brief Class to implement a ARFIMA(p, d, q) model
	*/
	class cArfima : public cAbstCondMean
	{
	private:
		cDVector mvAr; ///< Vector of AR coefficients. 
		cDVector mvMa; ///< Vector of MA coefficients.
		double mvFracD; ///< Value of the fractal d parameter
	private:
		uint mvNTruncLag; ///< Number of lags for AR representation
		cPolynome mvPolMa;
		cPolynome* mvGradPolMa;
		cPolynome** mvHessPolMa;

	public:
		cArfima(uint theNAr=0, uint theNMa=0, double theFracD=0, uint theNTruncLag = 50); ///< A simple constructor
		cArfima(const cDVector& theAr, const cDVector& theMa, double theFracD, uint theNTruncLag = 20); ///< Another constructor
		cArfima(const cArfima& theArfima); /// Recopy constructor
		virtual ~cArfima();  ///< aA simple destructor
		//cAbstCondMean* PtrCopy(void); /// < Return a copy of *this
		void Delete(void); ///< Delete mvAr
	#ifndef _RDLL_
		void Print(ostream& theOut = cout) const; ///< print the parameters
	#else
		void Print(void);
	#endif //_RDLL_
		void SetDefaultInitPoint(double theMean, double theVar);
		void SetDefaultInitPoint(cRegArchValue& theValue);
		void Set(const double theValue, const uint theIndex = 0, const uint theNumParam = 0); ///< Set model parameters.
		void Set(const cDVector& theVectParam, const uint theNumParam = 0); ///< Set model parameters.
		double Get(const uint theIndex, const uint theNumParam);
		cDVector& Get(const uint theNumParam);
		void ReAlloc(const uint theSize, const uint theNumParam = 0); ///< Allocation of the model parameters
		void ReAlloc(const cDVector& theVectParam, const uint theNumParam = 0); ///< Allocation of the model parameters
		cArfima& operator=(const cArfima& theSrc); ///< Standard affectation
		void ReAllocProxyMeanParameters(uint theOldNParam = 0);
		void UpdateProxyMeanParameters(void); /// Compute the truncate gradient 
		double ComputeMean(uint theDate, const cRegArchValue& theData) const; ///< Compute the conditional mean value
		uint  GetNAr(void) const;
		uint  GetNMa(void) const;
		uint GetNParam(void) const;
		uint GetNLags(void) const;
		void ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, uint theBegIndex);
		void ComputeHess(uint theDate, const cRegArchValue& theData, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex);
		void ComputeGradAndHess(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData, cRegArchHessien& theHessData, uint theBegIndex);
		void RegArchParamToVector(cDVector& theDestVect, uint theIndex);
		void VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex = 0);
#ifdef _DEBUG
		#ifndef _RDLL_
			void PolMaPrint(void) { cout << "mvPolMa = ";
									mvPolMa.Print(); 
									};
			void PolGradPrint(void) {
									uint myNParam = mvAr.GetSize() + mvMa.GetSize() + 1;
										for (uint i = 0; i < myNParam; i++)
										{
											cout << "mvGradPolMa[" << i << "]=";
											mvGradPolMa[i].Print();
										}
									};
		#endif //_RDLL_
#endif //_DEBUG
		void GetParamName(uint theIndex, char** theName);
		void GetParamName(uint theIndex, string theName[]);
	};

}
#endif // __CARFFIMA_H_

