#pragma once 
#ifndef _CFUNCCONDMEANANDVAR_H_
#define _CFUNCCONDMEANANDVAR_H_

#include "RegArchDef.h"
#include "cRegArchValue.h"
#include "cRegArchGradient.h"
#include "cRegArchHessien.h"

/*!
	\file cAbstCondVar.h
	\brief Definition of the abstract class to implement conditional variance

	\author Jean-Baptiste DURAND, Ollivier TARAMASCO
	\date dec-18-2006 - last change feb-18-2011
*/

namespace RegArchLib {


	/*!
	 * \class cFuncMeanAndVar
	 * \brief class to implement functional representation of the derivatives of conditional mean and variance.
	 */
	class cFuncMeanAndVar
	{
	public :
		uint mNtheta;
		uint mNu;
		uint mNh;
		bool mCondMeanFuncBool;
		cDVector mdFx;
		cDVector mdFu;
		cDVector mdFh;
		cDMatrix md2Fxx;
		cDMatrix md2Fxu;
		cDMatrix md2Fuu;
		cDMatrix md2Fxh;
		cDMatrix md2Fuh;
		cDMatrix md2Fhh;

	public :
		cFuncMeanAndVar();
		cFuncMeanAndVar(uint theNtheta, bool theCondMeanBool = true, uint theNu = 0, uint theNh = 0);
		virtual ~cFuncMeanAndVar();
		void ComputeGradFunc(cRegArchGradient& theGradData);
		void ComputeHessFunc(const cRegArchGradient& theGradData, cRegArchHessien& theHessData);
		void ComputeGradAndHessFunc(cRegArchGradient& theGradData, cRegArchHessien& theHessData);

	};


}
#endif
