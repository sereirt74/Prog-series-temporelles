// stdafx.h�: fichier Include pour les fichiers Include syst�me standard,
// ou les fichiers Include sp�cifiques aux projets qui sont utilis�s fr�quemment,
// et sont rarement modifi�s
//

#pragma once

#include <cmath>
#ifndef _RDLL_
	#include <iostream>
#endif //_RDLL_

#include "StdAfxError.h"
#include "StdAfxVectorAndMatrix.h"

#include "RegArchDef.h"
#include "cAbstCondMean.h"
#include "cAbstCondVar.h"
#include "cAbstResiduals.h"

#include "cPolynome.h"
#include "DerivativeTools.h"
#include "cNumericDerivative.h"

#include "cConst.h"
#include "cAr.h"
#include "cMa.h"
#include "cStdDevInMean.h"
#include "cVarInMean.h"
#include "cLinReg.h"
#include "cArfima.h"

#include "cConstCondVar.h"
#include "cArch.h"
#include "cGarch.h"
#include "cNgarch.h"
#include "cTarch.h"
#include "cEgarch.h"
#include "cAparch.h"
#include "cFigarch.h"
#include "cNgarch.h"
#include "cUgarch.h"
#include "cGtarch.h"
#include "cTsgarch.h"
#include "cNagarch.h"
#include "cSqrgarch.h"
#include "cStgarch.h"
#include "cLoggarch.h"

#include "SomeDistribution.h"
#include "cNormResiduals.h"
#include "cStudentResiduals.h"
#include "cGedResiduals.h"
#include "cMixNormResiduals.h"

#include "GslAndNloptOptim.h"

#include "cCondMean.h"
#include "RegArchCompute.h"
#include "cRegArchModel.h"
#include "cRegArchValue.h"
#include "cRegArchGradient.h"
#include "cRegArchHessien.h"

#include "cTabOfString.h"


#define WIN32_LEAN_AND_MEAN             // Exclure les en-t�tes Windows rarement utilis�s

