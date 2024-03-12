#include "StdAfxRegArchLib.h"
#include "cAparch.h"
#include <gsl/gsl_math.h>
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
        mvArch.ReAlloc(theNArch);
        mvGarch.ReAlloc(theNGarch);
        mvGamma.ReAlloc(theNArch);
        mvCste = 0.0;
        MESS_CREAT("cAparch");
    }

    cAparch::cAparch(const cAparch& theAparch)
        :cAbstCondVar(eAparch)  // call constructor of cAbstCondVar with type eAparch
    {
        *this = theAparch;
        MESS_CREAT("cAparch");
    }

    cAparch::~cAparch()
    {
        Delete();
        MESS_DESTR("cAparch");
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
        mvArch.Delete();
        mvGamma.Delete();
        mvGarch.Delete();
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
    }
#endif // _RDLL_

    void cAparch::SetDefaultInitPoint(double theMean, double theVar)
    {


    }

    void cAparch::SetDefaultInitPoint(cRegArchValue& theValue)
    {
    }

    void cAparch::ReAlloc(const uint theSize, const uint theNumParam)
    {
        switch (theNumParam)
        {
        case 1:
            mvArch.ReAlloc(theSize);
            break;
        case 2:
            mvGamma.ReAlloc(theSize);
            break;
        case 3:
            mvGarch.ReAlloc(theSize);
        default:
            //			throw cError("cAparch::ReAlloc - theNumParam must be in 1, 2, 3") ;
            break;
        }
    }

    void cAparch::ReAlloc(const cDVector& theVectParam, const uint theNumParam)
    {
        switch (theNumParam)
        {
        case 0: // mvConst
            if (theVectParam.GetSize() > 0)
                mvCste = theVectParam[0];
            else
                throw cError("cAparch::ReAlloc - Size of theVectParam must be > 0");
            break;
        case 1: // mvConst
            if (theVectParam.GetSize() > 0)
                mvDelta = theVectParam[0];
            else
                throw cError("cAparch::ReAlloc - Size of theVectParam must be > 0");
            break;
        case 2: // mvArch
            mvArch = theVectParam;
            break;
        case 3:
            mvGamma = theVectParam;
            break;
        case 4: // mvGarch
            mvGarch = theVectParam;
            break;
        default:
            throw cError("cAparch::ReAlloc - theNumParam must be in 0, 1, 2, 3, 4.");
            break;
        }
    }

    void cAparch::Set(const cDVector& theDVector, const uint thePlace)
    {
        switch (thePlace)
        {
        case 0:
            if (theDVector.GetSize() > 0)
                mvCste = theDVector[0];
            else
                throw cError("cAparch::Set - Size of theVectParam must be > 0");
            break;
        case 1:
            if (theDVector.GetSize() > 0)
                mvDelta = theDVector[0];
            else
                throw cError("cAparch::Set - Size of theVectParam must be > 0");
            break;
        case 2:
            mvArch = theDVector;
            break;
        case 3:
            mvGamma = theDVector;
            break;
        case 4:
            mvGarch = theDVector;
            break;
        default:
            throw cError("cAparch::Set - theNumParam must be in 0, 1, 2, 3, 4.");
            break;
        }
    }

    void cAparch::Set(const double theValue, const uint theIndex, const uint theNumParam)
    {
        switch (theNumParam)
        {
        case 0:
            mvCste = theValue;
            break;
        case 1:
            mvDelta = theValue;
            break;
        case 2:
            if (theIndex < mvArch.GetSize())
                mvArch[theIndex] = theValue;
            else
                throw cError("cAparch::Set - wrong index");
            break;
        case 3:
            if (theIndex < mvGamma.GetSize())
                mvGamma[theIndex] = theValue;
            else
                throw cError("cAparch::Set - wrong index");
            break;
        case 4:
            if (theIndex < mvGarch.GetSize())
                mvGarch[theIndex] = theValue;
            else
                throw cError("cAparch::Set - wrong index");
            break;
        default:
            throw cError("cAparch::Set - theNumParam must be in 0, 1, 2, 3, 4");
            break;
        }
    }

    double cAparch::Get(const uint theIndex, const uint theNumParam)
    {
        switch (theNumParam)
        {
        case 0:
            return mvCste;
            break;
        case 1:
            return mvDelta;
            break;
        case 2:
            return mvArch[theIndex];
            break;
        case 3:
            return mvGamma[theIndex];
            break;
        case 4:
            return mvGarch[theIndex];
            break;
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
        mvCste = theSrc.mvCste;
        mvDelta = theSrc.mvDelta;
        mvArch = theSrc.mvArch;
        mvGamma = theSrc.mvGamma;
        mvGarch = theSrc.mvGarch;
        return *this;
    }

    double cAparch::ComputeVar(uint theDate, const cRegArchValue& theValue) const
    {
        uint myp = mvArch.GetSize(),
            myq = mvGarch.GetSize();
        double myRes = mvCste;
        for (uint i = 1; i <= MIN(myp, theDate); i++)
            myRes += mvArch[i - 1] * pow((abs(theValue.mUt[theDate - i]) - mvGamma[i - 1] * theValue.mUt[theDate - i]), mvDelta);
        for (uint j = 1; j <= MIN(myq, theDate); j++)
            myRes += mvGarch[j - 1] * pow(theValue.mHt[theDate - j], mvDelta / 2.0);
        return myRes;
    }

    uint cAparch::GetNParam(void) const
    {
        return 2 + mvArch.GetSize() + mvGamma.GetSize() + mvGarch.GetSize();
    }

    uint cAparch::GetNLags(void) const
    {
        return MAX(mvArch.GetSize(), mvGarch.GetSize());
    }

    void cAparch::RegArchParamToVector(cDVector& theDestVect, uint theIndex)
    {
        uint mySize = GetNParam();
        if (theDestVect.GetSize() < mySize + theIndex)
            throw cError("Wrong size");
        theDestVect[theIndex] = mvCste;
        theDestVect[theIndex + 1] = mvDelta;
        mvArch.SetSubVectorWithThis(theDestVect, theIndex + 2);
        mvGamma.SetSubVectorWithThis(theDestVect, theIndex + 2 + mvArch.GetSize());
        mvGarch.SetSubVectorWithThis(theDestVect, theIndex + 2 + mvArch.GetSize() + mvGamma.GetSize());

    }

    void cAparch::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
    {
        uint mySize = theSrcVect.GetSize();
        if (GetNParam() + theIndex > mySize)
            throw cError("Wrong size");
        mvCste = theSrcVect[theIndex];
        mvDelta = theSrcVect[theIndex + 1];
        mvArch.SetThisWithSubVector(theSrcVect, theIndex + 2);
        mvGamma.SetThisWithSubVector(theSrcVect, theIndex + 2 + mvArch.GetSize());
        mvGarch.SetThisWithSubVector(theSrcVect, theIndex + 2 + mvArch.GetSize() + mvGamma.GetSize());

    }

    void cAparch::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, cAbstResiduals* theResiduals)
    {
        double condVar = cAparch::ComputeVar(theDate, theValue);
        uint myp = mvArch.GetSize(),
            myq = mvGarch.GetSize(),
            myBegIndex = theGradData.GetNMeanParam();
        theGradData.mCurrentGradVar = 0.0L;
        theGradData.mCurrentGradVar[myBegIndex] = 1.0;
        //ARCH and Gamma
        for (uint i = 1; i <= MIN(myp, theDate); i++)
            theGradData.mCurrentGradVar[myBegIndex + 1 + i] = pow((abs(theValue.mUt[theDate - i]) - mvGamma[i - 1] * theValue.mUt[theDate - i]), mvDelta);
        for (uint i = 1; i <= MIN(myp, theDate); i++)
            theGradData.mCurrentGradVar[myBegIndex + 1 + myp + i] = mvArch[i - 1] * mvDelta * (-theValue.mUt[theDate - i]) * pow((abs(theValue.mUt[theDate - i]) - mvGamma[i - 1] * theValue.mUt[theDate - i]), mvDelta - 1);
        for (uint i = 1; i <= MIN(myp, theDate); i++)
            theGradData.mCurrentGradVar -= mvDelta * mvArch[i - 1] * pow((abs(theValue.mUt[theDate - i]) - mvGamma[i - 1] * theValue.mUt[theDate - i]), mvDelta - 1) * (theGradData.mGradMt[i - 1] * (GSL_SIGN(theValue.mUt[theDate - i] - mvGamma[i - 1])));
        //GARCH
        for (uint j = 1; j <= MIN(myq, theDate); j++)
            theGradData.mCurrentGradVar[myBegIndex + 1 + 2 * myp + j] += pow(theValue.mHt[theDate - j], mvDelta / 2.0);
        for (uint j = 1; j <= MIN(myq, theDate); j++)
            theGradData.mCurrentGradVar += mvGarch[j - 1] * (mvDelta / 2.0) * pow(theValue.mHt[theDate - j], (mvDelta / 2.0) - 1) * theGradData.mGradHt[j - 1];
        ///
        theGradData.mCurrentGradVar *= 2.0 / mvDelta;
        //Delta
        theGradData.mCurrentGradVar[myBegIndex + 1] = 0;
        for (uint i = 1; i <= MIN(myp, theDate); i++) {
            double term = (abs(theValue.mUt[theDate - i]) - mvGamma[i - 1] * theValue.mUt[theDate - i]);
            theGradData.mCurrentGradVar[myBegIndex + 1] += mvArch[i - 1] * pow(term, mvDelta) * log(term) + (mvDelta / term) * (theGradData.mGradMt[i - 1][1] * (GSL_SIGN(theValue.mUt[theDate - i] - mvGamma[i - 1])));
        }
        for (uint j = 1; j <= MIN(myq, theDate); j++) {
            theGradData.mCurrentGradVar[myBegIndex + 1] += mvGarch[j - 1] * pow(theValue.mHt[theDate - j], mvDelta / 2.0) * (0.5 *
                log(theValue.mHt[theDate - j]) + (0.5 * mvDelta / pow(theValue.mHt[theDate - j], mvDelta / 2.0)) * theGradData.mGradHt[j - 1][1]);
        }
        theGradData.mCurrentGradVar[myBegIndex + 1] *= 2.0 / (mvDelta * pow(condVar, 0.5 * mvDelta));
        theGradData.mCurrentGradVar[myBegIndex + 1] -= (2.0 / (mvDelta * mvDelta)) * 0.5 * mvDelta * log(condVar);
        theGradData.mCurrentGradVar[myBegIndex + 1] *= condVar;
    }

    void cAparch::ComputeHess(uint theDate, const cRegArchValue& theValue, const cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
    {
    }

    void cAparch::ComputeGradAndHess(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData, cRegArchHessien& theHessData, cAbstResiduals* theResiduals)
    {

    }

    void cAparch::GetParamName(uint theIndex, char** theName)
    {

    }


    void cAparch::GetParamName(uint theIndex, string theName[])
    {

    }

}//namespace