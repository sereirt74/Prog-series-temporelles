#include "StdAfxVectorAndMatrix.h"
/*!
	\file cGSLMatrix.cpp
	\brief cGSLMatrix class functions definitions.
	\author Ollivier TARAMASCO 
	\date oct-13-2008
*/
namespace VectorAndMatrixNameSpace
{

	/*!
		\fn cGSLMatrix::cGSLMatrix
		\param uint theNRow, uint theNCol, double theVal
		\brief standard constructor: initializes cGSLMatrix size to theSize and affects theVal to all components
	*/
	cGSLMatrix::cGSLMatrix(int theNRow, int theNCol, double theVal)
	{
		if ( (theNRow == 0) && (theNCol == 0) )
			mvMat = (gsl_matrix *)NULL ;
		else
		{	if ( (theNRow > 0) && (theNCol > 0) )
			{	mvMat = gsl_matrix_alloc(theNRow, theNCol) ;	
				gsl_matrix_set_all(mvMat, theVal) ;
			}
			else
				throw cError("NRow and NCol must be strictly positive") ;
		}
	}

	/*!
		\fn cGSLMatrix::cGSLMatrix
		\param const cGSLVector &theVect
		\brief Constructor: creates a one column matrix = theVect
	*/
	cGSLMatrix::cGSLMatrix(const cGSLVector& theVect) 
	{
	uint myNRow = theVect.GetSize() ;
		if (myNRow != 0)
		{	mvMat = gsl_matrix_alloc(theVect.GetSize(), 1) ;
			for (int i = 0 ; i < (int)myNRow ; i++)
				gsl_matrix_set(mvMat, i, 0, theVect[i]) ;
		}
		else
			mvMat = NULL ;
	}

	/*!
		\fn cGSLMatrix::cGSLMatrix
		\param const cGSLMatrix& theMat
		\brief Constructor: creates and affects theMat to *this
	*/
	cGSLMatrix::cGSLMatrix(const cGSLMatrix& theMat)
	{
	int myNRow = (int)theMat.GetNRow() ;
	int myNCol = (int)theMat.GetNCol() ;
		if ( (myNRow > 0) && (myNCol > 0) )
		{	mvMat = gsl_matrix_alloc(myNRow, myNCol) ;
	
			for (int i = 0 ; i < myNRow ; i++)
				for (int j = 0 ; j < myNCol ; j++)
					gsl_matrix_set(mvMat, i, j, theMat[i][j]) ;
		}
		else
			mvMat = NULL ;
	}

	/*!
		\fn cGSLMatrix::~cGSLMatrix
		\param none
		\brief standard destructor
	*/
	cGSLMatrix::~cGSLMatrix()
	{
		if (mvMat != NULL)
		{	gsl_matrix_free(mvMat) ;
			mvMat = NULL ;
		}
	}

	/*!
		\fn void cGSLMatrix::Delete
		\param void
		\brief destructor
	*/
	void cGSLMatrix::Delete(void)
	{
		if (mvMat != NULL)
		{	gsl_matrix_free(mvMat) ;
			mvMat = NULL ;
		}
	}

	/*!
		\fn void cGSLMatrix::ReAlloc
		\param uint theNRow, uint theNCol, double theVal
		\brief Reallocates the matrix to a theNRow x theNCol matrix and affects theVal to all matrix components
	*/
	void cGSLMatrix::ReAlloc(uint theNRow, uint theNCol, double theVal)
	{
		Delete() ;
		if ( (theNRow > 0) && (theNCol > 0) )
		{	mvMat = gsl_matrix_alloc(theNRow, theNCol) ;	
			gsl_matrix_set_all(mvMat, theVal) ;
		}
		else
			throw cError("NRow and NCol must be strictly positive") ;
	
	}

	/*!
		\fn void cGSLMatrix::ReAlloc
		\param const cGSLVector& theVect
		\brief Reallocates the matrix to a one colum matrix and affects theVect to this matrix
	*/
	void cGSLMatrix::ReAlloc(const cGSLVector& theVect)
	{
	int myNRow = (int)GetNRow() ;
	int myNCol = (int)GetNCol() ;
		if ( (myNRow != (int)theVect.GetSize()) || (myNCol != 1) )
		{	myNRow = (int)theVect.GetSize() ;
			Delete() ;
			mvMat = gsl_matrix_alloc(myNRow, 1) ;
		}
		for (int i = 0 ; i < myNRow ; i++)
			gsl_matrix_set(mvMat, i, 0, theVect[i]) ;
	}

	/*!
		\fn void cGSLMatrix::ReAlloc
		\param const cGSLMatrix& theMat
		\brief Reallocates the matrix and *this=theMat
	*/
	void cGSLMatrix::ReAlloc(const cGSLMatrix& theMat)
	{
	uint myNRow = GetNRow() ;
	uint myNCol = GetNCol() ;
		if ( (myNRow != theMat.GetNRow()) || (myNCol != theMat.GetNCol()) )
		{	myNRow = theMat.GetNRow() ;
			myNCol = theMat.GetNCol() ;
			Delete() ;
			mvMat = gsl_matrix_alloc(myNRow, myNCol) ;
		}
		for (uint i = 0 ; i < myNRow ; i++)
			for (uint j = 0 ; j < myNCol ; j++)
				gsl_matrix_set(mvMat, i, j, theMat[i][j]) ;
	}

	/*!
		\fn double* cGSLMatrix::operator []
		\param uint theNRow
		\brief Returns the adress of theMat[theNRow]. 
	*/
	double* cGSLMatrix::operator [](int theNRow) const
	{
		if (theNRow < (int)GetNRow())
		{	
//	#ifdef _RDLL_
//			return (double *)&(mvMat->mData[theNRow * mvMat->mNCol]) ;
//	#else
			return (double *)&(mvMat->data[theNRow * mvMat->tda]) ;
//	#endif // _RDLL_
		}
		else
			throw cError("bad index") ;
	}

	void cGSLMatrix::Set(double theValue, int theRow, int theCol)
	{
		if ( (theRow < (int)GetNRow()) && (theCol < (int)GetNCol()))
			gsl_matrix_set(mvMat, theRow, theCol, theValue);
		else
			throw cError("bad index");
	}

	/*!
		\fn uint cGSLMatrix::GetNRow
		\param void
		\brief Returns the number of rows of the matrix
	*/
	uint cGSLMatrix::GetNRow(void) const
	{	if (mvMat == NULL)
			return 0 ;
		else
//		#ifdef _RDLL_
//			return mvMat->mNRow ;
//		#else
			return (uint)(mvMat->size1) ;
//		#endif // _RDLL_
	}

	/*!
		\fn uint cGSLMatrix::GetNCol
		\param void
		\brief Returns the number of columns of the matrix
	*/
	uint cGSLMatrix::GetNCol(void) const
	{ 	if (mvMat == NULL)
			return 0 ;
		else
//		#ifdef _RDLL_
//			return mvMat->mNCol ;
//		#else
			return (uint)(mvMat->size2) ;
//		#endif // _RDLL_
	}

	/*!
		\fn cGSLMatrix& cGSLMatrix::operator =
		\param const cGSLMatrix& theSrcMatrix
		\brief *this=theSrcMatrix
	*/
	cGSLMatrix& cGSLMatrix::operator =(const cGSLMatrix& theSrcMatrix)
	{
		Delete();
		ReAlloc(theSrcMatrix);
		return *this ;
	}

	/*!
		\fn cGSLMatrix& cGSLMatrix::operator =
		\param const cGSLVector& theVect
		\brief *this=one column matrix=theVect
	*/
	cGSLMatrix& cGSLMatrix::operator =(const cGSLVector& theVect)
	{
		Delete();
		ReAlloc(theVect);
		return *this ;
	}

	/*!
		\fn cGSLMatrix& cGSLMatrix::operator =
		\param theVal
		\brief *this[i][j] = theVal for all i,j 
	*/
	cGSLMatrix& cGSLMatrix::operator =(double theVal)
	{
		if ( mvMat != NULL )
		{	gsl_matrix_set_all(mvMat, theVal) ;
		}
		return *this ;
	}

	/*!
		\fn cGSLMatrix& cGSLMatrix::operator +=
		\param const cGSLMatrix& theMatrix
		\brief *this[i][j] += theMatrix[i][j] for all i,j
	*/
	cGSLMatrix& cGSLMatrix::operator +=(const cGSLMatrix& theMatrix)
	{
		*this = *this + theMatrix ;
		return *this ;
	}

	/*!
		\fn cGSLMatrix& cGSLMatrix::operator +=
		\param double theVal
		\brief *this[i][j] += theVal for all i, j
	*/
	cGSLMatrix& cGSLMatrix::operator +=(double theVal)
	{
		*this = *this + theVal ;
		return *this ;
	}
	
	/*!
		\fn cGSLMatrix& cGSLMatrix::operator -=
		\param const cGSLMatrix& theMatrix
		\brief *this[i][j] += theMatrix[i][j] for all i, j

		*/
	cGSLMatrix& cGSLMatrix::operator -=(const cGSLMatrix& theMatrix)
	{
		*this = *this -theMatrix ;
		return *this ;
	}
	
	/*!
		\fn cGSLMatrix& cGSLMatrix::operator -=
		\param double theVal
		\brief *this[i][j] += theVal for all i, j
	*/
	cGSLMatrix& cGSLMatrix::operator -=(double theVal)
	{
		*this = *this - theVal ;
		return *this ;
	}

	/*!
		\fn cGSLMatrix& cGSLMatrix::operator *=
		\param double theMatrix
		\brief *this *= theMatrix. Beware to matrices sizes
	*/
	cGSLMatrix& cGSLMatrix::operator *=(const cGSLMatrix& theMatrix)
	{	
		*this = *this * theMatrix ;
		return *this ;
	}

	/*!
		\fn cGSLMatrix& cGSLMatrix::operator *=
		\param double theVal
		\brief *this[i][j] *= theVal for all i, j
	*/
	cGSLMatrix& cGSLMatrix::operator *=(double theVal)
	{
		*this = *this * theVal ;
		return *this ;
	}

	/*!
		\fn cGSLMatrix& cGSLMatrix::operator /=
		\param double theVal
		\brief **this[i][j] /= theVal for all i, j
	*/
	cGSLMatrix& cGSLMatrix::operator /=(double theVal)
	{
		*this = *this / theVal ;
		return *this ;
	}

#ifndef _RDLL_
/*!
		\fn void cGSLMatrix::Print(std::ostream& theStream)
		\param std::ostream& theStream (default cout) 
		\brief Pront the matrix
	*/
	void cGSLMatrix::Print(std::ostream& theStream)
	{	theStream << *this ;
	}
#else
	void cGSLMatrix::Print(void)
	{
	uint myNCol = GetNCol() - 1;
	uint myNRow = GetNRow();
		for (uint i = 0; i < myNRow; i++)
		{
			for (uint j = 0; j < myNCol; j++)
				Rprintf("%f\t", (*this)[i][j]);
			Rprintf("%f\n", (*this)[i][myNCol]);
		}
	}
#endif // _RDLL_
	void cGSLMatrix::SetRow(uint theRowNumber, const cGSLVector& theVector)
	{
		if ((theRowNumber < GetNRow()) && (theVector.GetSize() == GetNCol()))
			for (uint j = 0; j < GetNCol(); j++)
				gsl_matrix_set(mvMat, theRowNumber, j, theVector[j]);
		else
			throw cError("Worng dimension in cGSLMatrix::SetRow");
	}
	
	void cGSLMatrix::SetColumn(uint theColumnNumber, const cGSLVector& theVector)
	{
		if ((theColumnNumber < GetNCol()) && (theVector.GetSize() == GetNRow()))
			for (uint i = 0; i < GetNRow(); i++)
				gsl_matrix_set(mvMat, i, theColumnNumber, theVector[i]);
		else
			throw cError("Worng dimension in cGSLMatrix::SetColumn");
	}

	/*!
		\fn cGSLMatrix operator+
		\param const cGSLMatrix& theLeft, const cGSLMatrix& theRight
		\brief returns matrix operation: theLeft+theRight
	*/
	cGSLMatrix operator+(const cGSLMatrix& theLeft, const cGSLMatrix& theRight) 
	{
	uint myNCol = theLeft.GetNCol() ;
	uint myNRow = theLeft.GetNRow() ;
	cGSLMatrix* myTmpMat = new cGSLMatrix(theLeft);

		if ( (myNRow == theRight.GetNRow()) && (myNCol == theRight.GetNCol()) )
		{	for (uint i = 0 ; i < myNRow ; i++)
				for (uint j = 0 ; j < myNCol ; j++)
					(*myTmpMat)[i][j] += theRight[i][j] ;
			return *myTmpMat ;
		}
		else
			throw cError("wrong size") ;
	}

	/*!
		\fn cGSLMatrix operator+
		\param const cGSLMatrix& theLeft, double theVal
		\brief returns matrix theLeft[i][j] + theVal for all i,j
	*/
	cGSLMatrix operator+(const cGSLMatrix& theLeft, double theVal) 
	{
	uint myNCol = theLeft.GetNCol() ;
	uint myNRow = theLeft.GetNRow() ;
	cGSLMatrix* myTmpMat = new cGSLMatrix(theLeft);

		for (uint i = 0 ; i < myNRow ; i++)
			for (uint j = 0 ; j < myNCol ; j++)
				(*myTmpMat)[i][j] += theVal ;
		return *myTmpMat ;
	}

	/*!
		\fn cGSLMatrix operator+
		\param double theVal, const cGSLMatrix& theRight
		\brief returns matrix theVal + theRight[i][j] for all i,j
	*/
	cGSLMatrix operator+(double theVal, const cGSLMatrix& theRight) 
	{
		return theRight + theVal ;
	}

	/*!
		\fn cGSLMatrix operator-
		\param const cGSLMatrix& theLeft, const cGSLMatrix& theRight
		\brief returns matrix theLeft[i][j]-theRight[i][j] for all i,j
	*/
	cGSLMatrix operator-(const cGSLMatrix& theLeft, const cGSLMatrix& theRight) 
	{
	uint myNCol = theLeft.GetNCol() ;
	uint myNRow = theLeft.GetNRow() ;
	cGSLMatrix* myTmpMat = new cGSLMatrix(theLeft);

		if ( (myNRow == theRight.GetNRow()) && (myNCol == theRight.GetNCol()) )
		{	for (uint i = 0 ; i < myNRow ; i++)
				for (uint j = 0 ; j < myNCol ; j++)
					(*myTmpMat)[i][j] -= theRight[i][j] ;
			return *myTmpMat ;
		}
		else
			throw cError("wrong size") ;
	}

	/*!
		\fn cGSLMatrix operator-
		\param const cGSLMatrix& theLeft, const cGSLMatrix& theRight
		\brief returns matrix theLeft[i][j]-theRight[i][j] for all i,j
	*/
	cGSLMatrix operator-(const cGSLMatrix& theLeft, double theVal) 
	{
	uint myNCol = theLeft.GetNCol() ;
	uint myNRow = theLeft.GetNRow() ;
	cGSLMatrix* myTmpMat = new cGSLMatrix(theLeft);

		for (uint i = 0 ; i < myNRow ; i++)
				for (uint j = 0 ; j < myNCol ; j++)
					(*myTmpMat)[i][j] -= theVal ;
			return *myTmpMat ;
	}

	/*!
		\fn cGSLMatrix operator-
		\param double theVal, const cGSLMatrix& theRight
		\brief returns matrix theVal-theRight[i][j] for all i,j
	*/
	cGSLMatrix operator-(double theVal, const cGSLMatrix& theRight) 
	{
	int myNCol = (int)theRight.GetNCol() ;
	int myNRow = (int)theRight.GetNRow() ;
	cGSLMatrix* myTmpMat = new cGSLMatrix(myNRow, myNCol, theVal);
		for (int i = 0 ; i < myNRow ; i++)
			for (int j = 0 ; j < myNCol ; j++)
				(*myTmpMat)[i][j] -= theRight[i][j] ;

		return *myTmpMat ;
	}

	/*!
		\fn cGSLMatrix operator*
		\param const cGSLMatrix& theLeft, const cGSLMatrix &theRight
		\brief returns matrix multiplication theLeft * theRight
	*/
	cGSLMatrix operator *(const cGSLMatrix& theLeft, const cGSLMatrix &theRight)
	{	
		cGSLMatrix* myTmpMat = new cGSLMatrix(theLeft.GetNRow(), theRight.GetNCol(), 0.0L);
	
		if ( (theLeft.GetNCol() == theRight.GetNRow()) )
		{	for (uint i = 0 ; i < theLeft.GetNRow() ; i++)
				for (uint j = 0 ; j < theRight.GetNCol() ; j++)
					for (uint k = 0 ; k < theLeft.GetNCol() ; k++)
						(*myTmpMat)[i][j] += theLeft[i][k] * theRight[k][j] ;
			return *myTmpMat ;
		}
		else
			throw cError("wrong matrices size") ;
	}

	/*!
		\fn cGSLMatrix operator*
		\param const cGSLMatrix& theLeft, const cGSLVector& theVect
		\brief returns vector theLeft * theVect
	*/
	cGSLVector operator *(const cGSLMatrix& theLeft, const cGSLVector& theVect)
	{
	uint myNCol = theLeft.GetNCol() ;
	uint myNRow = theLeft.GetNRow() ;
		if (myNCol == theVect.GetSize())
		{	
		cGSLVector* myTmpVect = new cGSLVector(myNRow, 0.0L);
			for (uint i = 0 ; i < myNRow ; i++)
				for (uint j = 0 ; j < myNCol ; j++)
					(*myTmpVect)[i] += theLeft[i][j] * theVect[j] ;
			return *myTmpVect ;
		}
		else
			throw cError("wrong matrix or vector size") ;
	}

	/*!
		\fn cGSLMatrix operator*
		\param const cGSLVector& theVect, const cGSLMatrix& theRight
		\brief returns  matrix theVect * theRight
	*/
	cGSLMatrix operator *(const cGSLVector& theVect, const cGSLMatrix& theRight)
	{
	uint myNCol = theRight.GetNCol() ;
	uint myNRow = theVect.GetSize() ;
		if  (theRight.GetNRow() == 1)
		{	
		cGSLMatrix* myTmpMat = new cGSLMatrix(myNRow, myNCol);
			for (uint i = 0 ; i < myNRow ; i++)
				for (uint j = 0 ; j < myNCol ; j++)
					(*myTmpMat)[i][j] += theVect[i]*theRight[0][j] ;
			return *myTmpMat;
		}
		else
			throw cError("wrong matrix or vector size") ;
	}

	/*!
		\fn cGSLMatrix operator*
		\param const cGSLMatrix& theMatrix, double theLambda
		\brief returns matrix theLambda * theMatrix
	*/
	cGSLMatrix operator *(const cGSLMatrix& theMatrix, double theLambda)
	{	
	cGSLMatrix* myTmpMat = new cGSLMatrix(theMatrix.GetNRow(), theMatrix.GetNCol()) ;
		for (uint i = 0 ; i < theMatrix.GetNRow() ; i++)
			for (uint j = 0 ; j< theMatrix.GetNCol() ; j++)
				(*myTmpMat)[i][j] = theLambda*theMatrix[i][j] ;
		return *myTmpMat;
	}

	/*!
		\fn cGSLMatrix operator*
		\param double theLambda, const cGSLMatrix& theMatrix
		\brief returns matrix theLambda * theMatrix
	*/
	cGSLMatrix operator *(double theLambda, const cGSLMatrix& theMatrix)
	{
		cGSLMatrix* myTmpMat = new cGSLMatrix(theMatrix.GetNRow(), theMatrix.GetNCol());
		for (uint i = 0; i < theMatrix.GetNRow(); i++)
			for (uint j = 0; j < theMatrix.GetNCol(); j++)
				(*myTmpMat)[i][j] = theLambda * theMatrix[i][j];
		return *myTmpMat;
	}

	/*!
		\fn cGSLMatrix operator/
		\param const cGSLMatrix& theMatrix, double theLambda
		\brief returns matrix theMatrix/theLambda
	*/
	cGSLMatrix operator /(const cGSLMatrix& theMatrix, double theLambda)
	{	
		if (theLambda == 0.0L)
			throw cError("Division by zero") ;
	cGSLMatrix* myTmpMat = new cGSLMatrix(theMatrix);
		for (uint i = 0; i < theMatrix.GetNRow(); i++)
			for (uint j = 0; j < theMatrix.GetNCol(); j++)
				(*myTmpMat)[i][j] = theMatrix[i][j] / theLambda;
	return *myTmpMat ;
	}

#ifndef _RDLL_
	/*!
		\fn std::ostream& operator <<
		\param std::ostream& theStream, const cGSLMatrix& theMat
		\brief prints matrix theMat
	*/
	ostream& operator <<(ostream& theStream, const cGSLMatrix& theMat)
	{
	uint myNRow = theMat.GetNRow() ;
	uint myNColm1 = theMat.GetNCol() - 1 ;

		for (uint i = 0 ; i < myNRow ; i++)
		{	for (uint j = 0 ; j < myNColm1 ; j++)
				theStream << theMat[i][j] << "\t" ;
			theStream << theMat[i][myNColm1] << std::endl ;
		}
		return theStream ;
	}
#endif // _RDLL_

	/*!
		\fn void Svd
		\param const cGSLMatrix& theMatrix, cGSLMatrix& theU, cGSLVector& theS, cGSLMatrix& theV
		\brief singular value decomposition of theMatrix
	*/
	void Svd(const cGSLMatrix& theMatrix, cGSLMatrix& theU, cGSLVector& theS, cGSLMatrix& theV) 
	{
	int myNRow = theMatrix.GetNRow() ;
	int myNCol = theMatrix.GetNCol() ;

	gsl_matrix* myA = gsl_matrix_alloc(myNRow, myNCol) ;
	gsl_vector* myWork = gsl_vector_alloc(myNCol) ;
	gsl_matrix* myV = gsl_matrix_alloc(myNRow, myNCol) ;
	gsl_vector* myS = gsl_vector_alloc(myNCol) ;

		for (int i = 0 ; i < myNRow ; i++)
			for (int j = 0 ; j < myNCol ; j++)
				gsl_matrix_set(myA, i, j, theMatrix[i][j]) ;
		gsl_linalg_SV_decomp(myA, myV, myS, myWork) ; 
		theS.ReAlloc(myS) ;
		theU.ReAlloc(myNRow, myNCol) ;
		theV.ReAlloc(myNRow, myNCol) ;
		for (int i = 0 ; i < myNRow ; i++)
			for (int j = 0 ; j < myNCol ; j++)
			{	theU[i][j] = gsl_matrix_get(myA, i, j) ;
				theV[i][j] = gsl_matrix_get(myV, i, j) ;
			}
		gsl_matrix_free(myA) ;
		gsl_matrix_free(myV) ;
		gsl_vector_free(myWork) ;
		gsl_vector_free(myS);
	}

	cGSLMatrix Diag(cGSLVector& theVect)
	{
	uint myN = theVect.GetSize();
	cGSLMatrix* myTmpMat = new cGSLMatrix(myN, myN, 0.0);
		for (uint i = 0; i < myN; i++)
			(*myTmpMat)[i][i] = theVect[i];
		return*myTmpMat;

	}

	cGSLMatrix Diag(cGSLMatrix& theMat)
	{
	uint myN = theMat.GetNCol();
	cGSLVector* myTmpVect = new cGSLVector(myN);
		
		for (uint i = 0; i < myN; i++)
			(*myTmpVect)[i] = theMat[i][i];
		return *myTmpVect;
	}

	cGSLMatrix Transpose(const cGSLMatrix& theMatrix)
	{
	uint myN = theMatrix.GetNRow(), myP = theMatrix.GetNCol();
	cGSLMatrix* myTmpMat = new cGSLMatrix(myP, myN);
		
		for (uint i = 0; i < myN; i++)
			for (uint j = 0; j < myP; j++)
				(*myTmpMat)[j][i] = theMatrix[i][j];
		return *myTmpMat;
	}

	cGSLMatrix Transpose(const cGSLVector& theVector)
	{
	uint myN = theVector.GetSize();
	cGSLMatrix* myTmpMat = new cGSLMatrix(1, myN);

		for (uint i = 0; i < myN; i++)
			(*myTmpMat)[0][i] = theVector[i];
		return *myTmpMat;
	}

	/*!
		\fn cGSLMatrix Inv
		\param const cGSLMatrix& theMatrix
		\brief Returns the inverse of theMatrix
	*/
	cGSLMatrix Inv(const cGSLMatrix& theMatrix)
	{
	int myNRow = theMatrix.GetNRow() ;
	int myNCol = theMatrix.GetNCol() ;
		if (myNRow != myNCol)
			throw cError("not a square matrix") ;
		if (myNRow == 0)
			throw cError("wrong size") ;

	cGSLMatrix myU, myV ;
	cGSLVector myS ;
		Svd(theMatrix, myU, myS, myV) ;

		for (int i = 0 ; i < myNCol ; i++)
			if (myS[i] == 0)
				throw cError("Not inversible matrix") ;
			else
				myS[i] = 1.0/myS[i] ;
	cGSLMatrix mySInvMat = Diag(myS);
	cGSLMatrix* myTmpMat = new cGSLMatrix(myV);
//		myTmpMat = myV * mySInvMat * Transpose(myU);
		*myTmpMat *= mySInvMat;
		*myTmpMat *= Transpose(myU);
	
		return *myTmpMat ;
	}

	cGSLMatrix Inv(const cGSLMatrix& theMatrix, int& theError)
	{
		theError = 0;
		int myNRow = theMatrix.GetNRow();
		int myNCol = theMatrix.GetNCol();
		if (myNRow != myNCol)
			throw cError("not a square matrix");
		if (myNRow == 0)
			throw cError("wrong size");

		cGSLMatrix myU, myV;
		cGSLVector myS;
		Svd(theMatrix, myU, myS, myV);

		for (int i = 0; i < myNCol; i++)
			if (myS[i] == 0)
			{
				theError = 1;
			}
			else
				myS[i] = 1.0 / myS[i];
		cGSLMatrix mySInvMat = Diag(myS);

		cGSLMatrix* myTmpMat = new cGSLMatrix(myV);
		//		myTmpMat = myV * mySInvMat * Transpose(myU);
		*myTmpMat *= mySInvMat;
		*myTmpMat *= Transpose(myU);
		return *myTmpMat;
	}

	/*!
		\fn void ClearMatrix
		\param cGSLMatrix& theMatrix
		\brief Deletes matrix theMatrix
	*/
	void ClearMatrix(cGSLMatrix& theMatrix)
	{
		theMatrix.Delete() ;
	}

	cGSLMatrix Abs(const cGSLMatrix& theMat)
	{
	cGSLMatrix* myTmpMat = new cGSLMatrix(theMat);
		for (uint i = 0; i < theMat.GetNRow(); i++)
			for (uint j = 0 ; j < theMat.GetNCol() ; j++)
				(*myTmpMat)[i][j] = fabs(theMat[i][j]);
		return *myTmpMat ;
	}

	cGSLMatrix Zeros(uint theN, uint theP)
	{
	cGSLMatrix* myTmpMat = new cGSLMatrix(theN, theP, 0.0L);
		return *myTmpMat;
	}

	cGSLMatrix Identity(uint theN)
	{
	cGSLMatrix* myTmpMat = new cGSLMatrix(theN, theN, 0.0L);

		for (uint i = 0; i < theN; i++)
			(*myTmpMat)[i][i] = 1.0;

		return *myTmpMat;
	}
	
	double Maxi(const cGSLMatrix& theMat)
	{
	double myMaxi = theMat[0][0];
		for (uint i = 0; i < theMat.GetNRow(); i++)
			for (uint j = 0; i < theMat.GetNCol(); i++)
				if (theMat[i][j] > myMaxi)
					myMaxi = theMat[i][j];
		return(myMaxi);
	}

	double Mini(const cGSLMatrix& theMat)
	{
	double myMini = theMat[0][0];
		for (uint i = 0; i < theMat.GetNRow(); i++)
			for (uint j = 0; i < theMat.GetNCol(); i++)
				if (theMat[i][j] < myMini)
					myMini = theMat[i][j];
		return(myMini);
	}

	bool IsNaN(const cGSLMatrix& theMatrix)
	{
		for (uint i = 0; i < theMatrix.GetNRow(); i++)
			for (uint j = 0; j < theMatrix.GetNCol(); j++)
				if (!isfinite(theMatrix[i][j]))
					return true;
		return false;
	}

} // namespace
