#include "StdAfxError.h"
/*! 
	\file cError.cpp
	\brief cError class functions declaration
*/
namespace ErrorNameSpace
{
	/*!
	 \fn cError
	 \param const string& theMess
	 \brief prints theMess and exits program
	*/
#ifndef _RDLL_
	string gMessage;
#endif //_RDLL_

#ifdef _RDLL_
	cError::cError(const char* theMess)
	{	error(theMess) ;
	}
#else
	cError::cError(const char* theMess)
	{
		if (theMess != NULL)
			printf("%s\n", theMess);
		abort() ;
	}

	cError::cError(const string& theMess)
	{
		if (theMess.length() != 0) 
			cout << theMess.c_str() << endl ;
		gMessage = theMess;
		abort();
	}

#endif // _RDLL_
} //namespace
