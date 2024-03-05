#include "StdAfxRegArchLib.h"

namespace RegArchLib {

	cTabOfString::cTabOfString(int theSize)
	{
		if (theSize > 0)
		{
			mSize = theSize;
			mTab = new string[theSize];
		}
		else
		{
			mSize = 0;
			mTab = NULL;
		}
	}

	cTabOfString::~cTabOfString()
	{
		if (mSize > 0)
		{
			delete[] mTab;
			mSize = 0;
			mTab = NULL;
		}
		else
		{
			mSize = 0;
			mTab = NULL;
		}
	}

	void cTabOfString::Delete(void)
	{
		if (mSize > 0)
		{
			delete[] mTab;
			mSize = 0;
			mTab = NULL;
		}
		else
		{
			mSize = 0;
			mTab = NULL;
		}
	}

	void cTabOfString::ReAlloc(int theSize)
	{
		Delete();
		if (theSize > 0)
			mTab = new string[theSize];
	}
}
