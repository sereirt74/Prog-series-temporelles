#pragma once

#ifndef _CTABOFSTRING_H_
#define _CTABOFSTRING_H_
#include <string>
using namespace std;

namespace RegArchLib {
	/*!
	 * \class cTabOfString
	 * \brief  Class to implement Tab of string
	 */
	class cTabOfString
	{
	public:
		int mSize;
		string* mTab;
	public: 
		cTabOfString(int theSize=0);
		virtual ~cTabOfString();
		void Delete(void);
		void ReAlloc(int theSize);

	};

}
#endif // _CTABOFSTRING_H_

