// stdafx.cpp�: fichier source incluant simplement les fichiers Include standard
// VectorAndMatrix.pch repr�sente l'en-t�te pr�compil�
// stdafx.obj contient les informations de type pr�compil�es

#include "StdAfxVectorAndMatrix.h"

// TODO: faites r�f�rence aux en-t�tes suppl�mentaires n�cessaires dans STDAFX.H
// absents de ce fichier

#ifdef _LOG_
namespace VectorAndMatrixNameSpace
{
	using namespace std;

	cLogFile gLogFile;


	cLogFile::cLogFile()
	{
			mvLogFile = fopen("C:\\Users\\Ollivier\\AppData\\Local\\Temp\\RegArch23\\RegArch23.log", "w");
	}
		
	cLogFile::~cLogFile()
	{
		fclose(mvLogFile);
	}
		
	void cLogFile::WriteMessage(const char* theMessage)
	{
		fprintf(mvLogFile, "%s\n", theMessage);
	}
}
#endif // _LOG_