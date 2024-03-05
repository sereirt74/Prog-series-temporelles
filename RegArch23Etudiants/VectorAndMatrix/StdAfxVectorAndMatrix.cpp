// stdafx.cpp : fichier source incluant simplement les fichiers Include standard
// VectorAndMatrix.pch représente l'en-tête précompilé
// stdafx.obj contient les informations de type précompilées

#include "StdAfxVectorAndMatrix.h"

// TODO: faites référence aux en-têtes supplémentaires nécessaires dans STDAFX.H
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