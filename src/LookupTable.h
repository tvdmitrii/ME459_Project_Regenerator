#pragma once
#include <fstream>
#include <string>
#include <map>
#include "lib_util.h"

using namespace std;

/*!
*  \brief     Base class for LookupTable_1D and LookupTable_2D.
*  \details   Loads lookup tables stored in CSV format. Empty entries in the table are treated as 0. Each line MUST end with a comma ","
*  \author    Dmitrii Turygin 	https://github.com/tvdmitrii
*  \version   1.0
*  \date      1/1/2018
*/
class LookupTable
{
protected:

	//! Pointer to a matrix that stores lookup table entries.
	util::matrix_t<double>* data;

	//! Stores colomn names in [COLOMN_NAME]->COLOMN_NUMBER format.
	map<string, int> colNames;

	/*! \brief Puts colomn names in a map.
	*	Puts colomn names into colNames map in [COLOMN_NAME]->COLOMN_NUMBER format.
	*	\sa colNames
	*/
	virtual void parseColNames() = 0;

	string headers;

public:

	//! Returns a string containing colomn names in [COLOMN_NAME]->COLOMN_NUMBER format.
	string getTitles();

	/*! Constructor that takes path to lookup table CSV file
	*
	*	\param FilePath String containing absolute path to the CSV file
	*/
	LookupTable(string filepath);

	~LookupTable();
};

