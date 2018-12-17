/*!
*  \brief     Used to look up a functions of one variable.
*  \details   Base class is LookupTable. Loads lookup tables stored in CSV format, where first line needs to contain colomn names separated with ','. Example: "Re","f","j",
*				Each line MUST end with a comma ","
*  \author    Dmitrii Turygin	https://github.com/tvdmitrii
*  \version   1.1
*  \date      1/9/2018
*/

#include "LookupTable_1D.h"

void LookupTable_1D::parseColNames()
{
	string delimiter = ",";

	size_t pos = 0;
	int colNumber = 0;
	string token;

	while ((pos = headers.find(delimiter)) != string::npos) {
		token = headers.substr(0, pos);

		colNames.insert(pair<string, int>(token, colNumber));

		headers.erase(0, pos + delimiter.length());

		colNumber++;
	}
}

double LookupTable_1D::getValue(string colY, string colX, double x3)
{
	if (colNames.empty()) {
		throw exception();
	}

	//Convert colomn names to colomn numbers using generated map 
	int colXNum = colNames.find(colX)->second;
	int colYNum = colNames.find(colY)->second;

	return util::linterp_col((*data), colXNum, x3, colYNum);
}

LookupTable_1D::LookupTable_1D(string filepath) : LookupTable(filepath)
{	
	parseColNames();
}


LookupTable_1D::~LookupTable_1D()
{
}
