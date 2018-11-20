/*!
*  \brief     Used to look up a function of two variables.
*  \details   Base class is LookupTable. Loads lookup tables stored in CSV format, where first line needs to contain "ROWS_NAME/COLs_NAME".
*  \author    Dmitrii Turygin	https://github.com/tvdmitrii
*  \version   1.0
*  \date      1/1/2018
*/

#include "LookupTable_2D.h"

void LookupTable_2D::parseColNames() 
{
	string delimiter = "/";
	size_t pos = headers.find(delimiter);

	if (pos == string::npos) {
		colNames.insert(pair<string, int>(headers, 0));
		return;
	}

	string token = headers.substr(0, pos);
	colNames.insert(pair<string, int>(token, 0));
	
	headers.erase(0, pos + delimiter.length());
	headers.erase(headers.length()-1, headers.length());
	colNames.insert(pair<string, int>(headers, 1));
}

double LookupTable_2D::getValue(double rowValue, double colValue)
{
	return util::bilinear(rowValue, colValue, (*data));
}

LookupTable_2D::LookupTable_2D(string filepath) : LookupTable(filepath)
{	
	parseColNames();
}


LookupTable_2D::~LookupTable_2D()
{
}
