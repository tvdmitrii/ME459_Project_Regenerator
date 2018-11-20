/*!
*  \brief     Base class for LookupTable_1D and LookupTable_2D.
*  \details   Loads lookup tables stored in CSV format. Empty entries in the table are treated as 0. Each line MUST end with a comma ","
*  \author    Dmitrii Turygin 	https://github.com/tvdmitrii
*  \version   1.0
*  \date      1/1/2018
*/

#include "LookUpTable.h"
#include <iostream>
#include <sstream>

string LookupTable::getTitles()
{
	if (colNames.empty()) {
		return "No titles found!";
	}

	map<string, int>::iterator it = colNames.begin();
	string output = "";
	string colName;
	int colNum;

	while (it != colNames.end())
	{
		colName = it->first;
		colNum = it->second;

		output += "[" + colName + "]->" + to_string(colNum) + " ";
		it++;
	}

	return output;
}

LookupTable::LookupTable(string filepath)
{
	ifstream inputFile(filepath);

	if (!inputFile.good()) {
		throw invalid_argument("Path: '" + filepath + "' is not valid!");
	}

	vector<vector<double>> fields;
	if (inputFile) {
		
		string line;

		//Reads colomn titles
		getline(inputFile, headers);

		while (getline(inputFile, line)) {
			stringstream sep(line);
			string field;
			fields.push_back(vector<double>());
			while (getline(sep, field, ',')) {
				if (field.empty()) {
					field.assign("0");
				}
				fields.back().push_back(stod(field));
			}
		}
	}

	inputFile.close();

	size_t numRows = fields.size();
	size_t numColomns = fields.at(0).size();

	//Fills matrix_t<double> 'data' structure
	data = new util::matrix_t<double>(numRows, numColomns);
	int i = 0, j = 0;
	for (auto row : fields) {
		j = 0;
		for (auto field : row) {
			(*data)(i, j) = field;
			j++;
		}
		i++;
	}
}

LookupTable::~LookupTable()
{
	delete data;
}
