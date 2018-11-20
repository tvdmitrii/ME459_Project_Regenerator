#pragma once
#include "LookupTable.h"

/*!
*  \brief     Used to look up a functions of one variable.
*  \details   Base class is LookupTable. Loads lookup tables stored in CSV format, where first line needs to contain colomn names separated with ','. Example: "Re","f","j",
*				Each line MUST end with a comma ","
*  \author    Dmitrii Turygin	https://github.com/tvdmitrii
*  \version   1.1
*  \date      1/9/2018
*/
class LookupTable_1D : public LookupTable
{
protected:

	/*! \brief Puts colomn names in a map.
	*	Puts colomn names into colNames map in [COLOMN_NAME]->COLOMN_NUMBER format.
	*	\sa colNames
	*/
	virtual void parseColNames() override;
public:

	/*! \brief Performs linear interpolation to retrieve the value.
	*
	*	\param colY Name of the Y colomn
	*	\param colX Name of the X colomn
	*	\param x3 Interpolation value (x)
	*	\return Interpolated value.
	*/
	double getValue(string colY, string colX, double x3);

	/*! Uses constructor of the base class LookupTable.*/
	LookupTable_1D(string filepath);
	~LookupTable_1D();
};

