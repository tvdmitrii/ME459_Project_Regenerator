#pragma once
#include "LookupTable.h"

/*!
*  \brief     Used to look up a function of two variables.
*  \details   Base class is LookupTable. Loads lookup tables stored in CSV format, where first line needs to contain "ROWS_NAME/COLs_NAME".
*  \author    Dmitrii Turygin	https://github.com/tvdmitrii
*  \version   1.0
*  \date      1/1/2018
*/
class LookupTable_2D : public LookupTable
{
protected:

	/*! \brief Puts colomn names in a map.
	*	Puts colomn names into colNames map in [COLOMN_NAME]->COLOMN_NUMBER format.
	*	\sa colNames
	*/
	virtual void parseColNames() override;
public:

	/*! \brief Performs bilinear interpolation to retrieve the value.
	*
	*	\param rowValue Value of a variable that which entries go verticaly
	*	\param colomnValue Value of a variable that which entries go horizontaly
	*	\return Interpolated value from the lookup table
	*	\sa linearInterpolation()
	*/
	double getValue(double rowValue, double colomnValue);

	/*! Uses constructor of the base class LookupTable.*/
	LookupTable_2D(string filepath);
	~LookupTable_2D();
};

