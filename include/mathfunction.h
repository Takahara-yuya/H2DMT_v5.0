#pragma once
#ifndef MATH_FUNCTION_H_
#define MATH_FUNCTION_H_

#include "CommonHeaders.h"

double sign(double a)
{
	if (a > 0)
		return 1.0;
	else if (a < 0)
		return -1.0;
	else
		return 0.0;
}


#endif