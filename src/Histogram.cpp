//===========================================================================
//===========================================================================
//===========================================================================
//==     PdfTransfer.cpp   ==   Author: Diana-Maria POPA      ==
//===========================================================================
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#include "stdafx.h"
#include "Histogram.h"
//===========================================================================
//===========================================================================
Histogram::Histogram( ) :
_inferior(0.0),_superior(0.0),_counter (0),_globalCounter(0)
{
}
Histogram::Histogram( float limit1 , float limit2 ) :
_counter (0),_globalCounter(0)
{
	if(limit1 < limit2 )
	{
		_inferior = limit1;
		_superior = limit2;
	}
	else if (limit1 > limit2  )
	{
		_inferior = limit2;
		_superior = limit1;
	}
	else  
	{
		//do not know
	}
}
Histogram::~Histogram()
{
	//freeeee things
}
unsigned int Histogram::getCounter()
{
	return _counter;
}
unsigned int Histogram::getGlobalCounter()
{
	return _globalCounter;
}
float Histogram::GetInferiorLimit()
{
	return _inferior ;
}
float Histogram::GetSuperiorLimit()
{
	return _superior;
}
void Histogram::incrementCounter()
{
	++_counter;
}
void Histogram::incrementGlobalCounter()
{
	++_globalCounter;
}
void Histogram::SetInferiorLimit(float infLimit)
{
	_inferior = infLimit;
}
void Histogram::SetSuperiorLimit(float supLimit)
{
	_superior = supLimit;
}
void Histogram::SetOverallMax(float Max)
{
	_overallMax = Max;
}
void Histogram::SetOverallMin(float Min)
{
	_overallMin = Min;
}
float Histogram::GetOverallMax()
{
	return _overallMax;
}
float Histogram::GetOverallMin()
{
	return _overallMin;
}