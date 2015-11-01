//===========================================================================
//===========================================================================
//===========================================================================
//==      Color2Gray.h   ==   Author: Diana-Maria Popa       ==
//===========================================================================
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#ifndef __COLOR2GRAY__H__
#define __COLOR2GRAY__H__
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#include "Color_Models_Conversions.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <cmath> //this is for M_PI
#include <vector>
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#define _CRT_SECURE_NO_WARNINGS
class Color2Gray
{
protected:
	int _iterations ;

public:

	enum  GREY_TYPE
	{
		BASIC_GREY = 0 ,
		GOOCH_GREY = 1 
	};

	//constructor
	Color2Gray( KImage *image )
	{
		_image = image;
		_iterations = 30;
	}

	~Color2Gray();
	

	KImage * localC2G( const char * strParams , GREY_TYPE Mode = GOOCH_GREY );

	void setMaxIter( unsigned int maxIter )
	{
		_iterations = maxIter;
	}

private:
	//variables
	double _theta;
	double _alpha;
	int _radius;
	//precompute these ones for better performance
	double cosTheta;
	double sinTheta;
	//pointer variables
	CIE_Lab * _data;
	double * _dataOutput;
	KImage *_image;
	KImage *_goochImage;
	KImage *_basicImage;
	double *_deltas ;
	
	KRGBColor * _greyOutput;
	KRGBColor * _basicOutput;

	GREY_TYPE _enumType;
	//private functions declaration
	bool getParams(const char * );
	
	bool convertToLab();

	double  crunch(double );

	double computeDelta( int  , int  );

	void  traverseCIELab();

	void computeOutput();

	void convertToRgb();

	void saveGreyImage();

	void setEnumType( GREY_TYPE  x )
	{
		_enumType = x;
	}
};
#endif