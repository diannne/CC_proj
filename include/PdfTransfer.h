//===========================================================================
//===========================================================================
//===========================================================================
//==      PdfTransfer.h   ==   Author: Diana-Maria Popa       ==
//===========================================================================
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#ifndef __PDFTRANSFER__H__
#define __PDFTRANSFER__H__
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#include "KImage.h"
#include "MatrixModule.h"
#include <stdio.h>
#include <stdlib.h>
#include "Histogram.h"
#include <vector>
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
class PdfTransfer
{
public:
	PdfTransfer(unsigned int,unsigned int,unsigned int,unsigned int,int,int);
	//do not forget about freeing stuff !!!
	~PdfTransfer();
	bool NDPdfTransfer( float * , float * , float **  );
private:
	unsigned int _iterations;
	unsigned int _projections;
	unsigned int _inputSize;

	unsigned int _paletteSize;
	int _height;
	int _width;
	std::vector<float>  SinglePdfTransfer(Histogram * , Histogram * );
	MatrixModule *_matrixModule;
	void CreateHistograms(float *,  unsigned int , 
									float *, unsigned int  , unsigned int , 
									Histogram **, Histogram ** );
	std::vector<float> linearInterpolation(std::vector<float>,
											std::vector<float>,
											std::vector<float>);
	//to be deleted
	std::vector<float> oldLinearInterpolation(std::vector<float>,
											std::vector<float>,
											std::vector<float>);
};
#endif