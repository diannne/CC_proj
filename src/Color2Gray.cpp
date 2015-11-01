//===========================================================================
//===========================================================================
//===========================================================================
//==     Color2Gray.cpp   ==   Author: Diana-Maria POPA      ==
//===========================================================================
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#include "stdafx.h"
#include "Color2Gray.h"
#include <time.h>
#define _USE_MATH_DEFINES
#include <cmath> //this is for M_PI
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#pragma comment(lib, "../FreeImage/")
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================

Color2Gray::~Color2Gray()
{
	if( _data )
	{
		delete [] _data;
	}
	if(_dataOutput)
	{
		//code style for EA ???
		delete  [] _dataOutput;
	}
	if( _image ) 
	{
		delete [] _image ;
	}
	if( _goochImage ) 
	{
		delete [] _goochImage ;

	}
	if( _basicImage ) 
	{
		delete [] _basicImage ;
	}
	if( _deltas ) 
	{
		delete [] _deltas ;
	}
	if( _greyOutput ) 
	{
		delete [] _greyOutput ;
	}

}

//returns the greyscaled image
bool Color2Gray::getParams(const char * strParams)
{
	char * pch;
	pch = strtok ((char *)strParams," ,.-");
	int itr = 1;
	while (pch != NULL)
	{
		char *ptr;
		if( itr == 1 )
		{
			ptr = strstr(pch , "theta=");
			if( ptr  )
			{
				_theta = atof(ptr+6);
				//TODO use defines
				_theta = _theta * M_PI / 180.0;
				cosTheta = cos(_theta);
				sinTheta = sin(_theta);
			}
			else 
			{	
				printf("Expected parameters <Input_Image_File_Name (24BPP True-Color)> \
							theta=... alpha=... radius=... \n");
				return false;
			}

		}
		if( itr == 2 )
		{
			ptr = strstr(pch,"alpha=");
			if( ptr )
			{
				_alpha = atof(ptr+6);
			}
			else 
			{	
				printf("Expected parameters <Input_Image_File_Name (24BPP True-Color)> \
							theta=... alpha=... radius=... \n");
				return false;
			}

		}
		if( itr == 3 )
		{
			ptr = strstr(pch,"radius=") ;
			if( ptr )
			{
				_radius = atoi(ptr+7);
			}
			else 
			{	
				printf("Expected parameters <Input_Image_File_Name (24BPP True-Color)> \
							theta=... alpha=... radius=... \n");
				return false;
			}

		}
		pch = strtok (NULL, " ,.-");
		itr++;
	}
	return true;
}

bool Color2Gray::convertToLab() //this populates the _data
{
	if(!_image)
	{
		printf("There is no source image \n");
		return false;
	}
	//allocate the CIE_Lab
	//TODO: use multiple allocations(probl with big mem chunk)
	_data = new CIE_Lab[_image->GetWidth()*_image->GetHeight()];
	_image->BeginDirectAccess();
	int itr = 0  ;
	
	for( int i = _image->GetWidth() - 1 ; i >= 0 ; i-- )
	{
		for( int j = _image->GetHeight() - 1 ; j >= 0 ; j-- )
		{
			KRGBColor color ;
			_image->Get24BPPPixel(i,j,&color);
			CIE_Lab temp(color);
			_data[itr].a = temp.a;
			_data[itr].b = temp.b;
			_data[itr].L = temp.L;
			itr++;
		}
	}

	/*output data with a,b =0 to an image*/
	/*measure how much time it took*/
	_image->EndDirectAccess();

	return true;
}

double  Color2Gray::crunch(double x)
{
	return _alpha == 0 ? 0 : (_alpha*tanh(x/_alpha));
}

double Color2Gray::computeDelta( int i , int j )
{
	double dL = _data[i].L - _data[j].L;
	//TODO use a*a+b*b
	//double dC = sqrt( pow(_data[i].a-_data[j].a,2)  + pow(_data[i].b - _data[j].b,2) );
	double a = _data[i].a -_data[j].a;
	double b = _data[i].b - _data[j].b;
	double dC = sqrt( a*a  + b*b );
	//TODO if possible precalc crunch(dc)
	double crunchX = crunch(dC);
	if( fabs(dL) > crunchX )
	{
		return dL;
	}
	//TODO precalc cos and sin
	//double vTheta = cos(_theta)*(_data[i].a - _data[j].a) + sin(_theta)*(_data[i].b - _data[j].b);
	//double vTheta = cos(_theta)*a + sin(_theta)*b;
	double vTheta = cosTheta*a + sinTheta*b;
	if( (vTheta*crunchX) > 0 )
		return crunchX;
	return -crunchX;
}

void  Color2Gray::traverseCIELab()
{
	int N = _image->GetWidth() * _image->GetHeight();

	_deltas = new double[N];
	/*free me!!!*/
	
	//initialize to zeroed values
	for( int i = 0;i< N ; ++i )
		_deltas[i] = 0.0;

	if( _radius )
	{
		for( int x = _image->GetWidth() - 1 ; x >= 0 ; x-- )
		{
			for( int y = _image->GetHeight() - 1 ; y >= 0 ; y-- )
			{
				int i = x + y*_image->GetWidth(); //position of first relative pixel
				int x_start  = x - _radius;
				int x_end = x + _radius;
				int y_start = y - _radius;
				int y_end = y + _radius;
				for(int xx = x_start ; xx <= x_end ; xx++ )
				{
					if( xx >= _image->GetWidth() || xx < 0 )
						continue;
					for(int yy = y_start; yy <= y_end; yy++)
					{
						if( yy >= _image->GetHeight() || yy < 0 )
							continue;
						int j = xx + yy*_image->GetWidth();
						double little_delta = computeDelta(i,j);
						_deltas[i] += little_delta;
						_deltas[j]-= little_delta;
					}
				}
			}
		}
	}
	else
	{
		//TODO reverse for (x = N-2; x >= 0; x--)
		//for(int x = 0 ; x < N - 1 ; ++x )
		for(int x = N-1 ; x >= 0 ; x-- )
		{
			for( int y = x+1 ; y < N ; ++y )
			{
				double little_delta = computeDelta(x,y);
				_deltas[x] += little_delta;
				_deltas[y]-= little_delta;
			}
		}
	}
}

void Color2Gray::computeOutput()
{
	int N = _image->GetWidth() * _image->GetHeight();
	//allocate the CIE_Lab output 
	_dataOutput = new double[N];
	/*show intermediate progress based on iterations*/
	/*initializing the gray image to the luminance values*/
	for( int i = 0; i < N ; ++i )
		_dataOutput[i] = _data[i].L;
	if( _radius )
	{
		for (int iter = this->_iterations; iter != 0; iter--)
		{
			for( int x = _image->GetWidth() - 1 ; x >= 0 ; x-- )
			{
				for( int y = _image->GetHeight() - 1 ; y >= 0 ; y-- )
				{
					int i = x + y*_image->GetWidth(); //position of the relative pixel
					int x_start  = x - _radius;
					int x_end = x + _radius;
					int y_start = y - _radius;
					int y_end = y + _radius;
					double sum = 0;
					int count = 0;
					for(int xx = x_start ; xx <= x_end ; ++xx )
					{
						if( xx >= _image->GetWidth() || xx < 0 )
							continue;
						for(int yy = y_start; yy <= y_end; ++yy )
						{
							if( yy >= _image->GetHeight() || yy < 0 )
								continue;
							int j = xx + yy*_image->GetWidth();
							sum += _dataOutput[j];
							count++;
						}
					}
					_dataOutput[i]  =  ( _deltas[i] +  sum ) / (double) count;
				}
			}
		}
	}
	else 
	{
		for ( int x = 1 ; x < N ; ++ x )
		{
			_dataOutput[x]  =  ( _deltas[x] -  _deltas[x-1] + N*_dataOutput[x-1] ) / (double) N ;
		}
	}
}

void  Color2Gray::convertToRgb()
{
	unsigned int N  = _image->GetWidth() * _image->GetHeight() ;
	
	if( _enumType == GOOCH_GREY )
	{
		_greyOutput = new KRGBColor[N];
		for( unsigned int i = 0; i < N ; ++i )
		//for( int i = N-1 ; i >= 0 ; --i )
		{
			_greyOutput[i] = CIE_Lab(_dataOutput[i],0,0).toRGB();
			
		}
	}
	else if( _enumType == BASIC_GREY )
	{
		_basicOutput = new KRGBColor[N];
		for( unsigned int i = 0; i < N ; ++i )
		{
			_basicOutput[i] = CIE_Lab(_data[i].L,0,0).toRGB();
			
		}
	}
}

void Color2Gray::saveGreyImage()
{
	int N = _image->GetWidth() * _image->GetHeight();
	int w = _image->GetWidth();
	int h = _image->GetHeight();

	switch ( _enumType ) 
	{
	case GOOCH_GREY :
		{
			_goochImage = new KImage( w, h , _image->GetBPP() );
			_goochImage->BeginDirectAccess();

			int itr = 0;
			for( int i = _image->GetWidth() - 1 ; i >= 0 ; i-- )
			{
				for( int j = _image->GetHeight() - 1 ; j >= 0 ; j-- )
				{
					_goochImage->Put24BPPPixel( i , j  , &_greyOutput[itr] ) ;
					itr++;
				}
			}
			_goochImage->EndDirectAccess();
			break;
		}
	case BASIC_GREY :
		{
			_basicImage = new KImage( w, h , _image->GetBPP() );
			_basicImage->BeginDirectAccess();

			int itr = 0;
			for( int i = _image->GetWidth() - 1 ; i >= 0 ; i-- )
			{
				for( int j = _image->GetHeight() - 1 ; j >= 0 ; j-- )
				{
					_basicImage->Put24BPPPixel( i , j  , &_basicOutput[itr] ) ;
					itr++;
				}
			}
			_basicImage->EndDirectAccess();
			break;
		}
	default:
		{
			break;
		}
	}


}
KImage * Color2Gray::localC2G(const char * strParams, GREY_TYPE Mode )
{
	FILE *time_file = fopen("time.txt","a");
	assert(time_file);
	getParams(strParams);
	clock_t start = clock();
	convertToLab();
	
	switch (Mode) 
	{
	case GOOCH_GREY : 
		{
			clock_t start1,end1;
			double time = 0.0f;
			
			traverseCIELab();
			//here we compute the output data from the grayscale 
			computeOutput();
			
			setEnumType(this->GOOCH_GREY);
			//transform the above output data to CIELab
			convertToRgb();
			//..and finally save the image
			saveGreyImage();

			clock_t end = clock();
			time = (double)(end-start)/CLOCKS_PER_SEC;
			fprintf(time_file, "it took %.2f for GOOCH with %lu pixels and radius %d \n" , time, _image->GetWidth()*
				_image->GetHeight(),_radius);
			fclose(time_file);
			return _goochImage;
		}
	case BASIC_GREY : 
		{
			setEnumType(this->BASIC_GREY);
			//transform the above output data to CIELab
			convertToRgb();
			//..and finally save the image
			saveGreyImage();
			clock_t end = clock();
			unsigned long time = (end-start)/CLOCKS_PER_SEC;
			fclose(time_file);
			return _basicImage;
		}
	default :
		{
			break;
		}
	}
	//should never , ever, ever reach this code 
	return NULL;
}