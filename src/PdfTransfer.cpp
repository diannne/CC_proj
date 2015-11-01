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
#include "PdfTransfer.h"
#include "MatrixModule.h"
#include <math.h>
#include <algorithm>
#include <time.h>
#include <cmath>
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#pragma comment(lib, "./FreeImage/FreeImage.lib")
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#define NUM_CHANNELS 3
const int bin_count = 300;
static float _root ;

static int sBigPos;
static int sSmallPos;
static int sEqual ;
bool sEqualPos;
//#define DEBUG
PdfTransfer::PdfTransfer(unsigned int iter , unsigned int proj , unsigned int size1, 
	unsigned int size2,int _h, int _w)
{
	_iterations = iter;
	_projections = proj;
	_inputSize = size1;
	_paletteSize = size2;
	_width = _w;
	_height = _h;
	_matrixModule = new MatrixModule();
	assert(_matrixModule);
}
PdfTransfer::~PdfTransfer()
{

}
void printMatrix(float *Matrix, unsigned int lines, unsigned int cols, char *mname )
{
	char name[100];
	sprintf(name, "test_%s.txt", mname);
	FILE *fp = fopen(name,"w");
	assert(fp);
	
	for(unsigned int j = 0 ; j < cols ; j++ )
	{
		for (unsigned int i = 0 ; i < lines ; i++ )
		{
			fprintf( fp, "%.5f",matrix_access(Matrix , lines,cols, i , j )*255);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}
void printHistogr( Histogram **input, int it)
{
	char filename[20];
	sprintf(filename,"%s%d.txt","histo1",it);
	FILE *fp = fopen(filename,"w");
	assert(fp);
	
	for( unsigned int j = 0 ; j < bin_count + 1; ++j )
	{
		fprintf( fp, "%d %d\n",j,input[0][j].getCounter() );
	}
	fprintf( fp,"\n");
	fclose(fp);
	sprintf(filename,"%s%d.txt","histo2",it);
	fp = fopen(filename,"w");
	for( unsigned int j = 0 ; j < bin_count + 1; ++j )
	{
		fprintf( fp, "%d %d\n",j,input[1][j].getCounter() );
	}
	fprintf( fp,"\n");
	fclose(fp);
	sprintf(filename,"%s%d.txt","histo3",it);
	fp = fopen(filename,"w");
	for( unsigned int j = 0 ; j < bin_count + 1; ++j )
	{
		fprintf( fp, "%d %d\n",j,input[2][j].getCounter() );
	}
	fprintf( fp,"\n");
	

	fclose(fp);
}
bool PdfTransfer::NDPdfTransfer( float  *input , float *palette  , float **rotations)
{
	clock_t start , end ;
	FILE *fp = fopen("time.txt","w");
	unsigned long int var_time = 0 ;
	
	float *finalPic = new float[_projections*_inputSize];

	for ( unsigned int i = 0 ; i < _iterations ; ++ i )
	{
		start = clock();
		float * ptrRotation = rotations[i];
		assert(ptrRotation);
		
		float *inputRotation = new float[_projections*_inputSize];
		float *paletteRotation = new float[_projections*_paletteSize];
		
		_matrixModule->matrix_multiply( ptrRotation , _projections , NUM_CHANNELS , 
			input , NUM_CHANNELS , _inputSize , inputRotation , _projections, _inputSize ) ;
#ifdef DEBUG
		if( !i )
		{
			printMatrix( ptrRotation , _projections, NUM_CHANNELS ,"ptrRotation" ,_height,_width) ;
			printMatrix( input , NUM_CHANNELS, _inputSize , "input",_height,_width) ;
			printMatrix( inputRotation , _projections, _inputSize , "inputRotation",_height,_width) ; 
		}
#endif
		_matrixModule->matrix_multiply( ptrRotation , _projections , NUM_CHANNELS , 
			palette , NUM_CHANNELS , _paletteSize , paletteRotation , _projections, _paletteSize ) ;

		#ifdef DEBUG
		if( !i )
		{
			printMatrix( palette , NUM_CHANNELS, _paletteSize , "palette",_height,_width) ;
			printMatrix( paletteRotation , _projections, _paletteSize , "paletteRotation",_height,_width) ; 
		}
#endif

		Histogram ** histogramsInput = new Histogram*[_projections];
		Histogram ** histogramsPalette = new Histogram*[_projections];

		//Step1: Creation of Histograms
		CreateHistograms( inputRotation , _inputSize , 
			paletteRotation , _paletteSize , _projections , histogramsInput , histogramsPalette) ;


#ifdef DEBUG
		printHistogr(histogramsInput,i);
#endif

		std::vector<float> *inter = new  std::vector<float>[_projections];
		for (unsigned int j = 0 ; j < _projections ; ++ j )
		{
			std::vector<float> ptr = SinglePdfTransfer( histogramsInput[j] , histogramsPalette[j]);

			std::vector<float> newPtr (ptr.begin()+1,ptr.end()-1);
			float scale = bin_count / (histogramsInput[j]->GetOverallMax() - 
				histogramsInput[j]->GetOverallMin());
			
			//populate X points before interpolation 
			std::vector<float> X(newPtr.size());
			for( unsigned int k = 0 ; k < X.size() ; ++k )
			{
				X[k] = (float)k;
			}
			std::vector<float> Xq(_inputSize);
			unsigned int itr = 0; 
			for (unsigned k = j* _inputSize; k < (j+1)* _inputSize ; ++ k, ++itr )
			{
				Xq[itr] = inputRotation[k];
				Xq[itr] =  (Xq[itr] - histogramsInput[j]->GetOverallMin()) * scale;
			}


			inter[j] = linearInterpolation( X , newPtr , Xq );

			for( unsigned int k = 0 ; k < inter[j].size() ; ++k )
			{
				inter[j][k] = inter[j][k]/scale + histogramsInput[j]->GetOverallMin();
				inter[j][k] -= matrix_access( inputRotation , _projections, _inputSize, j, k );
				matrix_access( finalPic , _projections, _inputSize, j, k ) = inter[j][k];

			}
		}
		
#ifdef DEBUG
		char filename2[10];
		sprintf(filename2,"finalPic%d",i);
		printMatrix(finalPic, _projections, _inputSize, filename2,_height,_width );

#endif
		//resolve linear system of ptrRotation*x = inter
		_matrixModule->matrix_inverse(ptrRotation, _projections , NUM_CHANNELS );
		float *aux  = new float[ _projections * _inputSize ];
		_matrixModule->matrix_multiply(ptrRotation, _projections , NUM_CHANNELS,
			finalPic , _projections, _inputSize ,
			aux, _projections , _inputSize);
		_matrixModule->matrix_add(  input , aux , input , _projections , _inputSize );


#ifdef DEBUG
		//if(!i)
		//{
			char filename[10];
			sprintf(filename,"iter%d",i);
			printMatrix(input, _projections, _inputSize, filename,_height,_width);
		//}
#endif

		end = clock();
		var_time = (end - start)/(CLOCKS_PER_SEC);
		fprintf(fp, "It took %lu seconds to complete %d iteration\n",var_time, i);

		//freeing memory
		delete [] inputRotation;
		delete [] paletteRotation;

	}

	//freeing heap
	delete [] finalPic;

	//used for logging time
	fclose(fp);
	return true;
}

std::vector<float> PdfTransfer::SinglePdfTransfer( Histogram *input , Histogram *palette) 
{
	/*float *original = new float[ bin_count + 3];
	float *target = new float[ bin_count + 3];*/
	std::vector<float> original (bin_count + 1);
	std::vector<float> target (bin_count + 1);
	unsigned int globalIC = input->getGlobalCounter();
	unsigned int globalTC = palette->getGlobalCounter();
	std::vector<float> finalO (bin_count + 3);
	std::vector<float> finalT (bin_count + 3);
	std::vector<float> values (bin_count + 3);
	//to help maintain the interpolation correct
	finalO[0] = 0;
	original[0] = (float)input[0].getCounter();
	finalO[1] = original[0]/globalIC;
	finalO[bin_count+2] = (float) bin_count + 1;

	finalT[0] = 0;
	target[0] = (float)palette[0].getCounter();
	finalT[1] = target[0]/globalTC;
	finalT[bin_count+2] = (float) bin_count + 1;

	for( unsigned int i = 1 ; i <= bin_count ; ++i )
	{
		original[i] = original[i-1] + input[i].getCounter();
		finalO[i+1] = original[i]/globalIC;

		target[i] = target[i-1] +  palette[i].getCounter();
		finalT[i+1] = target[i]/globalTC;
	}

	//populate values
	values[0] = 0.0;
	for( unsigned int i = 1 ; i < bin_count + 3 ; ++i )
	{
		values[i] = (float)i-1;
	}
	return linearInterpolation( finalT, values, finalO );
}

bool IsGreater( float thisOne )
{
	sBigPos++;
	if( thisOne == _root )
	{
		sEqual = sBigPos ; 
		return true;
	}
	return thisOne > _root ;
}

bool IsSmaller( float thatOne )
{
	sSmallPos --;
	if( thatOne == _root )
	{
		sEqual = sSmallPos;
		return true;
	}
	return thatOne < _root ; 
}


std::vector<float>::iterator upper_bound (std::vector<float>::iterator first, 
	std::vector<float>::iterator last, const float& val)
{
	std::vector<float>::iterator it;
	std::iterator_traits<std::vector<float>::iterator>::difference_type count, step;
	count = std::distance(first,last);
	while ( count > 0 )
	{
		it = first; 
		step = count/2;
		std::advance (it,step);
		if ( !(val<*it) )
		{ 
			first = it;
			break;
		}
		else 
		{
			count = step;
		}
	}
	std::vector<float>::iterator from = first;
	while ( from != last - 1)
	{
		if( *(from+1) <= val )
		{
			from++;
		}
		else
		{
			break;
		}

	}
	if( *from == val)
	{
		sEqualPos = true;
	}
	return from;
}

std::vector<float> PdfTransfer::linearInterpolation( std::vector<float> X ,std::vector<float> V ,
	std::vector<float> Xq )
{
	std::vector<float>::iterator it = Xq.begin();
	std::vector<float>::iterator from = Xq.end();

	std::vector<float> result( Xq.size());
	unsigned int pos = 0 ;

	while( it != from )
	{
		_root = *it;
		sEqualPos = false;
		
		////optimized version
		std::vector<float>::iterator low,up;
		low = std::lower_bound (X.begin() , X.end(), _root);
		up = upper_bound (X.begin() , X.end(), _root);
		int posS,posB;
		posS = low - X.begin();
		posB = up - X.begin();
		if( sEqualPos == true )
		{
			result[pos++] = V[posS];
		}
		else
		{
			//get last value in case no other bigger value has been found
			//there is no risk for the other  one
			if( low == X.end() )
			{
				result[pos++] = (float) V[V.size()-1];
			}

			else 
			{
				result[pos++] = (float) V[posB] + ((float)(V[posS] - V[posB])*
					( _root - *up )/(*low - *up));
			}
		}
		it++;
	}
	return result;
}

std::vector<float> PdfTransfer::oldLinearInterpolation( std::vector<float> X ,std::vector<float> V ,
	std::vector<float> Xq )
{
	std::vector<float>::iterator it = Xq.begin();
	std::vector<float>::iterator from = Xq.end();

	std::vector<float> result( Xq.size());
	unsigned int pos = 0 ;

	while( it != from )
	{
		_root = *it;
		sBigPos = -1; 
		sEqual = -3;
		sSmallPos = X.size() ;
		std::vector<float>::iterator itBig = std::find_if ( X.begin() , X.end() , IsGreater);
		std::vector<float>::reverse_iterator itSmall = std::find_if ( X.rbegin() , X.rend() , IsSmaller);

		////optimized version
		std::vector<float>::iterator low,up;
		low = std::lower_bound (X.begin() , X.end(), _root);
		//up =  upper_bound (X.begin() , X.end(), _root);
		int posS,posB;
		posS = low-X.begin();
		if( low != X.end() && itBig!= X.end() && *itBig != *low && posS != sSmallPos )
		{
			printf("Something...\n");
			posS = low-X.begin();
			posB = up-X.begin();
		}
		if( sSmallPos < 0 )
			sSmallPos = 0;
		if( sEqual != -3 )
		{
			result[pos++] = V[sBigPos];
			//result[pos] = *low;
		}
		else
		{
			/*result[pos++] = (float) sSmallPos + ((float)(sBigPos - sSmallPos)*
			( _root - *itSmall )/(*itBig - *itSmall));*/
			if( itBig == X.end() )
			{
				result[pos++] = (float) V[V.size()-1]; //get last value
			}

			else 
			{
				result[pos++] = (float) V[sSmallPos] + ((float)(V[sBigPos] - V[sSmallPos])*
					( _root - *itSmall )/(*itBig - *itSmall));
			}
		}
		it++;


	}
	return result;
}
static inline float round(float val)
{    
	return floor(val + 0.5f);
}
void PdfTransfer::CreateHistograms( float *input ,  unsigned int iSize , 
	float *palette , unsigned int pSize , 
	unsigned int projections , 
	Histogram  **orig , Histogram  **target )
{
	static int var = 0 ;
	assert(orig);
	assert(target);

	for ( unsigned int i = 0 ; i < projections ; ++i )
	{
		float minimum ,  maximum ;
		minimum = maximum = matrix_access( input , projections, iSize, i, 0 )  ;

		//finding minumum and maximum per line for input matrix
		for ( unsigned int j = 0 ; j < iSize ; ++j )
		{
			if( matrix_access( input , projections, iSize, i, j )  < minimum )
			{
				minimum = matrix_access( input , projections , iSize, i, j )  ;
			}
			if (matrix_access( input , projections, iSize, i, j )  > maximum )
			{
				maximum = matrix_access( input , projections, iSize, i, j )  ;
			}
		}

		//finding minumum and maximum per line for palette matrix
		for (unsigned int j = 0 ; j < pSize ; ++j )
		{
			if( matrix_access( palette , projections,pSize, i, j )  < minimum )
			{
				minimum = matrix_access( palette , projections, pSize, i, j )  ;
			}
			if (matrix_access( palette , projections , pSize, i, j )  > maximum )
			{
				maximum = matrix_access( palette , projections, pSize, i, j )  ;
			}
		}


		float step = (maximum - minimum)/(float) bin_count;
		orig[i] = new Histogram[ bin_count + 1];

		for (unsigned int j = 0 ; j < iSize ; ++j )
		{
			float element = matrix_access( input , projections, iSize, i, j );
			float dblPos =  ( element - minimum ) / step;
			//int pos = static_cast<int>(floor(dblPos));
			int pos = static_cast<int>(round(dblPos));
			
			assert( pos < bin_count + 1 );
			assert (pos >= 0 );
			orig[i][pos].incrementCounter();
			orig[i][pos].SetInferiorLimit( minimum + pos*step );
			orig[i][pos].SetSuperiorLimit( minimum + (pos+1)*step );
			orig[i]->incrementGlobalCounter();
			orig[i]->SetOverallMax(maximum);
			orig[i]->SetOverallMin(minimum);

		}

		target[i] = new Histogram[ bin_count + 1];
		for (unsigned int j = 0 ; j < pSize ; ++j )
		{
			float dblPos =  ( matrix_access( palette , projections, pSize, i, j ) - minimum ) / step;
			//int pos = static_cast<int>(floor(dblPos));
			int pos = static_cast<int>(round(dblPos));

			assert( pos < bin_count + 1 );
			assert (pos >= 0 );
			target[i][pos].incrementCounter();
			target[i][pos].SetInferiorLimit( minimum + pos*step );
			target[i][pos].SetSuperiorLimit( minimum + (pos+1)*step );
			target[i]->incrementGlobalCounter();
			target[i]->SetOverallMax(maximum);
			target[i]->SetOverallMin(minimum);

		}

#ifdef DEBUG

		if( !var)
		{
			FILE *fp = NULL ;
			char name[100];
			char name2[100];
			sprintf(name, "histo_orig%d.txt",i);
			sprintf(name2, "histo_palette%d.txt",i);
			fp = fopen(name , "w");
			for (unsigned int j = 0 ; j < 301 ; ++j )
			{
				fprintf(fp, "%d\n" , orig[i][j].getCounter() );
			}
			fclose(fp);
			fp = fopen(name2 , "w");
			for (unsigned int j = 0 ; j < 301 ; ++j )
			{
				fprintf(fp, "%d\n" , target[i][j].getCounter() );
			}
			fclose(fp);

		}
#endif

	}
	var++;
}
