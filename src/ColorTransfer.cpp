//===========================================================================
//===========================================================================
//===========================================================================
//==     ColorTransfer.cpp   ==   Author: Diana-Maria POPA      ==
//===========================================================================
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#include "ColorTransfer.h"
#include <time.h>
//===========================================================================
//===========================================================================
#define NUM_CHANNELS	3
#define PROJECTIONS		3
static bool _firstTime = true;

//#define DEBUG

const float firstRot[18] = {0.666667f , 0.666667f , -0.333333f ,
							0.666667f ,-0.333333f , 0.666667f ,
							-0.333333f , 0.666667f , 0.666667f};

	ColorTransfer::ColorTransfer( KImage *input , KImage *palette , unsigned int iter )
	{
		unsigned int inputOffset = input->GetWidth() * input->GetHeight();
		unsigned int paletteOffset = palette->GetWidth() * palette->GetHeight();
		_inputImage = new float[3*inputOffset];
		_palette = new float[3*paletteOffset];

		input->BeginDirectAccess();

		for ( int i = 0 ; i < input->GetHeight() ; ++i )
		{
			for ( int j = 0 ; j < input->GetWidth() ; ++j )
			{
				KRGBColor color ;
				input->Get24BPPPixel(j,i,&color);
				_inputImage[ j * input->GetHeight() + i ] = static_cast<float> (color.r/255.0);
				_inputImage[inputOffset + j* input->GetHeight() + i ] = static_cast<float> (color.g/255.0);
				_inputImage[2 *inputOffset +  j* input->GetHeight() + i ] = static_cast<float> (color.b/255.0);

			}
		}
		input->EndDirectAccess();

		palette->BeginDirectAccess();

		for ( int i = 0 ; i < palette->GetHeight() ; ++i )
		{
			for ( int j = 0 ; j < palette->GetWidth() ; ++j )
			{
				KRGBColor color ;
				palette->Get24BPPPixel(j,i,&color);
				_palette[j* palette->GetHeight() + i ] = static_cast<float> (color.r/255.0);
				_palette[paletteOffset +  j* palette->GetHeight() + i ] = static_cast<float> (color.g/255.0);
				_palette[2*paletteOffset +  j* palette->GetHeight() + i ] = static_cast<float> (color.b/255.0);
			}
		}

		palette->EndDirectAccess();

		_iterations = iter;
		_input = input;
		_paletteImage = palette;

		if( _firstTime == true )
		{
			//instantiate the pdf transfer class
			_pdfTransfer  = new PdfTransfer( _iterations, PROJECTIONS, inputOffset,
				paletteOffset, input->GetHeight() ,input->GetWidth() );
			_matrixModule = new MatrixModule();

			_firstTime = false;


		}
	}

	KImage * ColorTransfer::applyRecoloring()
	{
		FILE *fp = fopen("total_time.txt","w");
		clock_t end,start;
		double time = 0.0f;
		start = clock();
		generateRotationMatrix();
		assert(_pdfTransfer); //the one that was allocated in constructor(see above!)
		if( _pdfTransfer->NDPdfTransfer( _inputImage, _palette, _rotations ) != true )
		{
                    printf("Error in applying NDPdfTransfer for image. Exiting...\n");
                    exit(1);
		}
		if( this->saveImage() == false )
		{
                    printf("Error in  saving image for color transfer. Exiting...\n");
                    exit(1);
		}
		end = clock();
		time = (double)(end-start)/CLOCKS_PER_SEC;
		fprintf(fp, "It took %.2f for the color transfer with: \n input ->>> height = %d width = %d \n palette ->>> height = %d width = %d \n " ,
			time,_input->GetHeight(),_input->GetWidth(), _paletteImage->GetHeight(), _paletteImage->GetWidth() );
		fclose(fp);
		return this->_gradedImage;
	}

	void ColorTransfer::generateRotationMatrix()
	{
		_rotations = new float*[_iterations];
		assert(_matrixModule);


		for (unsigned int i = 0 ; i < _iterations ; ++ i )
		{
			_rotations[i]  = new float[NUM_CHANNELS * PROJECTIONS ];
			if( !i )
			{
				//_rotations[i] = firstRot;
				_matrixModule->matrix_copy(_rotations[i],firstRot,NUM_CHANNELS*PROJECTIONS);
			}
			else
			{
				float * rand = _matrixModule->matrix_generate_random(9);
				float * Q = new float[9];
				float * R = new float[9];

				for (unsigned int k = 0 ; k < 9 ; ++ k )
				{
					Q[k] = R[k]=0;
				}

				_matrixModule->matrix_mygram( rand , 3 , Q ,R );

#ifdef DEBUG
				char filename[30];
				sprintf(filename,"gram-schmidt_%d.txt",i);
				FILE *fp = fopen(filename,"w");
				for (unsigned int k = 0 ; k < 3 ; ++ k )
				{
					for (unsigned int l = 0 ; l < 3 ; ++ l  )
					{
						fprintf ( fp , " %.5f %.5f %.5f  : " , Q[k*3 + l]  ,R[k * 3 + l], rand[k * 3 + l]) ;
					}
					fprintf ( fp , "\n" ) ;
				}
				fclose(fp);
#endif
				_matrixModule->matrix_multiply( _rotations[0] , PROJECTIONS ,
					NUM_CHANNELS , Q , NUM_CHANNELS , NUM_CHANNELS , _rotations[i] ,
					PROJECTIONS , NUM_CHANNELS );
				delete [] rand;
				delete [] Q;
			}

		}

#ifdef DEBUG
		FILE *fp = fopen("rotations.txt","w");
		for (unsigned int i = 0 ; i < _iterations ; ++ i )
		{
			for (unsigned int j = 0 ; j < NUM_CHANNELS*PROJECTIONS ; ++ j  )
			{
				fprintf ( fp , "%.5f " , _rotations[i][j] ) ;
			}
			fprintf ( fp , "\n" ) ;
		}
		fclose(fp);
#endif
	}

	bool ColorTransfer::saveImage()
	{
		assert(_inputImage);
		assert(_input);
		unsigned int offset = _input->GetWidth() * _input->GetHeight();
		_gradedImage = new KImage( _input->GetWidth() , _input->GetHeight() , _input->GetBPP() );
		_gradedImage->BeginDirectAccess();
                printf ("height of input %d:%d:", _input->GetHeight(), _input->GetWidth());
		for ( int i = 0; i < _input->GetHeight(); ++i )
		{
			for ( int j = 0 ; j < _input->GetWidth() ; ++j )

			{
				BYTE R = (BYTE)adjustPixel((_inputImage[j * _input->GetHeight() + i ]*255.0f));
				BYTE G = (BYTE)adjustPixel((_inputImage[ offset + j * _input->GetHeight() + i] *255.0f));
				BYTE B = (BYTE)adjustPixel((_inputImage[ 2 * offset +  j * _input->GetHeight() + i ] * 255.0f));
				KRGBColor color(R,G,B);
				_gradedImage->Put24BPPPixel( j , i , &color ) ;
			}
		}

		_gradedImage->EndDirectAccess();
		return true;

	}
        
	float ColorTransfer::adjustPixel(float toAdjustPixel)
	{
		if( toAdjustPixel > 255.0 )
			return 255.0;
		if( toAdjustPixel < 0 )
			return 0.0;
		return toAdjustPixel;
	}