//===========================================================================
//===========================================================================
//===========================================================================
//==      MatrixModule.h   ==   Author: Diana-Maria Popa       ==
//===========================================================================
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#ifndef __MATRIXMODULE__H__
#define __MATRIXMODULE__H__
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#include <stdio.h>
#include <stdlib.h>
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#define matrix_access(X,R,C,r,c) ((X)[(r)*(C)+(c)])
class MatrixModule
{
public :
	float * matrix_generate_random(unsigned int);
	void matrix_gram_schmidt( float * , unsigned int  , unsigned int ,float * );
	float * matrix_multiply(	float *  ,unsigned int , unsigned int, 
								float * , unsigned int, unsigned int, 
								float * , unsigned int , unsigned int );
	void matrix_copy(float * , const float [] , unsigned int );

	void matrix_GJelimination(float *, unsigned int , unsigned int);
	void matrix_solvelinear(float *, unsigned int , unsigned int , 
							float *, unsigned int , unsigned int , 
							float *);
	void matrix_pivot( float * _X, unsigned int _XR, unsigned int _XC, 
	unsigned int _r, unsigned int _c);
	void matrix_swaprows(float * _X, unsigned int _XR, unsigned int _XC, 
	unsigned int _r1, unsigned int _r2);
	void matrix_inverse( float * _X, unsigned int _XR, unsigned int _XC);
	void matrix_add(float *  , float * , float * , unsigned int , unsigned int );
	void matrix_gram(float *a, unsigned int size, float *q ,float *r);
	void matrix_mygram(float *a, unsigned int size, float *q ,float *r);
	
};
#endif