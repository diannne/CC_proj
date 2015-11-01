//===========================================================================
//===========================================================================
//===========================================================================
//==     MatrixModule.cpp   ==   Author: Diana-Maria POPA      ==
//===========================================================================
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#include "stdafx.h"
#include "MatrixModule.h"
#include <math.h>
#include <complex>
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#pragma comment(lib, "./FreeImage/FreeImage.lib")
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================

#define randd() ((float) rand() / (float) RAND_MAX)



float * MatrixModule::matrix_generate_random(unsigned int size)
{
	float * rand_matrix  = new float [size];
	for ( unsigned int i = 0 ;i < size  ; ++i )
	{
		rand_matrix[i] = randd();
	}
	return rand_matrix;
}

void MatrixModule::matrix_gram_schmidt(float *A, unsigned int lines, unsigned int cols, 
	float *Q   )
{
	// Orthnormalization using the Gram-Schmidt algorithm

	// validate input
	if ( lines  == 0 || cols  == 0) 
	{
		fprintf(stderr,"Invalid arguments in Gram Schmidt\n");
		exit(1);
	}

	unsigned int i;
	unsigned int j;
	unsigned int k;

	// copy into new matrix
	memmove( Q , A , lines * cols * sizeof(float));

	unsigned int n = lines ;   // dimensionality of each vector
	float *proj_ij = new float [n];

	for (unsigned int j = 0 ; j < cols ; ++j ) 
	{
		for ( unsigned int i = 0; i < j ; ++i ) 
		{

			// compute proj(v_i, v_j)
			float vij = 0.;     // dotprod(v_i, v_j)
			float vii = 0.;     // dotprod(v_i, v_i)
			float ti;
			float tj;
			for (unsigned int k = 0; k < n; ++k ) 
			{
				ti = matrix_access(Q, lines, cols, k, i);
				tj = matrix_access(Q, lines, cols, k, j);
				vij += ti * tj;
				vii += ti * ti;
			}
			// TODO : vii should be 1.0 from normalization step below
			float g = vij / vii;

			// complete projection
			for (unsigned int k = 0 ; k < n; ++k )
			{
				proj_ij[k] = matrix_access(  Q , lines, cols, k, i) * g;
			}

			// subtract projection from v_j
			for ( unsigned int k = 0; k < n; ++k ) 
			{
				matrix_access( Q , lines, cols, k, j ) -= proj_ij[k];
			}
		}

		// normalize v_j
		float vjj = 0.;     // dotprod(v_j, v_j)
		float tj  = 0.;
		for ( unsigned int k = 0; k < n; ++k ) 
		{
			tj = matrix_access( Q , lines, cols, k, j);
			vjj += tj * tj;
		}
		// TODO : check magnitude of vjj
		float g = 1. / sqrt( vjj );
		for ( unsigned int k = 0; k < n; ++k ) 
		{
			matrix_access( Q , lines, cols , k , j) *= g;
		}
	}
	delete [] proj_ij;
}

float * MatrixModule::matrix_multiply ( float * _X, unsigned int _XR, unsigned int _XC,
	float * _Y, unsigned int _YR, unsigned int _YC,
	float * _Z, unsigned int _ZR, unsigned int _ZC )
{
	// ensure lengths are valid
	if ( _ZR != _XR || _ZC != _YC || _XC != _YR ) {
		fprintf(stderr,"invalid dimensions in multiplication \n");
		exit(1);
	}

	unsigned int r, c, i;

	for (unsigned int r=0; r<_ZR; ++r ) 
	{
		for ( unsigned int c = 0 ; c < _ZC; ++c ) 
		{
			float sum = 0.0f;
			for ( unsigned int i = 0; i < _XC; ++i ) 
			{
				sum += matrix_access(_X,_XR,_XC,r,i) *
					matrix_access(_Y,_YR,_YC,i,c);
			}
			matrix_access(_Z,_ZR,_ZC,r,c) = sum ;
		}
	}
}
void MatrixModule::matrix_copy(float * toCopy , const float * fromCopy , unsigned int size ) 
{
	for(unsigned int i = 0 ; i < size ; ++ i )
	{
		toCopy[i] = fromCopy[i];
	}
}
//matrix inverse of a square matrix
void MatrixModule::matrix_inverse( float * _X, unsigned int _XR, unsigned int _XC)
{
	// ensure lengths are valid
	if ( _XR != _XC )
	{
		fprintf(stderr, "error: matrix_inv(), invalid dimensions\n");
		exit(1);
	}

	// X:
	//  x11 x12 ... x1n
	//  x21 x22 ... x2n
	//  ...
	//  xn1 xn2 ... xnn

	// allocate temporary memory
	float *x = new float[2*_XR*_XC];
	unsigned int xr = _XR;
	unsigned int xc = _XC*2;

	// x:
	//  x11 x12 ... x1n 1   0   ... 0
	//  x21 x22 ... x2n 0   1   ... 0
	//  ...
	//  xn1 xn2 ... xnn 0   0   ... 1
	for (unsigned int r = 0; r < _XR ; ++r ) 
	{
		// copy matrix elements
		for (unsigned int c = 0; c<_XC; ++c )
		{
			matrix_access(x,xr,xc,r,c) = matrix_access(_X,_XR,_XC,r,c);
		}

		// append identity matrix
		for (unsigned int c = 0; c < _XC; ++c )
		{
			matrix_access(x,xr,xc,r,_XC+c) = (r==c) ? 1 : 0;
		}
	}

	// perform Gauss-Jordan elimination on x
	// x:
	//  1   0   ... 0   y11 y12 ... y1n
	//  0   1   ... 0   y21 y22 ... y2n
	//  ...
	//  0   0   ... 1   yn1 yn2 ... ynn

	matrix_GJelimination(x,xr,xc);
	// copy result from right half of x
	for (unsigned int r=0; r<_XR; ++r ) 
	{
		for (unsigned int c=0; c<_XC; ++c )
		{
			matrix_access(_X,_XR,_XC,r,c) = matrix_access(x,xr,xc,r,_XC+c);
		}
	}
	delete [] x;
}
// Gauss-Jordan elmination
void MatrixModule::matrix_GJelimination(float * _X, unsigned int _XR, unsigned int _XC)
{
	// choose pivot rows based on maximum element along column
	float v;
	float v_max=0.;
	unsigned int r_opt=0;

	for ( unsigned int r =  0; r < _XR; ++r ) 
	{

		// check values along this column and find the optimal row
		for (unsigned int r_hat = r ; r_hat < _XR;  r_hat++) 
		{
			v = fabs( matrix_access(_X,_XR,_XC,r_hat,r) );
			// swap rows if necessary
			if ( v > v_max || r_hat == r ) 
			{
				r_opt = r_hat;
				v_max = v;
			}
		}
		// if the maximum is zero, matrix is singular
		if ( v_max == 0.0f ) 
		{
			fprintf(stderr,"warning: matrix_gjelim(), matrix singular to machine precision\n");
		}

		// if row does not match column (e.g. maximum value does not
		// lie on the diagonal) swap the rows
		if (r != r_opt) 
		{
			matrix_swaprows(_X,_XR,_XC,r,r_opt);
		}

		// pivot on the diagonal element
		matrix_pivot(_X,_XR,_XC,r,r);
	}
	// scale by diagonal
	float g;
	for (unsigned int r = 0 ; r < _XR; ++r ) 
	{
		g = 1 / matrix_access(_X,_XR,_XC,r,r);
		for (unsigned int c=0; c < _XC; c++ )
		{
			matrix_access(_X,_XR,_XC,r,c) *= g;
		}
	}
}
// pivot on element _r, _c
void MatrixModule::matrix_pivot( float * _X, unsigned int _XR, unsigned int _XC, 
	unsigned int _r, unsigned int _c)
{
	float v = matrix_access(_X,_XR,_XC,_r,_c);
	if ( v == 0 ) 
	{
		fprintf(stderr, "warning: matrix_pivot(), pivoting on zero matrix_inv 607 \n");
		return;
	}
	unsigned int r,c;

	// pivot using back-substitution
	float g;			// multiplier
	for (unsigned int r = 0; r < _XR; ++r ) 
	{

		// skip over pivot row
		if (r == _r)
			continue;

		// compute multiplier
		g = matrix_access(_X,_XR,_XC,r,_c) / v;
		// back-substitution
		for (unsigned int c = 0; c < _XC; ++c )  
		{
			matrix_access(_X,_XR,_XC,r,c) = g*matrix_access(_X,_XR,_XC,_r,c) -
				matrix_access(_X,_XR,_XC, r,c);
		}
	}
}


void MatrixModule::matrix_swaprows(float * _X, unsigned int _XR, unsigned int _XC, 
	unsigned int _r1, unsigned int _r2)
{
	unsigned int c;
	for (unsigned int c=0; c<_XC; ++c ) 
	{
		float v_tmp = matrix_access(_X,_XR,_XC,_r1,c);
		matrix_access(_X,_XR,_XC,_r1,c) = matrix_access(_X,_XR,_XC,_r2,c);
		matrix_access(_X,_XR,_XC,_r2,c) = v_tmp;
	}
}

void MatrixModule::matrix_solvelinear(float * A, unsigned int mA, unsigned int nA, 
	float *b, unsigned int mB, unsigned int nB, 
	float *)
{

}
void MatrixModule::matrix_add(float * _X,
	float * _Y,
	float * _Z,
	unsigned int _R,
	unsigned int _C)
{

	for (unsigned int i=0; i< (_R*_C) ; ++i )
		_Z[i] = _X[i] + _Y[i];
}
void MatrixModule::matrix_gram(float *a, unsigned int size, 
	float *q , float *r)
{
	for (int k=0; k<size; k++)
	{
		//r[k][k]=0; // equivalent to sum = 0
		matrix_access(r,size,size,k,k) = 0;

		for (int i=0; i<size; i++)
		{
			//r[k][k] = r[k][k] + a[i][k] * a[i][k]; // rkk = sqr(a0k) + sqr(a1k) + sqr(a2k) 
			matrix_access(r,size,size,k,k) = matrix_access(r,size,size,k,k) + 
				matrix_access(a,size,size,i,k)*matrix_access(a,size,size,i,k);
		}
		
		//r[k][k] = sqrt(r[k][k]);  // ||a||
		matrix_access(r,size,size,k,k) = sqrt(matrix_access(r,size,size,k,k));

		//cout << endl << "R"<<k<<k<<": " << r[k][k];

		for (int i=0; i<size; i++) 
		{
			//q[i][k] = a[i][k]/r[k][k];
			//cout << " q"<<i<<k<<": "<<q[i][k] << " ";

			matrix_access(q,size,size,i,k) = matrix_access(a,size,size,i,k)/matrix_access(r,size,size,k,k);
		}

		for(int j=k+1; j<size; j++) 
		{
			//r[k][j]=0;
			matrix_access(r,size,size,k,j) = 0 ; 

			for(int i = 0; i < size; i++ )
			{
				//r[k][j] += q[i][k] * a[i][j];
				matrix_access(r,size,size,k,j) += matrix_access(q,size,size,i,k)*matrix_access(a,size,size,i,j);
			}
			//cout << endl << "r"<<k<<j<<": " <<r[k][j] <<endl;

			for ( int i = 0; i < size; i++ )
			{
				//a[i][j] = a[i][j] - r[k][j]*q[i][k];

				matrix_access(a,size,size,i,j) = matrix_access(a,size,size,i,j) - 
					matrix_access(r,size,size,k,j)*matrix_access(q,size,size,i,k);
			}

			//for ( int i=0; i<3; i++ )
				//cout << "a"<<j<<": " << a[i][j]<< " ";
		}
	}
}

void MatrixModule::matrix_mygram(float *a, unsigned int size, 
	float *q , float *r)
{
	for ( int k = 0 ; k < size; k++ )
	{
		matrix_access(r,size,size,k,k) = 0;

		for (int i=0; i<size; i++)
		{
			matrix_access(r,size,size,k,k) = matrix_access(r,size,size,k,k) + 
				matrix_access(a,size,size,i,k)*matrix_access(a,size,size,i,k);
		}
		
		matrix_access(r,size,size,k,k) = sqrt(matrix_access(r,size,size,k,k));

		for (int i=0; i<size; i++) 
		{
			matrix_access(q,size,size,i,k) = matrix_access(a,size,size,i,k)/matrix_access(r,size,size,k,k);
		}
		
		for(int j = k+1; j < size; j++) 
		{
			matrix_access(r,size,size,k,j) = 0 ; 

			for(int i = 0; i < size; i++ )
			{
				matrix_access(r,size,size,k,j) += matrix_access(q,size,size,i,k)*matrix_access(a,size,size,i,j);
			}
		}
		
		for(int j = k+1; j < size; j++) 
		{
			for ( int i = 0; i < size; i++ )
			{
				matrix_access(a,size,size,i,j) = matrix_access(a,size,size,i,j) - 
					matrix_access(r,size,size,k,j)*matrix_access(q,size,size,i,k);
			}
		
		}
	}
}