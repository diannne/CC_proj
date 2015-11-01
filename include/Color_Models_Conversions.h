//===========================================================================
//===========================================================================
//===========================================================================
//==      Color_Models_Conversions.h   ==   Author: Diana-Maria Popa       ==
//===========================================================================
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
#ifndef __COLOR_MODELS_CONVERSIONS__H__
#define __COLOR_MODELS_CONVERSIONS__H__
//===========================================================================
//===========================================================================

#include "KImage.h"
//===========================================================================
//===========================================================================

//===========================================================================
//===========================================================================
//RGB model to CIE XYZ color model
const double one_third = 0.33333333333333333333333333333333333333333333333333333333333333;
const double lab_variable = 0.008856;
inline double adjust( double input )
{
	if(input > 255.0 )
	{
		input = 255.0;
	}
	else if (input < 0.0)
	{
		input = 0.0;
	}
	return input;
}
struct CIE_XYZ 
{
	double X , Y, Z;

	CIE_XYZ (double x, double y , double z ) : X(x) , Y(y) , Z(z) {}

	CIE_XYZ (const KRGBColor &pixel)
	{
		double R_ = pixel.r/255.0;
		double G_ = pixel.g/255.0;
		double B_ = pixel.b/255.0;

		///X = (0.412453*pixel.r + 0.35758 *pixel.g + 0.180423*pixel.b)/255.0;

		X = 0.412453*R_ + 0.35758 *G_ + 0.180423*B_;
		Y = 0.212671*R_ + 0.71516 *G_ + 0.072169*B_;
		Z = 0.019334*R_ + 0.119193*G_ + 0.950227*B_;

	}

};
//CIE XYZ model to CIE LAB color model
struct CIE_Lab
{
	double L, a, b;

	CIE_Lab() {}
	CIE_Lab(double l, double a , double b ) : L(l) , a(a) , b(b) {}

	//direct conversion from rgb to lab
	CIE_Lab(const KRGBColor &pixel)
	{
		CIE_XYZ xyz(pixel);
		
		
		*this = CIE_Lab(xyz);
		
	}
	CIE_Lab(const CIE_XYZ &xyz)
	{
		double Xn, Zn, t1, t2, t3,ft1,ft2,ft3;
		
		// no need for using Yn as divider , as Yn = 1.000f
		//ISO13655 [32] states that the tristimulus values of the reference white
		//should be those of  illuminant D50
		Xn = 0.9513f;
		Zn=		1.08886f;		
		const double parameter = 16.0/116.0;
		t1 = xyz.X/Xn;
		t2 = xyz.Y;
		t3 = xyz.Z/Zn;
		
		
		/*if( t1 >  lab_variable )
		{
			t1 = pow( t1 , one_third );
		}
		else 
		{
			t1 =  7.787*t1 + parameter ;
		}

		if( t2 >  lab_variable )
		{
			t2 = pow( t2 , one_third );
			L = (116.0 * t2 ) - 16.0;			
		}
		else 
		{
			L = 903.3 * pow(t2 , one_third);
			t2 =  7.787*t2 + parameter;
		}
		if( t3 >  lab_variable )
		{
			t3 = pow( t3 , one_third );
		}
		else 
		{
			t3 =  7.787*t3 + parameter ;
		}*/
		
		if( t2 >  lab_variable )
		{
			L = (116.0 * pow( t2 , one_third ) ) - 16.0;			
		}
		else 
		{
			L = 903.3 * t2;
		}

		if( t1 >  lab_variable )
		{
			ft1 = pow( t1 , (double) one_third );
		}
		else 
		{
			ft1 =  7.787*t1 + parameter ;
		}

		if( t2 >  lab_variable )
		{
			ft2 = pow( t2 , (double) one_third );
		}
		else 
		{
			ft2 =  7.787*t2 + parameter ;
		}

		if( t3 >  lab_variable )
		{
			ft3 = pow( t3 , (double) one_third );
		}
		else 
		{
			ft3 =  7.787*t3 + parameter ;
		}

		a = 500.0f * ( ft1 - ft2  );
		b = 200.0f * ( ft2 - ft3 );
	}

	KRGBColor toRGB()
	{
		double P = (L +16.0)/116.0;
		double Xn =  0.950455;
		//no need to initialize and use it
		//double Yn =  1.0;
		double Zn =  1.088753;
		double Y = pow(P,3);
		double X = Xn * pow((P + a/500.0),3);
		double Z = Zn * pow((P - b/200.0),3);
		double aux = 3.240479 * X - 1.537150 * Y - 0.498535 * Z;
		if( aux > 1.0 )
			aux = 1.0;
		if( aux < 0.0 )
			aux = 0.0;
		aux =  aux*255.0;
		BYTE R = (BYTE)aux;
		
		aux = -0.969256 * X + 1.875992 * Y + 0.041556 * Z;
		if( aux > 1.0 )
			aux = 1.0;
		if( aux < 0.0 )
			aux = 0.0;
		aux = aux*255.0;
		BYTE G = (BYTE)aux;
		
		aux = 0.055648 * X - 0.204043 * Y + 1.057311 * Z;
		if(aux > 1.0 )
			aux = 1.0;
		if( aux < 0.0 )
			aux = 0.0;
		aux = aux*255.0;
		BYTE B = (BYTE)aux;
		
		KRGBColor color(R,G,B);
		return color;
	}

};
#endif