#ifndef _COMPLEX_H_
#define _COMPLEX_H_

#include <math.h>
#include <stdlib.h>

	
typedef struct complex
{
	double re;
	double im;
	} complex;
	


complex mult(double,complex);
complex initiallize(double,double);
complex adddouble(complex, double);
complex add(complex, complex);
complex clog(complex);
complex cmult(complex, complex);
complex csqrt(complex);
complex cdiv(double,complex);
complex cdivision(complex, complex);
complex cnegative(complex);
complex cpow(complex,double);

double creal(complex);
double carg(complex);	
double cabs(complex);

#endif