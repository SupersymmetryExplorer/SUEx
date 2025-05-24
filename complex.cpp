#include <math.h>
#include <stdlib.h>

#include "complex.h"

const double pi = 3.14159265359; 		


complex mult(double x,complex y)
{
	complex ans;
	ans.re = x*y.re;
	ans.im = x*y.im;
	return(ans);
	}	

complex initiallize(double re,double im)
{
	complex ans;
	ans.re = re;
	ans.im = im;
	return ans;
	}		

double carg(complex x)
{
	double arg;
	arg = atan(x.im/x.re);

	if(x.re < 0.0 && x.im >= 0.0)
		arg = pi+arg;
	else if(x.re < 0.0 && x.im < 0.0)
		arg = -pi+arg;
				
	return arg;
	}
	
complex adddouble(complex x, double y)
{
	x.re = x.re + y;
	return x;
	}

complex add(complex x, complex y)
{
	x.re = x.re + y.re;
	x.im = x.im + y.im;
	return x;
	}

complex clog(complex z)
{
	complex ans;
	ans.re = log(sqrt((z.re)*(z.re)+(z.im)*(z.im)));
	ans.im = carg(z);
	return ans;
	}

complex cmult(complex x, complex y)
{
	complex ans;
	ans.re = x.re*y.re - x.im*y.im;
	ans.im = x.re*y.im + x.im*y.re;
	return ans;
	}	
	
double cabs(complex x)
{
	double ans;
	ans = sqrt(x.re*x.re + x.im*x.im);
	return ans;
	}		
	

	
complex csqrt(complex x)
{
	complex ans;
	double modulus,argument;
	modulus = cabs(x);
	argument = carg(x);//atan(x.im/x.re);


	ans.re = sqrt(modulus)*cos(argument/2.0);
	ans.im = sqrt(modulus)*sin(argument/2.0);

	//if(x.im <= 1.0e-15 && x.re<0 && argument < 1.0e-15)
	//	{
	//		ans.re = 0.0;
	//		ans.im = sqrt(modulus);	
	//		}
		
	return ans;
	}	

complex cdiv(double x,complex z)
{
	complex ans;
	ans.re = x*z.re/(cabs(z)*cabs(z));
	ans.im = -x*z.im/(cabs(z)*cabs(z));
	return ans;	
	}	

complex cdivision(complex x, complex y)
{
	complex ans;
	ans.re = (cabs(x)/cabs(y))*cos(carg(x)-carg(y));	
	ans.im = (cabs(x)/cabs(y))*sin(carg(x)-carg(y));	
	return ans;
	}	

complex cnegative(complex x)
{
	x.re *=-1.0;
	x.im *=-1.0;
	return x;
	}
		
double creal(complex z)
{
	return z.re;
	}

complex cpow(complex z,double x)
{
	double Rz,thetaz;
	complex ans;
	
	Rz = cabs(z);
	thetaz = carg(z);
	Rz = pow(Rz,x);
	thetaz = thetaz*x;
	ans.re = Rz*cos(thetaz);
	ans.im = Rz*sin(thetaz);
	
	return ans;	
	}
