#include <math.h>
#include "numericx.h"

double max1(double a1,double a2,double a3,double a4)
{
	double a5,a6;
	a5 = (a1>a2)?a1:a2;
	a6 = (a3>a4)?a3:a4;
	a1 = (a5>a6)?a5:a6;
	return(a1);
	}

double min1(double a1,double a2,double a3,double a4)
{
	double a5,a6;
	a5 = (a1<a2)?a1:a2;
	a6 = (a3<a4)?a3:a4;
	a1 = (a5<a6)?a5:a6;
	return(a1);
	}
		
int sgn(double x)
{
	if(x>=0) return 1;
	return -1;
	}	

double sign(double x,double y)
{
	if(y>=0) return fabs(x);
	return -fabs(x);
	}			

/*    
int isnan(double x)
{
	return (x!=x);
	}
	
int isinf(double x) 
{ 
	return (!isnan(x) && isnan(x - x)); 
	}
*/
