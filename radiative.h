#ifndef _RADIATIVE_H_ 
#define _RADIATIVE_H_

#include <stdio.h>
#include <math.h>

void su_deriv1(double, double, double*);
void su_deriv2(double, double* ,double*);
int su_odeint(double*, int, double, double, double, double, double, int);
void su_ginocr(double, double, double, double, double* );
void su_sfbpmz(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double);

double su_rkqc(double*, double*, int, double*, double, double, double*, double*, double*, int);
void su_rk4(double*, double*, int, double, double, double*, void su_derivs(double, double *, double *));


#endif
