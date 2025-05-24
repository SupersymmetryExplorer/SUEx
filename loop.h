#ifndef _LOOP_H_
#define _LOOP_H_


#include<stdio.h>
#include "complex.h"

double su_a(double);
void roots(double, double, double, complex *, complex *, complex *, complex *, complex *);
complex xlogx(complex);
complex fpv(int, complex, complex );
complex su_b0(double, double, double);
complex su_b1(double, double, double);
complex su_bt22(double, double, double);
complex su_bh(double, double, double);
complex su_bf(double, double, double);
complex su_bg(double, double,double);
complex f0_hdec(double, double, double);

#endif
