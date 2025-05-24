#ifndef _ELECTROWEAK_H_
#define _ELECTROWEAK_H_

#include <stdio.h>
#include <math.h>


void su_pixx(double, double, double, double, double *, double *, double *, double);
void su_runningcp(double, double, double, double, double, double, double, double, double *, double *, double *);
void su_ginocr(double, double, double, double, double *);
void su_vloop2(double, double, double, double, double, double*, double*, double*);
int su_stopcr(double, double, double, double, double, double*, double*, double*);
void su_topmscr(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double*);
void su_radcino(double, double, double, double, double, double, double, double, double, double, double, double, double*, double*, double*);
double su_taumscr(double, double, double, double);
void su_sqcr(double, double, double, double*);
void su_bmsusycr(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double*);
#endif