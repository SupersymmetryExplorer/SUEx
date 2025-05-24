#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include <stdio.h>
#include <math.h>


double xitla(int, double, double);
double xiter(double, double, int, double, int, double);
void alsini(double);
double runm(double, int);
void su_gaugino(double, double, double, double, double, double, double*, double*, double*);
int su_sfermion(double, double, double, double, double, double, double, double, double, double*, double*, double*, double*, double*, double*, double*);
double falphas(double, int);
void su_delrho(double, double*, double*, double*, double, double, double, double, double*);
void su_gminus2(double, double, double, double, double, double, double, double*);
int chargino(double, double, double, double, double*, int);
void matching(int, int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double*, double*, double*, double*, double*, double*, double*, double* );
void su_bsg(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double*);
void su_finetune(double, double, double, double, double*, double*, double*, double*);

#endif