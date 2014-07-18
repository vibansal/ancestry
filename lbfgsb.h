/* lbfgsb.f -- translated by f2c (version 20090411). */
#ifndef LBFGSB_H 
#define LBFGSB_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "Lbfgsb.3.0/f2c.h"

extern double BFGS_DELTA;

/* code from wikipedia to calculate machine precision in constant time */
typedef union {
	  long long i64;
	    double d64;
} dbl_64;

double machine_eps ();

// replacement for fortran string copy function
int new_stringcopy(char* s1,char* s2,int sl); 

#define MVAL 20 // increasing this from 10 to 20 has beneficial effect
#define FACTR 1.0e7 // reduce this value to increase precision, if (fold-fnew) < FACTR*epsmch*max(1.0,abs(fold)), iteration stops 
#define PGTOL 1.0e-6 // if max{|proj g_i | i = 1, ..., n} <= pgtol, iteration stops, proj(g_i) is i-th component of projected gradient 


// function call to optimization using the bounded variable Limited Memory BFGS algorithm 

double optim_LBFGS_BOUNDED(int numvariables, double *initvector, double (*function)( double x[]), void (*dfunction)( double x[], double deriv[]),double *lowerbounds, double *upperbounds,int *nbd, int iprint,int* negdelta);


/* ============================================================================= */

int setulb_(integer *n, integer *m, doublereal *x, doublereal *l, doublereal *u, integer *nbd, doublereal *f, doublereal *g, doublereal *factr, doublereal *pgtol, doublereal *wa, integer *iwa, char *task, integer *iprint, char *csave, logical *lsave, integer *isave, doublereal *dsave);

#endif
