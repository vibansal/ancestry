/* driver1.f -- translated by f2c (version 20090411).*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "f2c.h"

#define MVAL 10
#define FACTR 1.0e8
#define PGTOL 1.0e-6 

// function call to optimization using the bounded variable Limited Memory BFGS algorithm 
double optim_LBFGS_BOUNDED(int numvariables, double *initvector, double (*function)( double x[]), void (*dfunction)( double x[], double deriv[]),double *lowerbounds, double *upperbounds,int *nbd, int iprint);


double optim_LBFGS_BOUNDED(int numvariables, double *initvector, double (*function)( double x[]), void (*dfunction)( double x[], double deriv[]),double *lowerbounds, double *upperbounds,int *nbd, int iprint)
{

	// we assume that the 'initvector' or 'X' variables have been initialized, the bounds have been set (upperbounds,lowerbounds,nbd) 
	// both function and dfunction are defined functions 

	/*     We specify the tolerances in the stopping criteria. */
	int m= MVAL; double factr = FACTR; double pgtol = PGTOL; 
	int i=0;
	/*     We specify the dimension n of the sample problem and the number */
	/*        m of limited memory corrections stored.  (n and m should not */
	/*        exceed the limits nmax and mmax respectively.) */
	int n = numvariables; // numvariables
	static char task[60]; static char csave[60]; static integer isave[44]; static logical lsave[4]; double dsave[29];

	double likelihood;
	double* gradient = calloc(sizeof(double),numvariables);
	double* wa = malloc(((2*m+4)*numvariables + 12*m*m + 12*m)*sizeof(double));
	int* iwa = malloc(3*numvariables*sizeof(int));

	new_stringcopy(task, "START",5); for (i=5; i<60; i++) task[i]=' ';

	if (dfunction == NULL) 
	{
		fprintf(stderr,"current implementation of BFGS-B requires the derivative to be provided \n"); 
		return -1;
	}
	
	likelihood = (*function)(initvector);
	(*dfunction)(initvector,gradient); 

	while (1)
	{
		setulb_(&numvariables, &m, initvector, lowerbounds, upperbounds, nbd, &likelihood,gradient, &factr, &pgtol, wa, iwa, task, &iprint, csave, lsave,isave, dsave);
		if (task[0] == 'F' && task[1] == 'G')
		{
			likelihood = (*function)(initvector);
			(*dfunction)(initvector,gradient);
			if (iprint > 1) fprintf(stdout,"likelihood = %.2f\n",likelihood);
			continue; 
		}
		else if (strncmp(task,"NEW_X",5)==0) continue;
		else break; 
	}

	if (iprint >1) fprintf(stdout,"\n");
	free(wa); free(iwa); free(gradient); 
	return -1*likelihood;
}


int main(int argc,char* argv[])
{
	/* Format strings */
	static char fmt_16[] = "(/,5x,\002Solving sample problem.\002,/,5x,\002 "
		"(f = 0.0 at the optimal solution.)\002,/)";

	/* System generated locals */

	integer i__1;
	doublereal d__1, d__2;
	/* Local variables */
	static doublereal f, g[1024];
	static integer i__;
	static doublereal l[1024];
	static integer m, n;
	static doublereal u[1024], x[1024], t1, t2, wa[43251];
	static integer nbd[1024], iwa[3072];
	static char task[60];
	static doublereal factr;
	static char csave[60];
	static doublereal dsave[29];
	static integer isave[44];
	static logical lsave[4];
	static doublereal pgtol;
	extern /* Subroutine */ int setulb_(integer *, integer *, doublereal *, 
			doublereal *, doublereal *, integer *, doublereal *, doublereal *,
			doublereal *, doublereal *, doublereal *, integer *, char *, 
			integer *, char *, logical *, integer *, doublereal *, ftnlen, 
			ftnlen);
	static integer iprint;

	/* Fortran I/O blocks */
	static cilist io___11 = { 0, 6, 0, fmt_16, 0 };


	/*     This simple driver demonstrates how to call the L-BFGS-B code to */
	/*       solve a sample problem (the extended Rosenbrock function */
	/*       subject to bounds on the variables). The dimension n of this */
	/*       problem is variable. */
	/*        nmax is the dimension of the largest problem to be solved. */
	/*        mmax is the maximum number of limited memory corrections. */
	/*     Declare the variables needed by the code. */
	/*       A description of all these variables is given at the end of */
	/*       the driver. */
	/*     Declare a few additional variables for this sample problem. */
	/*     We wish to have output at every iteration. */
	iprint = 1;
	/*     We specify the tolerances in the stopping criteria. */
	factr = 1e7;
	pgtol = 1e-5;
	/*     We specify the dimension n of the sample problem and the number */
	/*        m of limited memory corrections stored.  (n and m should not */
	/*        exceed the limits nmax and mmax respectively.) */
	n = 25;
	m = 5;
	/*     We now provide nbd which defines the bounds on the variables: */
	/*                    l   specifies the lower bounds, */
	/*                    u   specifies the upper bounds. */
	/*     First set bounds on the odd-numbered variables. */
	i__1 = n;
	for (i__ = 1; i__ <= i__1; i__ += 2) {
		nbd[i__ - 1] = 2;
		l[i__ - 1] = 1.;
		u[i__ - 1] = 100.;
		/* L10: */
	}
	/*     Next set bounds on the even-numbered variables. */
	i__1 = n;
	for (i__ = 2; i__ <= i__1; i__ += 2) {
		nbd[i__ - 1] = 2;
		l[i__ - 1] = -100.;
		u[i__ - 1] = 100.;
		/* L12: */
	}
	/*     We now define the starting point. */
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		x[i__ - 1] = 3.;
		/* L14: */
	}

	// s_wsfe(&io___11);    e_wsfe(); this needs to be properly replaced
	/*     We start the iteration by initializing task. */

	new_stringcopy(task,"START",5); //s_copy(task, "START", (ftnlen)60, (ftnlen)5);
	/*        ------- the beginning of the loop ---------- */
L111:
	/*     This is the call to the L-BFGS-B code. */
	setulb_(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task, &
			iprint, csave, lsave, isave, dsave, (ftnlen)60, (ftnlen)60);
	if (task[0]=='F' && task[1]=='G') {
		/*        the minimization routine has returned to request the */
		/*        function f and gradient g values at the current x. */
		/*        Compute function value f for the sample problem. */
		/* Computing 2nd power */
		d__1 = x[0] - 1.;
		f = d__1 * d__1 * .25;
		i__1 = n;
		for (i__ = 2; i__ <= i__1; ++i__) {
			/* Computing 2nd power */
			d__2 = x[i__ - 2];
			/* Computing 2nd power */
			d__1 = x[i__ - 1] - d__2 * d__2;
			f += d__1 * d__1;
			/* L20: */
		}
		f *= 4.;
		/*        Compute gradient g for the sample problem. */
		/* Computing 2nd power */
		d__1 = x[0];
		t1 = x[1] - d__1 * d__1;
		g[0] = (x[0] - 1.) * 2. - x[0] * 16. * t1;
		i__1 = n - 1;
		for (i__ = 2; i__ <= i__1; ++i__) {
			t2 = t1;
			/* Computing 2nd power */
			d__1 = x[i__ - 1];
			t1 = x[i__] - d__1 * d__1;
			g[i__ - 1] = t2 * 8. - x[i__ - 1] * 16. * t1;
			/* L22: */
		}
		g[n - 1] = t1 * 8.;
		/*          go back to the minimization routine. */
		goto L111;
	}

	if (task[0]=='N' && task[1]=='E' && task[2] == 'W' && task[3] == '_' && task[4] == 'X') {
		goto L111;
	}
	/*        the minimization routine has returned with a new iterate, */
	/*         and we have opted to continue the iteration. */
	/*           ---------- the end of the loop ------------- */
	/*     If task is neither FG nor NEW_X we terminate execution. */
	//s_stop("", (ftnlen)0);
	return 0;
} /* MAIN__ */

///* Main program alias */ int driver_ () { MAIN__ (); return 0; }


/*     -------------------------------------------------------------- */
/*             DESCRIPTION OF THE VARIABLES IN L-BFGS-B */
/*     -------------------------------------------------------------- */

/*     numvariables is an INTEGER variable that must be set by the user to the number of variables.  It is not altered by the routine. */

/*     m/MVAL is an INTEGER variable that must be set by the user to the number of corrections used in the limited memory matrix. */
/*       It is not altered by the routine.  Values of m < 3  are not recommended, and large values of m can result in excessive */
/*       computing time. The range  3 <= MVAL <= 20 is recommended. */

/*     x/initvector is a DOUBLE PRECISION array of length n.  On initial entry */
/*       it must be set by the user to the values of the initial */
/*       estimate of the solution vector.  Upon successful exit, it */
/*       contains the values of the variables at the best point */
/*       found (usually an approximate solution). */

/*     l is a DOUBLE PRECISION array of length n that must be set by */
/*       the user to the values of the lower bounds on the variables. If */
/*       the i-th variable has no lower bound, l(i) need not be defined. */

/*     u is a DOUBLE PRECISION array of length n that must be set by */
/*       the user to the values of the upper bounds on the variables. If */
/*       the i-th variable has no upper bound, u(i) need not be defined. */

/*     nbd is an INTEGER array of dimension n that must be set by the */
/*       user to the type of bounds imposed on the variables: */
/*       nbd(i)=0 if x(i) is unbounded, */
/*              1 if x(i) has only a lower bound, */
/*              2 if x(i) has both lower and upper bounds, */
/*              3 if x(i) has only an upper bound. */

/*     f is a DOUBLE PRECISION variable.  If the routine setulb returns */
/*       with task(1:2)= 'FG', then f must be set by the user to */
/*       contain the value of the function at the point x. */

/*     g is a DOUBLE PRECISION array of length n.  If the routine setulb */
/*       returns with taskb(1:2)= 'FG', then g must be set by the user to */
/*       contain the components of the gradient at the point x. */

/*     factr is a DOUBLE PRECISION variable that must be set by the user. */
/*       It is a tolerance in the termination test for the algorithm. */
/*       The iteration will stop when */

/*        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch */

/*       where epsmch is the machine precision which is automatically */
/*       generated by the code. Typical values for factr on a computer */
/*       with 15 digits of accuracy in double precision are: */
/*       factr=1.d+12 for low accuracy; */
/*             1.d+7  for moderate accuracy; */
/*             1.d+1  for extremely high accuracy. */
/*       The user can suppress this termination test by setting factr=0. */

/*     pgtol is a double precision variable. */
/*       On entry pgtol >= 0 is specified by the user.  The iteration */
/*         will stop when */

/*                 max{|proj g_i | i = 1, ..., n} <= pgtol */

/*         where pg_i is the ith component of the projected gradient. */
/*       The user can suppress this termination test by setting pgtol=0. */

/*     wa is a DOUBLE PRECISION  array of length */
/*       (2mmax + 5)nmax + 11mmax^2 + 8mmax used as workspace. */
/*       This array must not be altered by the user. */

/*     iwa is an INTEGER  array of length 3nmax used as */
/*       workspace. This array must not be altered by the user. */

/*     task is a CHARACTER string of length 60. */
/*       On first entry, it must be set to 'START'. */
/*       On a return with task(1:2)='FG', the user must evaluate the */
/*         function f and gradient g at the returned value of x. */
/*       On a return with task(1:5)='NEW_X', an iteration of the */
/*         algorithm has concluded, and f and g contain f(x) and g(x) */
/*         respectively.  The user can decide whether to continue or stop */
/*         the iteration. */
/*       When */
/*         task(1:4)='CONV', the termination test in L-BFGS-B has been */
/*           satisfied; */
/*         task(1:4)='ABNO', the routine has terminated abnormally */
/*           without being able to satisfy the termination conditions, */
/*           x contains the best approximation found, */
/*           f and g contain f(x) and g(x) respectively; */
/*         task(1:5)='ERROR', the routine has detected an error in the */
/*           input parameters; */
/*       On exit with task = 'CONV', 'ABNO' or 'ERROR', the variable task */
/*         contains additional information that the user can print. */
/*       This array should not be altered unless the user wants to */
/*          stop the run for some reason.  See driver2 or driver3 */
/*          for a detailed explanation on how to stop the run */
/*          by assigning task(1:4)='STOP' in the driver. */

/*     iprint is an INTEGER variable that must be set by the user. */
/*       It controls the frequency and type of output generated: */
/*        iprint<0    no output is generated; */
/*        iprint=0    print only one line at the last iteration; */
/*        0<iprint<99 print also f and |proj g| every iprint iterations; */
/*        iprint=99   print details of every iteration except n-vectors; */
/*        iprint=100  print also the changes of active set and final x; */
/*        iprint>100  print details of every iteration including x and g; */
/*       When iprint > 0, the file iterate.dat will be created to */
/*                        summarize the iteration. */

/*     csave  is a CHARACTER working array of length 60. */

/*     lsave is a LOGICAL working array of dimension 4. */
/*       On exit with task = 'NEW_X', the following information is */
/*         available: */
/*       lsave(1) = .true.  the initial x did not satisfy the bounds; */
/*       lsave(2) = .true.  the problem contains bounds; */
/*       lsave(3) = .true.  each variable has upper and lower bounds. */

/*     isave is an INTEGER working array of dimension 44. */
/*       On exit with task = 'NEW_X', it contains information that */
/*       the user may want to access: */
/*         isave(30) = the current iteration number; */
/*         isave(34) = the total number of function and gradient */
/*                         evaluations; */
/*         isave(36) = the number of function value or gradient */
/*                                  evaluations in the current iteration; */
/*         isave(38) = the number of free variables in the current */
/*                         iteration; */
/*         isave(39) = the number of active constraints at the current */
/*                         iteration; */

/*         see the subroutine setulb.f for a description of other */
/*         information contained in isave */

/*     dsave is a DOUBLE PRECISION working array of dimension 29. */
/*       On exit with task = 'NEW_X', it contains information that */
/*         the user may want to access: */
/*         dsave(2) = the value of f at the previous iteration; */
/*         dsave(5) = the machine precision epsmch generated by the code; */
/*         dsave(13) = the infinity norm of the projected gradient; */

/*         see the subroutine setulb.f for a description of other */
/*         information contained in dsave */

/*     -------------------------------------------------------------- */
/*           END OF THE DESCRIPTION OF THE VARIABLES IN L-BFGS-B */
/*     -------------------------------------------------------------- */
