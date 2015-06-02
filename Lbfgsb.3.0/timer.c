/* timer.f -- translated by f2c (version 20090411).
*/

#include <time.h>
//#include "f2c.h"

/* Subroutine */ int timer_(doublereal *ttime)
{
    //extern /* Subroutine */ int cpu_time__(real *);
    static real temp;
	*ttime = (float)clock()/CLOCKS_PER_SEC; return 0; 
    //temp = (real) (*ttime);  cpu_time__(&temp); *ttime = (doublereal) temp; return 0;

/*     This routine computes cpu time in double precision; it makes use of the intrinsic f90 cpu_time therefore a conversion type is needed. */

} /* timer_ */

