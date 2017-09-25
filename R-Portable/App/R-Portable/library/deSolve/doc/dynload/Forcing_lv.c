/* compile within R with system("R CMD SHLIB Forcing_lv.c") */

#include <R.h>

static double parms[6];
static double forc[1];

/* A trick to keep up with the parameters and forcings */
#define b parms[0]
#define c parms[1]
#define d parms[2]
#define e parms[3]
#define f parms[4]
#define g parms[5]

#define import forc[0]

/* initializers*/
void parmsc(void (* odeparms)(int *, double *))
{
    int N=6;
    odeparms(&N, parms);
}

void forcc(void (* odeforcs)(int *, double *))
{
    int N=1;
    odeforcs(&N, forc);
}

/* derivative function */
void derivsc(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    if (ip[0] <2) error("nout should be at least 2");
    ydot[0] = import - b*y[0]*y[1] + g*y[2];
    ydot[1] =          c*y[0]*y[1] - d*y[2]*y[1];
    ydot[2] =          e*y[1]*y[2] - f*y[2];

    yout[0] = y[0]+y[1]+y[2];
    yout[1] = import;
}

void event(int *n, double *t, double *y) {
   y[2] = y[2]*0.5;
}
