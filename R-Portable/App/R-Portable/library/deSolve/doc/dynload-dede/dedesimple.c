/* File dedesimple.c */
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

static double parms[2];
#define tau parms[0]
#define k parms[1]

/* Interface to dede utility functions in package deSolve */

void lagvalue(double T, int *nr, int N, double *ytau) {
  static void(*fun)(double, int*, int, double*) = NULL;
  if (fun == NULL)
    fun =  (void(*)(double, int*, int, double*))R_GetCCallable("deSolve", "lagvalue");
  return fun(T, nr, N, ytau);
}

void lagderiv(double T, int *nr, int N, double *ytau) {
  static void(*fun)(double, int*, int, double*) = NULL;
  if (fun == NULL)
    fun =  (void(*)(double, int*, int, double*))R_GetCCallable("deSolve", "lagderiv");
  return fun(T, nr, N, ytau);
}

/* Initializer  */
void initmod(void (* odeparms)(int *, double *)) {
  int N = 2;
  odeparms(&N, parms);
}

/* Derivatives */
void derivs (int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip) {

  if (ip[0] < 1) error("nout should be at least 1");

  int nr[1] = {0};            // which lags are needed?
                              // numbering starts from zero !
  double ytau[1] = {1.0};     // array; initialize with default values !
  double T = *t - tau;
  
  if (*t > tau) {
    lagvalue(T, nr, 1, ytau);
    //Rprintf("test %g %g %g \n", T, y[0], ytau[0]);
  }

  yout[0] = ytau[0];
  ydot[0] = k * ytau[0];

}
