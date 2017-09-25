/* File dedesimple.c                                                          */
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

static double parms[6];
#define f parms[0]
#define g parms[1]
#define e parms[2]
#define m parms[3]
#define tau1 parms[4]
#define tau2 parms[5]

/* Interface to dede utility functions in package deSolve                     */

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

/* Initializer                                                                */
void initmod(void (* odeparms)(int *, double *)) {
  int N = 6;
  odeparms(&N, parms);
}

/* Derivatives                                                                */
void derivs(int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip) {

  if (ip[0] < 2) error("nout should be at least 2");

  double N = y[0];
  double P = y[1];

  int nr[2] = {0, 1};             // which lags are needed? try: (0, 0)
                                  // numbering starts from zero !
  double ytau[2] = {1.0, 1.0};    // array; initialize with default values !

  double T1 = *t - tau1;
  double T2 = *t - tau2;
  
  if (*t >= fmax(tau1, tau2)) {
    //       time, lag ID, number of returned lags, return value
    lagvalue(T1, &nr[0], 1, &ytau[0]);
    lagvalue(T2, &nr[1], 1, &ytau[1]);
  }

  ydot[0] = f * N - g * N * P;
  ydot[1] = e * g * ytau[0] * ytau[1] - m * P;

  yout[0] = ytau[0];
  yout[1] = ytau[1];
}

/* ---------------------------------------------------------------------------*/

/* Version 2: A helper function and "derivs2"                                */
double getlag(double t0, double t, double tau, double ydef, int nr) {
  double T = t - tau;
  double y = ydef;
  if ((t - tau) >= t0) lagvalue(T, &nr, 1, &y);
  return y;
}

void derivs2(int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip) {

  if (ip[0] < 2) error("nout should be at least 2");

  double N = y[0];
  double P = y[1];
  double ytau[2] = {1.0, 1.0};

  ytau[0] = getlag (0, *t, tau1, ytau[0], 0);
  ytau[1] = getlag (0, *t, tau2, ytau[1], 1);

  ydot[0] = f * N - g * N * P;
  ydot[1] = e * g * ytau[0] * ytau[1] - m * P;

  yout[0] = ytau[0];
  yout[1] = ytau[1];
}
