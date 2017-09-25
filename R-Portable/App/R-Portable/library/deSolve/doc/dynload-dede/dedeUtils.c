/* File dedeUtils.c */
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* FORTRAN-callable interface to dede utility functions in package deSolve */

void F77_SUB(lagvalue)(double *T, int *nr, int *N, double *ytau) {
  static void(*fun)(double, int*, int, double*) = NULL;
  if (fun == NULL)
    fun =  (void(*)(double, int*, int, double*))R_GetCCallable("deSolve", "lagvalue");
  fun(*T, nr, *N, ytau);
  return;
}

void F77_SUB(lagderiv)(double *T, int *nr, int *N, double *ytau) {
  static void(*fun)(double, int*, int, double*) = NULL;
  if (fun == NULL)
    fun =  (void(*)(double, int*, int, double*))R_GetCCallable("deSolve", "lagderiv");
  fun(*T, nr, *N, ytau);
  return;
}

