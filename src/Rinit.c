#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

void bacon(double *y, double *w, double *medyR, double *madyR, int *nR, int *niterR, double *levelR, int *binnedR, int *verboseR,
           double *gibbsmu, double *gibbssig, double *gibbsp,
           double *alpha, double *beta, double *lambda, double *tau, double *gamma);

static const R_CMethodDef cMethods[] = {
  {"bacon", (DL_FUNC) &bacon, 17},
  NULL
};

void R_init_bacon(DllInfo *info) {
  R_registerRoutines(info, cMethods,  NULL, NULL, NULL);
};
