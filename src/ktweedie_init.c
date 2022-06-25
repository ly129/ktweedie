#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <stdbool.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(td_bfgs)(int *N, double *K, double *y, double *lamreg, int *nlam, bool *sk, double *initp, double *rho, double *ftol, double *partol, double *abstol, int *maxit, bool *verbose, double *resfn, double *resgradf, double *resparam, double *resKa, int *conv);
extern void F77_NAME(td_sk)(int *N, int *p, double *x, double *y, double *lambda1, double *lambda2, double *sigmas, int *nval, double *rho, double *fn_tol, double *param_tol, double *abs_tol, double *inner_ftol, double *inner_partol, int *max_it, int *inner_maxit, bool *verbose, double *fn_final, double *gradf_final, double *param_final, double *wt_final, int *convergence);

static const R_FortranMethodDef FortranEntries[] = {
    {"td_bfgs", (DL_FUNC) &F77_NAME(td_bfgs), 18},
    {"td_sk",   (DL_FUNC) &F77_NAME(td_sk),   22},
    {NULL, NULL, 0}
};

void R_init_ktweedie(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
