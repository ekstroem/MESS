#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP MESS_cmd(SEXP, SEXP);
extern SEXP MESS_filldown(SEXP);
extern SEXP MESS_mfastLmCpp(SEXP, SEXP, SEXP);
extern SEXP MESS_onemargintest(SEXP, SEXP);
extern SEXP MESS_qdiag(SEXP);
extern SEXP MESS_quadform(SEXP, SEXP, SEXP, SEXP);
extern SEXP MESS_repmat(SEXP, SEXP, SEXP);
extern SEXP MESS_tracemp(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"MESS_cmd",           (DL_FUNC) &MESS_cmd,           2},
  {"MESS_filldown",      (DL_FUNC) &MESS_filldown,      1},
  {"MESS_mfastLmCpp",    (DL_FUNC) &MESS_mfastLmCpp,    3},
  {"MESS_onemargintest", (DL_FUNC) &MESS_onemargintest, 2},
  {"MESS_qdiag",         (DL_FUNC) &MESS_qdiag,         1},
  {"MESS_quadform",      (DL_FUNC) &MESS_quadform,      4},
  {"MESS_repmat",        (DL_FUNC) &MESS_repmat,        3},
  {"MESS_tracemp",       (DL_FUNC) &MESS_tracemp,       2},
  {NULL, NULL, 0}
};

void R_init_MESS(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
