#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void abargehanfunC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void abarlogfunC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void abarpwfunC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gehan_s_obj(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void omegafun(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP _aftgee_log_ns_est(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _aftgee_log_s_est(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _aftgee_gehan_ns_est(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _aftgee_gehan_s_est(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _aftgee_gehan_s_wt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _aftgee_gehan_ns_wt(SEXP, SEXP, SEXP, SEXP);
extern SEXP _aftgee_getSuv(SEXP, SEXP, SEXP);
extern SEXP _aftgee_eResC(SEXP, SEXP, SEXP);
extern SEXP _aftgee_eResC2(SEXP, SEXP, SEXP);
extern SEXP _aftgee_gee(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _aftgee_est_No_Margin(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _aftgee_resampling_No_Margin(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


static const R_CMethodDef CEntries[] = {
	{"abargehanfunC", (DL_FUNC) &abargehanfunC, 12},
	{"abarlogfunC",   (DL_FUNC) &abarlogfunC,   12},
	{"abarpwfunC",    (DL_FUNC) &abarpwfunC,    12},
	{"gehan_s_obj",   (DL_FUNC) &gehan_s_obj,   12},
	{"omegafun",      (DL_FUNC) &omegafun,      10},
	{NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
	{"_aftgee_log_ns_est",            (DL_FUNC) &_aftgee_log_ns_est,           6},
	{"_aftgee_log_s_est",             (DL_FUNC) &_aftgee_log_s_est,            8},
	{"_aftgee_gehan_ns_est",          (DL_FUNC) &_aftgee_gehan_ns_est,         6},
	{"_aftgee_gehan_s_est",           (DL_FUNC) &_aftgee_gehan_s_est,          8},
	{"_aftgee_gehan_s_wt",            (DL_FUNC) &_aftgee_gehan_s_wt,           6},
	{"_aftgee_gehan_ns_wt",           (DL_FUNC) &_aftgee_gehan_ns_wt,          4},
	{"_aftgee_getSuv",                (DL_FUNC) &_aftgee_getSuv,               3},
	{"_aftgee_eResC",                 (DL_FUNC) &_aftgee_eResC,                3},
	{"_aftgee_eResC2",                (DL_FUNC) &_aftgee_eResC2,               3},
	{"_aftgee_gee",                   (DL_FUNC) &_aftgee_gee,                  8},
	{"_aftgee_est_No_Margin",         (DL_FUNC) &_aftgee_est_No_Margin,        9},
	{"_aftgee_resampling_No_Margin",  (DL_FUNC) &_aftgee_resampling_No_Margin, 10},
	{NULL, NULL, 0}
};


void R_init_aftgee(DllInfo *dll)
{
	R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
}
