#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern "C" SEXP _FastWilcoxTest_fc(SEXP, SEXP);
extern "C" SEXP _FastWilcoxTest_log_fc(SEXP, SEXP);
extern "C" SEXP _FastWilcoxTest_cpp_wilcox_test(SEXP, SEXP, SEXP);
extern "C" SEXP _FastWilcoxTest_stat_test(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_FastWilcoxTest_fc", (DL_FUNC) &_FastWilcoxTest_fc, 2},
  {"_FastWilcoxTest_log_fc", (DL_FUNC) &_FastWilcoxTest_log_fc, 2},
  {"_FastWilcoxTest_cpp_wilcox_test", (DL_FUNC) &_FastWilcoxTest_cpp_wilcox_test, 3},
  {"_FastWilcoxTest_stat_test", (DL_FUNC) &_FastWilcoxTest_stat_test, 6},
  {NULL, NULL, 0}
};

extern "C" void R_init_FastWilcoxTest(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
