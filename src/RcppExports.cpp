// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/FastWilcoxTest.h"
#include <RcppEigen.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// logFC
double logFC(std::vector<double> A, std::vector<double> B);
static SEXP _FastWilcoxTest_logFC_try(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type A(ASEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(logFC(A, B));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _FastWilcoxTest_logFC(SEXP ASEXP, SEXP BSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_FastWilcoxTest_logFC_try(ASEXP, BSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// fastWilcoxTest
std::vector<double> fastWilcoxTest(std::vector<double> x, std::vector<double> y, int type);
static SEXP _FastWilcoxTest_fastWilcoxTest_try(SEXP xSEXP, SEXP ySEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(fastWilcoxTest(x, y, type));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _FastWilcoxTest_fastWilcoxTest(SEXP xSEXP, SEXP ySEXP, SEXP typeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_FastWilcoxTest_fastWilcoxTest_try(xSEXP, ySEXP, typeSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// StatTest
NumericMatrix StatTest(Eigen::MappedSparseMatrix<double> X, std::vector<int> interest, std::vector<int> background, double logFCcut, double minPct, bool onlyPos);
static SEXP _FastWilcoxTest_StatTest_try(SEXP XSEXP, SEXP interestSEXP, SEXP backgroundSEXP, SEXP logFCcutSEXP, SEXP minPctSEXP, SEXP onlyPosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Eigen::MappedSparseMatrix<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type interest(interestSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type background(backgroundSEXP);
    Rcpp::traits::input_parameter< double >::type logFCcut(logFCcutSEXP);
    Rcpp::traits::input_parameter< double >::type minPct(minPctSEXP);
    Rcpp::traits::input_parameter< bool >::type onlyPos(onlyPosSEXP);
    rcpp_result_gen = Rcpp::wrap(StatTest(X, interest, background, logFCcut, minPct, onlyPos));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _FastWilcoxTest_StatTest(SEXP XSEXP, SEXP interestSEXP, SEXP backgroundSEXP, SEXP logFCcutSEXP, SEXP minPctSEXP, SEXP onlyPosSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_FastWilcoxTest_StatTest_try(XSEXP, interestSEXP, backgroundSEXP, logFCcutSEXP, minPctSEXP, onlyPosSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// ZScore
Eigen::SparseMatrix<double> ZScore(Eigen::SparseMatrix<double> data, bool display_progress);
RcppExport SEXP _FastWilcoxTest_ZScore(SEXP dataSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type data(dataSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(ZScore(data, display_progress));
    return rcpp_result_gen;
END_RCPP
}

// validate (ensure exported C++ functions exist before calling them)
static int _FastWilcoxTest_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("double(*logFC)(std::vector<double>,std::vector<double>)");
        signatures.insert("std::vector<double>(*fastWilcoxTest)(std::vector<double>,std::vector<double>,int)");
        signatures.insert("NumericMatrix(*StatTest)(Eigen::MappedSparseMatrix<double>,std::vector<int>,std::vector<int>,double,double,bool)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _FastWilcoxTest_RcppExport_registerCCallable() { 
    R_RegisterCCallable("FastWilcoxTest", "_FastWilcoxTest_logFC", (DL_FUNC)_FastWilcoxTest_logFC_try);
    R_RegisterCCallable("FastWilcoxTest", "_FastWilcoxTest_fastWilcoxTest", (DL_FUNC)_FastWilcoxTest_fastWilcoxTest_try);
    R_RegisterCCallable("FastWilcoxTest", "_FastWilcoxTest_StatTest", (DL_FUNC)_FastWilcoxTest_StatTest_try);
    R_RegisterCCallable("FastWilcoxTest", "_FastWilcoxTest_RcppExport_validate", (DL_FUNC)_FastWilcoxTest_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_FastWilcoxTest_logFC", (DL_FUNC) &_FastWilcoxTest_logFC, 2},
    {"_FastWilcoxTest_fastWilcoxTest", (DL_FUNC) &_FastWilcoxTest_fastWilcoxTest, 3},
    {"_FastWilcoxTest_StatTest", (DL_FUNC) &_FastWilcoxTest_StatTest, 6},
    {"_FastWilcoxTest_ZScore", (DL_FUNC) &_FastWilcoxTest_ZScore, 2},
    {"_FastWilcoxTest_RcppExport_registerCCallable", (DL_FUNC) &_FastWilcoxTest_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_FastWilcoxTest(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
