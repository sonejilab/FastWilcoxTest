// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/FastWilcoxTest.h"
#include <RcppEigen.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// correlationCoefficient
float correlationCoefficient(std::vector<double> X, std::vector<double> Y);
RcppExport SEXP _FastWilcoxTest_correlationCoefficient(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(correlationCoefficient(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// CorMatrixIDS
std::vector<double> CorMatrixIDS(Eigen::MappedSparseMatrix<double> X, std::vector<double> CMP, std::vector<int> ids);
RcppExport SEXP _FastWilcoxTest_CorMatrixIDS(SEXP XSEXP, SEXP CMPSEXP, SEXP idsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MappedSparseMatrix<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type CMP(CMPSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type ids(idsSEXP);
    rcpp_result_gen = Rcpp::wrap(CorMatrixIDS(X, CMP, ids));
    return rcpp_result_gen;
END_RCPP
}
// CorMatrix
std::vector<double> CorMatrix(Eigen::SparseMatrix<double> X, std::vector<double> CMP);
RcppExport SEXP _FastWilcoxTest_CorMatrix(SEXP XSEXP, SEXP CMPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type CMP(CMPSEXP);
    rcpp_result_gen = Rcpp::wrap(CorMatrix(X, CMP));
    return rcpp_result_gen;
END_RCPP
}
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
// minusOne
std::vector<int> minusOne(std::vector<int> X);
static SEXP _FastWilcoxTest_minusOne_try(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(minusOne(X));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _FastWilcoxTest_minusOne(SEXP XSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_FastWilcoxTest_minusOne_try(XSEXP));
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
// plusOne
std::vector<int> plusOne(std::vector<int> X);
static SEXP _FastWilcoxTest_plusOne_try(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(plusOne(X));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _FastWilcoxTest_plusOne(SEXP XSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_FastWilcoxTest_plusOne_try(XSEXP));
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
// cppWilcoxTest
std::vector<double> cppWilcoxTest(std::vector<double> x, std::vector<double> y, int type);
static SEXP _FastWilcoxTest_cppWilcoxTest_try(SEXP xSEXP, SEXP ySEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(cppWilcoxTest(x, y, type));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _FastWilcoxTest_cppWilcoxTest(SEXP xSEXP, SEXP ySEXP, SEXP typeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_FastWilcoxTest_cppWilcoxTest_try(xSEXP, ySEXP, typeSEXP));
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
// MEAN_STD
NumericMatrix MEAN_STD(Eigen::SparseMatrix<double> data);
RcppExport SEXP _FastWilcoxTest_MEAN_STD(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(MEAN_STD(data));
    return rcpp_result_gen;
END_RCPP
}
// SQRT
std::vector<double> SQRT(std::vector<double> data);
RcppExport SEXP _FastWilcoxTest_SQRT(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(SQRT(data));
    return rcpp_result_gen;
END_RCPP
}
// collapse
NumericMatrix collapse(Eigen::SparseMatrix<double> X, std::vector<int> ids, int type);
RcppExport SEXP _FastWilcoxTest_collapse(SEXP XSEXP, SEXP idsSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(collapse(X, ids, type));
    return rcpp_result_gen;
END_RCPP
}
// toColNums
std::vector<double> toColNums(Eigen::SparseMatrix<double> data);
static SEXP _FastWilcoxTest_toColNums_try(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(toColNums(data));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _FastWilcoxTest_toColNums(SEXP dataSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_FastWilcoxTest_toColNums_try(dataSEXP));
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

// validate (ensure exported C++ functions exist before calling them)
static int _FastWilcoxTest_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("double(*logFC)(std::vector<double>,std::vector<double>)");
        signatures.insert("std::vector<int>(*minusOne)(std::vector<int>)");
        signatures.insert("std::vector<int>(*plusOne)(std::vector<int>)");
        signatures.insert("std::vector<double>(*cppWilcoxTest)(std::vector<double>,std::vector<double>,int)");
        signatures.insert("NumericMatrix(*StatTest)(Eigen::MappedSparseMatrix<double>,std::vector<int>,std::vector<int>,double,double,bool)");
        signatures.insert("std::vector<double>(*toColNums)(Eigen::SparseMatrix<double>)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _FastWilcoxTest_RcppExport_registerCCallable() { 
    R_RegisterCCallable("FastWilcoxTest", "_FastWilcoxTest_logFC", (DL_FUNC)_FastWilcoxTest_logFC_try);
    R_RegisterCCallable("FastWilcoxTest", "_FastWilcoxTest_minusOne", (DL_FUNC)_FastWilcoxTest_minusOne_try);
    R_RegisterCCallable("FastWilcoxTest", "_FastWilcoxTest_plusOne", (DL_FUNC)_FastWilcoxTest_plusOne_try);
    R_RegisterCCallable("FastWilcoxTest", "_FastWilcoxTest_cppWilcoxTest", (DL_FUNC)_FastWilcoxTest_cppWilcoxTest_try);
    R_RegisterCCallable("FastWilcoxTest", "_FastWilcoxTest_StatTest", (DL_FUNC)_FastWilcoxTest_StatTest_try);
    R_RegisterCCallable("FastWilcoxTest", "_FastWilcoxTest_toColNums", (DL_FUNC)_FastWilcoxTest_toColNums_try);
    R_RegisterCCallable("FastWilcoxTest", "_FastWilcoxTest_RcppExport_validate", (DL_FUNC)_FastWilcoxTest_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_FastWilcoxTest_correlationCoefficient", (DL_FUNC) &_FastWilcoxTest_correlationCoefficient, 2},
    {"_FastWilcoxTest_CorMatrixIDS", (DL_FUNC) &_FastWilcoxTest_CorMatrixIDS, 3},
    {"_FastWilcoxTest_CorMatrix", (DL_FUNC) &_FastWilcoxTest_CorMatrix, 2},
    {"_FastWilcoxTest_logFC", (DL_FUNC) &_FastWilcoxTest_logFC, 2},
    {"_FastWilcoxTest_minusOne", (DL_FUNC) &_FastWilcoxTest_minusOne, 1},
    {"_FastWilcoxTest_plusOne", (DL_FUNC) &_FastWilcoxTest_plusOne, 1},
    {"_FastWilcoxTest_cppWilcoxTest", (DL_FUNC) &_FastWilcoxTest_cppWilcoxTest, 3},
    {"_FastWilcoxTest_StatTest", (DL_FUNC) &_FastWilcoxTest_StatTest, 6},
    {"_FastWilcoxTest_ZScore", (DL_FUNC) &_FastWilcoxTest_ZScore, 2},
    {"_FastWilcoxTest_MEAN_STD", (DL_FUNC) &_FastWilcoxTest_MEAN_STD, 1},
    {"_FastWilcoxTest_SQRT", (DL_FUNC) &_FastWilcoxTest_SQRT, 1},
    {"_FastWilcoxTest_collapse", (DL_FUNC) &_FastWilcoxTest_collapse, 3},
    {"_FastWilcoxTest_toColNums", (DL_FUNC) &_FastWilcoxTest_toColNums, 1},
    {"_FastWilcoxTest_RcppExport_registerCCallable", (DL_FUNC) &_FastWilcoxTest_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_FastWilcoxTest(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
