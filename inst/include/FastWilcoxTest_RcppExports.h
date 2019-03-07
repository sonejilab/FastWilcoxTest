// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_FastWilcoxTest_RCPPEXPORTS_H_GEN_
#define RCPP_FastWilcoxTest_RCPPEXPORTS_H_GEN_

#include <RcppEigen.h>
#include <Rcpp.h>

namespace FastWilcoxTest {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("FastWilcoxTest", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("FastWilcoxTest", "_FastWilcoxTest_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in FastWilcoxTest");
            }
        }
    }

    inline double logFC(std::vector<double> A, std::vector<double> B) {
        typedef SEXP(*Ptr_logFC)(SEXP,SEXP);
        static Ptr_logFC p_logFC = NULL;
        if (p_logFC == NULL) {
            validateSignature("double(*logFC)(std::vector<double>,std::vector<double>)");
            p_logFC = (Ptr_logFC)R_GetCCallable("FastWilcoxTest", "_FastWilcoxTest_logFC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_logFC(Shield<SEXP>(Rcpp::wrap(A)), Shield<SEXP>(Rcpp::wrap(B)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline std::vector<int> minusOne(std::vector<int> X) {
        typedef SEXP(*Ptr_minusOne)(SEXP);
        static Ptr_minusOne p_minusOne = NULL;
        if (p_minusOne == NULL) {
            validateSignature("std::vector<int>(*minusOne)(std::vector<int>)");
            p_minusOne = (Ptr_minusOne)R_GetCCallable("FastWilcoxTest", "_FastWilcoxTest_minusOne");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_minusOne(Shield<SEXP>(Rcpp::wrap(X)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<std::vector<int> >(rcpp_result_gen);
    }

    inline std::vector<int> plusOne(std::vector<int> X) {
        typedef SEXP(*Ptr_plusOne)(SEXP);
        static Ptr_plusOne p_plusOne = NULL;
        if (p_plusOne == NULL) {
            validateSignature("std::vector<int>(*plusOne)(std::vector<int>)");
            p_plusOne = (Ptr_plusOne)R_GetCCallable("FastWilcoxTest", "_FastWilcoxTest_plusOne");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_plusOne(Shield<SEXP>(Rcpp::wrap(X)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<std::vector<int> >(rcpp_result_gen);
    }

    inline std::vector<double> fastWilcoxTest(std::vector<double> x, std::vector<double> y, int type) {
        typedef SEXP(*Ptr_fastWilcoxTest)(SEXP,SEXP,SEXP);
        static Ptr_fastWilcoxTest p_fastWilcoxTest = NULL;
        if (p_fastWilcoxTest == NULL) {
            validateSignature("std::vector<double>(*fastWilcoxTest)(std::vector<double>,std::vector<double>,int)");
            p_fastWilcoxTest = (Ptr_fastWilcoxTest)R_GetCCallable("FastWilcoxTest", "_FastWilcoxTest_fastWilcoxTest");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_fastWilcoxTest(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(y)), Shield<SEXP>(Rcpp::wrap(type)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<std::vector<double> >(rcpp_result_gen);
    }

    inline NumericMatrix StatTest(Eigen::MappedSparseMatrix<double> X, std::vector<int> interest, std::vector<int> background, double logFCcut = 1.0, double minPct = 0.1, bool onlyPos = false) {
        typedef SEXP(*Ptr_StatTest)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_StatTest p_StatTest = NULL;
        if (p_StatTest == NULL) {
            validateSignature("NumericMatrix(*StatTest)(Eigen::MappedSparseMatrix<double>,std::vector<int>,std::vector<int>,double,double,bool)");
            p_StatTest = (Ptr_StatTest)R_GetCCallable("FastWilcoxTest", "_FastWilcoxTest_StatTest");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_StatTest(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(interest)), Shield<SEXP>(Rcpp::wrap(background)), Shield<SEXP>(Rcpp::wrap(logFCcut)), Shield<SEXP>(Rcpp::wrap(minPct)), Shield<SEXP>(Rcpp::wrap(onlyPos)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericMatrix >(rcpp_result_gen);
    }

}

#endif // RCPP_FastWilcoxTest_RCPPEXPORTS_H_GEN_
