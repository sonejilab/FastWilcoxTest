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

    inline float correlationCoefficient(std::vector<double> X, std::vector<double> Y) {
        typedef SEXP(*Ptr_correlationCoefficient)(SEXP,SEXP);
        static Ptr_correlationCoefficient p_correlationCoefficient = NULL;
        if (p_correlationCoefficient == NULL) {
            validateSignature("float(*correlationCoefficient)(std::vector<double>,std::vector<double>)");
            p_correlationCoefficient = (Ptr_correlationCoefficient)R_GetCCallable("FastWilcoxTest", "_FastWilcoxTest_correlationCoefficient");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_correlationCoefficient(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(Y)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<float >(rcpp_result_gen);
    }

    inline std::vector<double> CorMatrixIDS(Eigen::MappedSparseMatrix<double> X, std::vector<double> CMP, std::vector<int> ids) {
        typedef SEXP(*Ptr_CorMatrixIDS)(SEXP,SEXP,SEXP);
        static Ptr_CorMatrixIDS p_CorMatrixIDS = NULL;
        if (p_CorMatrixIDS == NULL) {
            validateSignature("std::vector<double>(*CorMatrixIDS)(Eigen::MappedSparseMatrix<double>,std::vector<double>,std::vector<int>)");
            p_CorMatrixIDS = (Ptr_CorMatrixIDS)R_GetCCallable("FastWilcoxTest", "_FastWilcoxTest_CorMatrixIDS");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_CorMatrixIDS(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(CMP)), Shield<SEXP>(Rcpp::wrap(ids)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<std::vector<double> >(rcpp_result_gen);
    }

    inline std::vector<double> CorMatrix(Eigen::SparseMatrix<double> X, std::vector<double> CMP) {
        typedef SEXP(*Ptr_CorMatrix)(SEXP,SEXP);
        static Ptr_CorMatrix p_CorMatrix = NULL;
        if (p_CorMatrix == NULL) {
            validateSignature("std::vector<double>(*CorMatrix)(Eigen::SparseMatrix<double>,std::vector<double>)");
            p_CorMatrix = (Ptr_CorMatrix)R_GetCCallable("FastWilcoxTest", "_FastWilcoxTest_CorMatrix");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_CorMatrix(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(CMP)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<std::vector<double> >(rcpp_result_gen);
    }

    inline NumericMatrix LinLang(Eigen::SparseMatrix<double> X, std::vector<int> Grouping, int nGroup, double minPct = 0.1, bool display_progress = true) {
        typedef SEXP(*Ptr_LinLang)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_LinLang p_LinLang = NULL;
        if (p_LinLang == NULL) {
            validateSignature("NumericMatrix(*LinLang)(Eigen::SparseMatrix<double>,std::vector<int>,int,double,bool)");
            p_LinLang = (Ptr_LinLang)R_GetCCallable("FastWilcoxTest", "_FastWilcoxTest_LinLang");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_LinLang(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(Grouping)), Shield<SEXP>(Rcpp::wrap(nGroup)), Shield<SEXP>(Rcpp::wrap(minPct)), Shield<SEXP>(Rcpp::wrap(display_progress)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericMatrix >(rcpp_result_gen);
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

    inline std::vector<double> cppWilcoxTest(std::vector<double> x, std::vector<double> y, int type) {
        typedef SEXP(*Ptr_cppWilcoxTest)(SEXP,SEXP,SEXP);
        static Ptr_cppWilcoxTest p_cppWilcoxTest = NULL;
        if (p_cppWilcoxTest == NULL) {
            validateSignature("std::vector<double>(*cppWilcoxTest)(std::vector<double>,std::vector<double>,int)");
            p_cppWilcoxTest = (Ptr_cppWilcoxTest)R_GetCCallable("FastWilcoxTest", "_FastWilcoxTest_cppWilcoxTest");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cppWilcoxTest(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(y)), Shield<SEXP>(Rcpp::wrap(type)));
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

    inline std::vector<double> toColNums(Eigen::SparseMatrix<double> data) {
        typedef SEXP(*Ptr_toColNums)(SEXP);
        static Ptr_toColNums p_toColNums = NULL;
        if (p_toColNums == NULL) {
            validateSignature("std::vector<double>(*toColNums)(Eigen::SparseMatrix<double>)");
            p_toColNums = (Ptr_toColNums)R_GetCCallable("FastWilcoxTest", "_FastWilcoxTest_toColNums");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_toColNums(Shield<SEXP>(Rcpp::wrap(data)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<std::vector<double> >(rcpp_result_gen);
    }

}

#endif // RCPP_FastWilcoxTest_RCPPEXPORTS_H_GEN_
