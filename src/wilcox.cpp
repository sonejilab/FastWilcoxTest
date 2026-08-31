#include <Rcpp.h>
#include <Rmath.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

using Rcpp::IntegerVector;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::S4;
using Rcpp::stop;

namespace {

inline double mean_or_nan(const NumericVector& x) {
  if (x.size() == 0) return NA_REAL;
  double sum = 0.0;
  for (double v : x) sum += v;
  return sum / static_cast<double>(x.size());
}

inline double historical_fc(const NumericVector& a, const NumericVector& b) {
  const double ma = mean_or_nan(a);
  const double mb = mean_or_nan(b);
  if (!R_FINITE(ma) || !R_FINITE(mb)) return NA_REAL;
  if (ma == 0.0) {
    if (mb == 0.0) return NA_REAL;
    return R_PosInf;
  }
  return mb / ma;
}

inline double historical_log_fc(const NumericVector& a, const NumericVector& b) {
  const double fc = historical_fc(a, b);
  return std::log(fc);
}

struct RankedResult {
  double rank_sum;
  double tie_coef;
};

RankedResult rank_sum_with_ties(const std::vector<double>& values, std::size_t n_interest) {
  const std::size_t n = values.size();
  if (n == 0 || n_interest > n) stop("invalid Wilcoxon input size");

  std::vector<std::pair<double, std::size_t>> ordered;
  ordered.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    ordered.emplace_back(values[i], i);
  }

  std::sort(
    ordered.begin(), ordered.end(),
    [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; }
  );

  std::vector<double> ranks(n, 0.0);
  long double tie_sum = 0.0L;

  std::size_t i = 0;
  while (i < n) {
    std::size_t j = i + 1;
    while (j < n && ordered[j].first == ordered[i].first) ++j;

    const double rank = (static_cast<double>(i + 1) + static_cast<double>(j)) / 2.0;
    for (std::size_t k = i; k < j; ++k) {
      ranks[ordered[k].second] = rank;
    }

    const long double t = static_cast<long double>(j - i);
    tie_sum += t * t * t - t;
    i = j;
  }

  double rank_sum = 0.0;
  for (std::size_t k = 0; k < n_interest; ++k) rank_sum += ranks[k];

  double tie_coef = 1.0;
  if (n > 1) {
    const long double nn = static_cast<long double>(n);
    const long double denom = nn * nn * nn - nn;
    if (denom > 0.0L) {
      tie_coef = 1.0 - static_cast<double>(tie_sum / denom);
    }
  }

  return {rank_sum, tie_coef};
}

double wmw_stat(double rank_sum, int n_interest, int n_total, double tie_coef, int type) {
  const int n_bg = n_total - n_interest;
  if (n_interest <= 0 || n_bg <= 0) stop("both Wilcoxon groups must be non-empty");

  const double u =
    static_cast<double>(n_interest) * n_bg +
    static_cast<double>(n_interest) * (n_interest + 1.0) * 0.5 -
    rank_sum;

  if (type == 3) return u;

  const double mu = static_cast<double>(n_interest) * n_bg * 0.5;
  const double sigma2 =
    static_cast<double>(n_interest) * n_bg * (n_total + 1.0) / 12.0 * tie_coef;

  if (!(sigma2 > 0.0) || !R_FINITE(sigma2)) return 1.0;

  double z = 0.0;
  double p = 1.0;

  if (type == 0 || type == 4) {
    z = (u + 0.5 - mu) / std::sqrt(sigma2);
    p = R::pnorm(z, 0.0, 1.0, /*lower_tail=*/false, /*log_p=*/false);
    return type == 0 ? p : -std::log10(std::max(p, std::numeric_limits<double>::min()));
  }

  if (type == 1 || type == 5) {
    z = (u - 0.5 - mu) / std::sqrt(sigma2);
    p = R::pnorm(z, 0.0, 1.0, /*lower_tail=*/true, /*log_p=*/false);
    return type == 1 ? p : std::log10(std::max(p, std::numeric_limits<double>::min()));
  }

  if (type == 2 || type == 6 || type == 7) {
    z = (u - mu - (u > mu ? 0.5 : -0.5)) / std::sqrt(sigma2);
    const double lower = R::pnorm(z, 0.0, 1.0, true, false);
    const double upper = R::pnorm(z, 0.0, 1.0, false, false);
    p = std::min(1.0, 2.0 * std::min(lower, upper));

    if (type == 2) return p;

    const double score = -std::log10(std::max(p, std::numeric_limits<double>::min()));
    if (type == 6) return score;
    return lower <= upper ? score : -score;
  }

  stop("type must be an integer from 0 to 7");
  return NA_REAL;
}

struct SparseCSC {
  int nrow;
  int ncol;
  IntegerVector p;
  IntegerVector i;
  NumericVector x;
};

SparseCSC parse_dgc_matrix(SEXP matrix_sexp) {
  if (!Rf_isS4(matrix_sexp) || !Rf_inherits(matrix_sexp, "dgCMatrix")) {
    stop("X must be a dgCMatrix");
  }

  S4 mat(matrix_sexp);
  IntegerVector dim = mat.slot("Dim");
  if (dim.size() != 2) stop("invalid dgCMatrix Dim slot");

  SparseCSC out{
    dim[0],
    dim[1],
    mat.slot("p"),
    mat.slot("i"),
    mat.slot("x")
  };

  if (out.p.size() != out.ncol + 1) stop("invalid dgCMatrix p slot");
  if (out.i.size() != out.x.size()) stop("invalid dgCMatrix i/x slots");
  return out;
}

std::vector<int> checked_zero_based_ids(const IntegerVector& ids, int nrow, const char* label) {
  std::vector<int> out;
  out.reserve(ids.size());
  for (int id : ids) {
    if (id == NA_INTEGER || id < 1 || id > nrow) {
      stop("%s row id %d is outside 1..%d", label, id, nrow);
    }
    out.push_back(id - 1);
  }
  return out;
}

} // namespace

extern "C" SEXP _FastWilcoxTest_fc(SEXP a_sexp, SEXP b_sexp) {
  BEGIN_RCPP
  NumericVector a(a_sexp), b(b_sexp);
  return Rcpp::wrap(historical_fc(a, b));
  END_RCPP
}

extern "C" SEXP _FastWilcoxTest_log_fc(SEXP a_sexp, SEXP b_sexp) {
  BEGIN_RCPP
  NumericVector a(a_sexp), b(b_sexp);
  return Rcpp::wrap(historical_log_fc(a, b));
  END_RCPP
}

extern "C" SEXP _FastWilcoxTest_cpp_wilcox_test(SEXP x_sexp, SEXP y_sexp, SEXP type_sexp) {
  BEGIN_RCPP
  NumericVector x(x_sexp), y(y_sexp);
  const int type = Rcpp::as<int>(type_sexp);

  std::vector<double> total;
  total.reserve(x.size() + y.size());
  total.insert(total.end(), x.begin(), x.end());
  total.insert(total.end(), y.begin(), y.end());

  const auto ranked = rank_sum_with_ties(total, static_cast<std::size_t>(x.size()));
  NumericVector result = NumericVector::create(
    Rcpp::_["rank.sum"] = ranked.rank_sum,
    Rcpp::_["p.value"] = wmw_stat(
      ranked.rank_sum,
      x.size(),
      total.size(),
      ranked.tie_coef,
      type
    )
  );
  return result;
  END_RCPP
}

extern "C" SEXP _FastWilcoxTest_stat_test(
    SEXP x_sexp,
    SEXP interest_sexp,
    SEXP background_sexp,
    SEXP logfc_cut_sexp,
    SEXP min_pct_sexp,
    SEXP only_pos_sexp) {
  BEGIN_RCPP

  const SparseCSC mat = parse_dgc_matrix(x_sexp);
  const IntegerVector interest_r(interest_sexp);
  const IntegerVector background_r(background_sexp);
  const double logfc_cut = Rcpp::as<double>(logfc_cut_sexp);
  const double min_pct = Rcpp::as<double>(min_pct_sexp);
  const bool only_pos = Rcpp::as<bool>(only_pos_sexp);

  const std::vector<int> interest = checked_zero_based_ids(interest_r, mat.nrow, "interest");
  const std::vector<int> background = checked_zero_based_ids(background_r, mat.nrow, "background");

  if (interest.empty()) stop("No values in interest group");
  if (background.empty()) stop("No values in background group");

  // row -> position within the compact test vector. -1 means row is unused.
  std::vector<int> interest_pos(mat.nrow, -1);
  std::vector<int> background_pos(mat.nrow, -1);

  for (std::size_t k = 0; k < interest.size(); ++k) {
    if (interest_pos[interest[k]] >= 0) stop("duplicate row id in interest");
    interest_pos[interest[k]] = static_cast<int>(k);
  }
  for (std::size_t k = 0; k < background.size(); ++k) {
    if (background_pos[background[k]] >= 0) stop("duplicate row id in background");
    if (interest_pos[background[k]] >= 0) stop("interest and background overlap");
    background_pos[background[k]] = static_cast<int>(k);
  }

  struct Candidate {
    int col;
    double logfc;
    double frac_in;
    double frac_out;
    double fc;
  };

  std::vector<Candidate> candidates;
  candidates.reserve(mat.ncol);

  for (int col = 0; col < mat.ncol; ++col) {
    double sum_in = 0.0;
    double sum_out = 0.0;
    int detected_in = 0;
    int detected_out = 0;

    for (int k = mat.p[col]; k < mat.p[col + 1]; ++k) {
      const int row = mat.i[k];
      const double value = mat.x[k];

      const int ip = interest_pos[row];
      if (ip >= 0) {
        sum_in += value;
        if (value > 0.0) ++detected_in;
        continue;
      }

      const int bp = background_pos[row];
      if (bp >= 0) {
        sum_out += value;
        if (value > 0.0) ++detected_out;
      }
    }

    const double frac_in = static_cast<double>(detected_in) / interest.size();
    const double frac_out = static_cast<double>(detected_out) / background.size();

    const double mean_in = sum_in / interest.size();
    const double mean_out = sum_out / background.size();

    double fc = NA_REAL;
    if (mean_in == 0.0) {
      if (mean_out > 0.0) fc = R_PosInf;
    } else {
      fc = mean_out / mean_in;
    }

    const double lfc = std::log(fc);
    const bool pct_pass = frac_in > min_pct || frac_out > min_pct;

    double direction_fc = fc;
    if (!only_pos && R_FINITE(direction_fc) && direction_fc < 1.0) {
      direction_fc = 1.0 / direction_fc;
    }

    // Preserve historical strict threshold semantics: >, not >=.
    const bool fc_pass =
      (R_FINITE(direction_fc) && direction_fc > std::exp(logfc_cut)) ||
      (direction_fc == R_PosInf);

    if (pct_pass && fc_pass) {
      candidates.push_back({col, lfc, frac_in, frac_out, fc});
    }
  }

  NumericMatrix result(candidates.size(), 7);
  Rcpp::colnames(result) = Rcpp::CharacterVector::create(
    "colID", "logFC", "fracExprIN", "fracExprOUT", "rank.sum", "p.value", "FC"
  );

  const std::size_t n_interest = interest.size();
  const std::size_t n_background = background.size();
  std::vector<double> total(n_interest + n_background, 0.0);

  for (std::size_t out_row = 0; out_row < candidates.size(); ++out_row) {
    const Candidate& cand = candidates[out_row];
    std::fill(total.begin(), total.end(), 0.0);

    for (int k = mat.p[cand.col]; k < mat.p[cand.col + 1]; ++k) {
      const int row = mat.i[k];
      const double value = mat.x[k];

      const int ip = interest_pos[row];
      if (ip >= 0) {
        total[static_cast<std::size_t>(ip)] = value;
        continue;
      }

      const int bp = background_pos[row];
      if (bp >= 0) {
        total[n_interest + static_cast<std::size_t>(bp)] = value;
      }
    }

    const auto ranked = rank_sum_with_ties(total, n_interest);
    const double p = wmw_stat(
      ranked.rank_sum,
      static_cast<int>(n_interest),
      static_cast<int>(total.size()),
      ranked.tie_coef,
      2
    );

    result(out_row, 0) = cand.col + 1;
    result(out_row, 1) = cand.logfc;
    result(out_row, 2) = cand.frac_in;
    result(out_row, 3) = cand.frac_out;
    result(out_row, 4) = ranked.rank_sum;
    result(out_row, 5) = p;
    result(out_row, 6) = cand.fc;
  }

  return result;
  END_RCPP
}
