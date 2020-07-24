// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace arma;

struct adfmout {
  arma::mat tests;
  arma::mat par;
  arma::mat res;
  arma::mat b_res;
};

struct adfvout {
  arma::vec tests;
  arma::vec par;
  arma::vec res;
  arma::vec b_res;
};

arma::mat lag_matrix(const arma::vec& x, const int& p, const bool& trim = true) {
  const int n = x.n_rows;
  const int k = x.n_cols;
  arma::mat lag_x = arma::mat(n, k * p, arma::fill::zeros);
  for (int j = 0; j < p; j++) {
    lag_x(arma::span(j + 1, n - 1), arma::span(k*j, k*(j+1) - 1)) = x.rows(0, n - j - 2);
  }
  return lag_x.rows(p * trim, n - 1);
}

arma::mat lag_matrix(const arma::mat& x, const int& p, const bool& trim = true) {
  const int n = x.n_rows;
  const int k = x.n_cols;
  arma::mat lag_x = arma::mat(n, k * p, arma::fill::zeros);
  for (int j = 0; j < p; j++) {
    lag_x(arma::span(j + 1, n - 1), arma::span(k*j, k*(j+1) - 1)) = x.rows(0, n - j - 2);
  }
  return lag_x.rows(p * trim, n - 1);
}

arma::vec diff(const arma::vec& x, const bool& trim = true, const double& c = 1.0) {
  const arma::vec d_x = x - c * lag_matrix(x, 1, false);
  return d_x.tail(x.n_rows - trim);
}

arma::mat diff(const arma::mat& x, const bool& trim = true, const double& c = 1.0) {
  const arma::mat d_x = x - c * lag_matrix(x, 1, false);
  return d_x.rows(trim, x.n_rows - 1);
}

double ols_cpp(const arma::vec& y, const arma::vec& x) {
  const double b = dot(x, y) / dot(x, x);
  return(b);
}

arma::vec ols_cpp(const arma::vec& y, const arma::mat& x) {
  const arma::vec b = arma::inv_sympd(x.t() * x) * x.t() * y;
  return(b);
}

arma::mat ols_cpp(const arma::mat& y, const arma::mat& x) {
  const arma::mat b = arma::inv_sympd(x.t() * x) * x.t() * y;
  return(b);
}

arma::vec de_trend(const arma::vec& y, const int& dc = 1, const bool& QD = false) {
  const int n = y.n_rows;
  const arma::vec cQD = {7, 13.5};
  arma::vec yd = y;
  if (dc > 0) {
    arma::mat d = zeros(n, dc);
    for (int i = 0; i < dc; i++) {
      d.col(i) = arma::pow(arma::linspace(1, n, n), i);
    }
    const double c_n = QD * (1 -  cQD(dc-1) / n);
    const arma::mat d_QD = diff(d, false, c_n);
    const arma::vec y_QD = diff(y, false, c_n);
    const arma::mat b = ols_cpp(y_QD, d_QD);
    yd = y - d * b;
  }
  return yd;
}

arma::mat de_trend(const arma::mat& y, const int& dc = 1, const bool& QD = false) {
  const int n = y.n_rows;
  const arma::vec cQD = {7, 13.5};
  arma::mat yd = y;
  if (dc > 0) {
    arma::mat d = arma::mat(n, dc, arma::fill::zeros);
    for (int i = 0; i < dc; i++) {
      d.col(i) = arma::pow(arma::linspace(1, n, n), i);
    }
    const double c_n = QD * (1 -  cQD(dc-1) / n);
    const arma::mat d_QD = diff(d, false, c_n);
    const arma::mat y_QD = diff(y, false, c_n);
    const arma::mat b = ols_cpp(y_QD, d_QD);
    yd = y - d * b;
  }
  return yd;
}

adfvout adf_cpp(const arma::vec& z, const int& p, const int& dc = 1, const bool& QD = false,
                const bool& trim = true, const int& trim_ic = 0) {
  arma::mat x;
  arma::vec y_ls;
  arma::mat x_ls;
  const int n = z.n_elem;

  const arma::vec y = de_trend(z, dc, QD);
  const arma::mat y_lag = lag_matrix(y, 1, false);
  const arma::vec y_dif = y - y_lag;

  if (p > 0) {
    const arma::mat y_dif_lags = lag_matrix(y_dif, p, false);
    x = arma::join_rows(y_lag, y_dif_lags);
  } else {
    x = y_lag;
  }
  const int trim_t = (trim_ic == 0) * p + trim_ic;
  if (trim) {
    x_ls = x.rows(trim_t + 1, n - 1);
    y_ls = y_dif.rows(trim_t + 1, n - 1);
  } else {
    x_ls = x;
    y_ls = y_dif;
  }
  const arma::mat xxi = arma::inv_sympd(x_ls.t() * x_ls);
  const arma::vec b = xxi * x_ls.t() * y_ls;
  const arma::vec e = y_dif - x * b;
  const arma::vec e_b = y_dif - y_lag * b(0);
  const arma::vec e_ls = e.subvec(trim_t + 1, n - 1);
  const double s2 = arma::dot(e_ls, e_ls) / e_ls.n_elem;
  const double t_adf = b(0) / std::sqrt(s2 * xxi(0, 0));
  const double c_adf = n * b(0) / (1 - arma::sum(b) + b(0));

  adfvout adf_out;
  adf_out.tests = {t_adf, c_adf};
  adf_out.par = b;
  adf_out.res = e;
  adf_out.b_res = e_b;
  return adf_out;
}

arma::vec npve_cpp(const arma::vec& z, const double& h){
  const int n = z.n_elem;
  const arma::vec r = linspace(1./n, 1, n);
  const arma::mat x = (repelem(r, 1, n) - repelem(r.t(), n, 1)) / h;
  const arma::mat k = normpdf(x);
  const arma::mat vm = k * pow(z, 2) / sum(k, 1);
  return(vm);
}

arma::vec rescale_cpp(const arma::vec& y, const double& h = 0.1, const int& p = 0, const int& dc = 1,
                      const bool& QD = false, const bool& trim = true, const int& trim_ic = 0){
  const adfvout adf_fit = adf_cpp(y, p, dc, QD, trim, trim_ic);
  const arma::vec u = adf_fit.res;
  arma::vec ydif = diff(y, false);
  arma::vec shat = sqrt(npve_cpp(u, h));
  arma::vec yscaled = cumsum(ydif / shat);
  return(yscaled);
}

double aic_cpp(const arma::vec& e, const int& k, const double& n, const double& b0, arma::mat& ylag){
  const double s2 = dot(e, e)/n;
  const double icvalue = log(s2) + k * 2. / n;
  return icvalue;
}

double bic_cpp(const arma::vec& e, const int& k, const double& n, const double& b0, arma::mat& ylag){
  const double s2 = dot(e, e)/n;
  const double icvalue = log(s2) + k * log(n) / n;
  return icvalue;
}

double maic_cpp(const arma::vec& e, const int& k, const double& n, const double& b0, arma::mat& ylag){
  const double s2 = dot(e, e)/n;
  const double tk = dot(ylag, ylag) * pow(b0, 2) / s2;
  const double icvalue = log(s2) + (k + tk) * 2. / n;
  return icvalue;
}

double mbic_cpp(const arma::vec& e, const int& k, const double& n, const double& b0, arma::mat& ylag){
  const double s2 = dot(e, e)/n;
  const double tk = dot(ylag, ylag) * pow(b0, 2) / s2;
  const double icvalue = log(s2) + (k + tk) * log(n) / n;
  return icvalue;
}

typedef double (*icFun) (const arma::vec&, const int&, const double&, const double&, arma::mat&);

icFun ic_function(const int& ic) {
  if (ic == 1) {
    return aic_cpp;
  } else if (ic == 2) {
    return bic_cpp;
  } else if (ic == 3) {
    return maic_cpp;
  } else if (ic == 4) {
    return mbic_cpp;
  } else {
    return NULL;
  }
}

adfvout adf_selectlags_cpp(const arma::vec& y, const int& pmin, const int& pmax, icFun ic_type,
                           const int& dc = 1, const bool& QD = false, const bool& ic_scale = false,
                           const double& h_rs = 0.1, const int& p_rs = 0, const bool& trim = true){
  const int n = y.n_elem;
  arma::vec ys = y;
  if (ic_scale) {
    ys = rescale_cpp(y, h_rs, p_rs, dc, false, true, 0);
  }
  arma::vec ylag = de_trend(ys, dc, false).subvec(pmax + 1, n - 2);
  arma::vec icvalue = zeros(pmax - pmin + 1);
  adfvout adfp;
  for(int ip = pmin; ip < (pmax + 1); ip++){
    adfp = adf_cpp(ys, ip, dc, false, trim, pmax);
    icvalue(ip - pmin) = ic_type(adfp.res.tail(n - pmax - 1), ip, n - pmax - 1, adfp.par(0), ylag);
  }

  const int p_opt = icvalue.index_min() + pmin;
  const adfvout ADFp = adf_cpp(y, p_opt, dc, QD, trim, 0);
  return ADFp;
}

arma::mat adf_tests_all_units_cpp(const arma::mat& y, const int& pmin, const int& pmax, icFun ic_type,
                                  const arma::vec& dc, const arma::vec& detr, const bool& ic_scale, const double& h_rs, const arma::umat& range){

  const int dclength = dc.size();
  const int N = y.n_cols;
  adfvout adf_OLS, adf_QD;
  arma::mat OLS_p = zeros(N, dclength);
  arma::mat tests_OLS = zeros(dclength, N);
  arma::mat tests_QD = zeros(dclength, N);
  arma::mat adftests;

  if (any(detr == 1)) {
    for (int iN = 0; iN < N; iN++) {
      for (int idc = 0; idc < dclength; idc++) {
        adf_OLS = adf_selectlags_cpp(y(span(range(0, iN), range(1, iN)), iN), pmin, pmax, ic_type, dc[idc], false, ic_scale, h_rs, 0, true);
        tests_OLS(idc, iN) = adf_OLS.tests(0);
        OLS_p(iN, idc) = adf_OLS.par.n_elem - 1;
      }
    }
    adftests = tests_OLS;
    if (any(detr == 2)) {
      for (int iN = 0; iN < N; iN++) {
        for (int idc = 0; idc < dclength; idc++) {
          adf_QD = adf_cpp(y(span(range(0, iN), range(1, iN)), iN), OLS_p(iN, idc), dc[idc], true, true, 0);
          tests_QD(idc, iN) = adf_QD.tests(0);
        }
      }
      adftests = join_cols(adftests, tests_QD);
    }
  } else {
    if (any(detr == 2)) {
      for (int iN = 0; iN < N; iN++) {
        for (int idc = 0; idc < dclength; idc++) {
          adf_QD = adf_selectlags_cpp(y(span(range(0, iN), range(1, iN)), iN), pmin, pmax, ic_type, dc[idc], true, ic_scale, h_rs, 0, true);
          tests_QD(idc, iN) = adf_QD.tests(0);
        }
      }
      adftests = tests_QD;
    }
  }
  return adftests;
}

// [[Rcpp::export]]
arma::mat adf_tests_panel_cpp(const arma::mat& y, const int& pmin, const int& pmax, const int& ic,
                              const arma::vec& dc, const arma::vec& detr, const bool& ic_scale, const double& h_rs, const arma::umat& range){
  icFun ic_type = ic_function(ic);
  const arma::mat adf_out = adf_tests_all_units_cpp(y, pmin, pmax, ic_type, dc, detr, ic_scale, h_rs, range);
  return(adf_out);
}

// [[Rcpp::export]]
Rcpp::List adf_panel_bootstrap_dgp_cpp(const arma::mat& y, const int& pmin, const int& pmax, const int& ic, const int& dc, const bool& QD, const bool& trim, const bool& ic_scale, const double& h_rs, const arma::umat& range){
  int N = y.n_cols;
  arma::mat tests = zeros(2, N); // row 1 contains t.ADF; row 2 contains c.ADF
  arma::mat e = datum::nan * ones(size(y));
  arma::mat eb = datum::nan * ones(size(y));
  arma::mat par = zeros(1 + pmax - pmin, N);
  arma::vec p = zeros(N);
  adfvout adf_fit;
  icFun ic_type = ic_function(ic);

  for (int iN = 0; iN < N; iN++){
    adf_fit = adf_selectlags_cpp(y(span(range(0, iN), range(1, iN)), iN), pmin, pmax, ic_type, dc, QD, ic_scale, h_rs, 0, trim);
    tests.col(iN) = adf_fit.tests;
    p[iN] = adf_fit.par.n_elem - 1;
    par.submat(0, iN, p[iN], iN) = adf_fit.par;
    eb(span(range(0, iN), range(1, iN)), iN) = adf_fit.b_res;
    e(span(range(0, iN), range(1, iN)), iN) = adf_fit.res;
  }

  return Rcpp::List::create(
    Rcpp::Named ("tests") = tests,
    Rcpp::Named ("res") = e,
    Rcpp::Named ("b_res") = eb,
    Rcpp::Named ("par") = par,
    Rcpp::Named ("p") = p);
}

arma::vec gen_AR_cpp(const arma::vec& x, const double& ar, const double& init = 0, const bool& include_init = false){
  const int T = x.size();
  arma::vec y = zeros(T+1);
  y[0] = init;
  for (int iT = 1; iT <= T; iT++) {
    y[iT] = x[iT-1] + ar * y[iT-1];
  }
  if (!include_init) {
    y = y.tail(T);
  }
  return(y);
}

arma::vec gen_AR_cpp(const arma::vec& x, const arma::vec& ar, const arma::vec& init = 0, const bool& include_init = false){
  const int T = x.size();
  const int p = ar.size();
  const arma::vec ar_rev = reverse(ar);
  arma::vec y = zeros(T + p);
  const int init_elem = init.n_elem;
  if (init_elem == p) {
    y.subvec(0, p-1) = init;
  } else {
    y.subvec(0, p-1).fill(init[0]);
  }

  for (int iT = p; iT < T + p; iT++) {
    y[iT] = x[iT-p] + dot(ar_rev, y.subvec(iT-p, iT-1));
  }
  if (!include_init) {
    y = y.tail(T);
  }
  return(y);
}

arma::mat MBB_cpp(const arma::mat& u, const arma::mat& e, const arma::vec& z, const arma::uvec& i, const int& l, const arma::mat& s, const double& ar, const arma::mat& ar_est, const arma::rowvec& y0){
  const int T = u.n_rows;
  const int N = u.n_cols;
  const int nb = ceil(double(T) / double(l));
  const arma::uvec startb = i.subvec(0, nb - 1);
  arma::mat u_star = zeros(nb*l + 1, N);
  u_star.row(0) = y0;
  for(int i = 0; i < nb; i++){
    u_star.rows(i*l + 1, i*l + l) = u.rows(startb[i], startb[i] + l - 1);
  }
  const arma::mat y_star = cumsum(u_star);
  return y_star.tail_rows(T);
}

arma::mat BWB_cpp(const arma::mat& u, const arma::mat& e, const arma::vec& z, const arma::uvec& i, const int& l, const arma::mat& s, const double& ar, const arma::mat& ar_est, const arma::rowvec& y0){
  const int T = u.n_rows;
  const int N = u.n_cols;
  const int nb = ceil(double(T) / double(l));
  const arma::mat xi_rep = repelem(z.subvec(0, nb - 1), l, N);
  const arma::mat u_star = join_cols(y0, u % xi_rep.head_rows(T));
  const arma::mat y_star = cumsum(u_star);
  return y_star.tail_rows(T);
}

arma::mat DWB_cpp(const arma::mat& u, const arma::mat& e, const arma::vec& z, const arma::uvec& i, const int& l, const arma::mat& s, const double& ar, const arma::mat& ar_est, const arma::rowvec& y0){
  const int T = u.n_rows;
  const int N = u.n_cols;
  const arma::mat xi = s * z.subvec(0, T - 1);
  const arma::mat xi_rep = arma::repelem(xi, 1, N);
  const arma::mat u_star = join_cols(y0, u % xi_rep);
  const arma::mat y_star = cumsum(u_star);
  return y_star.tail_rows(T);
}

arma::mat AWB_cpp(const arma::mat& u, const arma::mat& e, const arma::vec& z, const arma::uvec& i, const int& l, const arma::mat& s, const double& ar, const arma::mat& ar_est, const arma::rowvec& y0){
  const int T = u.n_rows;
  const int N = u.n_cols;
  const arma::vec zi = z.subvec(1, T - 1) * sqrt(1 - pow(ar, 2));
  const arma::vec xi = gen_AR_cpp(zi, ar, as_scalar(z(0)), true);
  const arma::mat xi_rep = arma::repelem(xi, 1, N);
  const arma::mat u_star = join_cols(y0, u % xi_rep);
  const arma::mat y_star = cumsum(u_star);
  return y_star.tail_rows(T);
}

arma::mat SB_cpp(const arma::mat& u, const arma::mat& e, const arma::vec& z, const arma::uvec& i, const int& l, const arma::mat& s, const double& ar, const arma::mat& ar_est, const arma::rowvec& y0){
  const int T = e.n_rows;
  const int N = e.n_cols;
  const arma::uvec index = i.subvec(0, T - 1);
  const arma::mat e_star = e.rows(index);
  arma::mat u_star = zeros(T, N);
  arma::vec init = zeros(ar_est.n_rows);
  for (int iN = 0; iN < N; iN++){
    u_star.col(iN) = gen_AR_cpp(e_star.col(iN), ar_est.col(iN), init, false);
  }
  u_star = join_cols(y0, u_star);
  const arma::mat y_star = cumsum(u_star);
  return y_star.tail_rows(T);
}

arma::mat SWB_cpp(const arma::mat& u, const arma::mat& e, const arma::vec& z, const arma::uvec& i, const int& l, const arma::mat& s, const double& ar, const arma::mat& ar_est, const arma::rowvec& y0){
  const int T = e.n_rows;
  const int N = e.n_cols;
  const arma::mat e_star = repelem(z, 1, N) % e;
  arma::mat u_star = zeros(T, N);
  arma::vec init = zeros(ar_est.n_rows);
  for (int iN = 0; iN < N; iN++){
    u_star.col(iN) = gen_AR_cpp(e_star.col(iN), ar_est.col(iN), init, false);
  }
  u_star = join_cols(y0, u_star);
  const arma::mat y_star = cumsum(u_star);
  return y_star.tail_rows(T);
}

typedef arma::mat (*bFun) (const arma::mat&, const arma::mat&, const arma::vec&, const arma::uvec&, const int&, const arma::mat&, const double&, const arma::mat&, const arma::rowvec&);

bFun boot_func(const int& boot) {
  if (boot == 1) {
    return MBB_cpp;
  } else if (boot == 2) {
    return BWB_cpp;
  } else if (boot == 3) {
    return DWB_cpp;
  } else if (boot == 4) {
    return AWB_cpp;
  } else if (boot == 5) {
    return SB_cpp;
  } else if (boot == 6) {
    return SWB_cpp;
  } else {
    return NULL;
  }
}

arma::mat bootstrap_tests_cpp(const arma::mat& u, const arma::mat& e, bFun boot_f, const arma::vec& z, const arma::uvec& i, const int& l, const arma::mat& s, const double& ar, const arma::mat& ar_est, const arma::mat& y0,
                              const int& pmin, const int& pmax, icFun ic_type, const arma::vec& dc, const arma::vec& detr, const bool& ic_scale, const double& h_rs, const arma::umat& range) {

  const arma::mat y_star = boot_f(u, e, z, i, l, s, ar, ar_est, y0);
  const arma::mat adf_btests = adf_tests_all_units_cpp(y_star, pmin, pmax, ic_type, dc, detr, ic_scale, h_rs, range);
  return adf_btests;
}

// [[Rcpp::export]]
arma::cube bootstrap_cpp(const int& B, const arma::mat& u, const arma::mat& e, const int& boot, const int& l, const arma::mat& s,
                         const double& ar, const arma::mat& ar_est, const arma::mat& y0, const int& pmin, const int& pmax, const int& ic, const arma::vec& dc,
                         const arma::vec& detr, const bool& ic_scale, const double& h_rs, const arma::umat& range, const bool& joint = true,
                         const bool& do_parallel = false, const int& nc = 1, const bool& show_progress = false){

  const int detrlength = detr.size();
  const int dclength = dc.size();
  const int N = u.n_cols;
  const int T = u.n_rows;
  int ub;

  const bFun boot_f = boot_func(boot);
  const icFun ic_type = ic_function(ic);
  arma::mat u0 = u, e0 = e;
  u0.replace(datum::nan, 0);
  e0.replace(datum::nan, 0);
  const arma::vec z_vec = Rcpp::rnorm(T * B, 0, 1);
  const arma::mat z = reshape(z_vec, T, B);

  if (boot == 6) {
    ub = 1;
  } else {
    ub = l;
  }

  const arma::uvec i_vec = Rcpp::RcppArmadillo::sample(linspace<arma::uvec>(0, T - ub, T - ub + 1), T * B, true);
  const arma::umat i = reshape(i_vec, T, B);
  arma::cube output = zeros(B, dclength * detrlength, N);

  #ifdef _OPENMP
    omp_set_num_threads(nc);
  #endif
  Progress prog(B, show_progress);
  bool break_flag = false;
  if (do_parallel) {
    if (joint) {
      #pragma omp parallel for schedule(static)
      for (int iB = 0; iB < B; iB++) {
        if ( ! Progress::check_abort() ) {
          output.subcube(iB, 0, 0, iB, dclength * detrlength - 1, N - 1) = bootstrap_tests_cpp(u0, e0, boot_f, z.col(iB), i.col(iB), l, s, ar, ar_est, y0, pmin, pmax, ic_type, dc, detr, ic_scale, h_rs, range);
          prog.increment();
        }
      }
    } else {
      #pragma omp parallel for schedule(static)
      for (int iB = 0; iB < B; iB++) {
        if ( ! Progress::check_abort() ) {
          for (int iN = 0; iN < N; iN++) {
            output.subcube(iB, 0, iN, iB, dclength * detrlength - 1, iN) = bootstrap_tests_cpp(u0.col(iN), e0.col(iN), boot_f, z.col(iB), i.col(iB), l, s, ar, ar_est.col(iN), y0.col(iN), pmin, pmax, ic_type, dc, detr, ic_scale, h_rs, range.col(iN));
          }
          prog.increment();
        }
      }
    }
  } else {
    if (joint) {
      for (int iB = 0; iB < B; iB++) {
        if ( ! Progress::check_abort() ) {
          output.subcube(iB, 0, 0, iB, dclength * detrlength - 1, N - 1) = bootstrap_tests_cpp(u0, e0, boot_f, z.col(iB), i.col(iB), l, s, ar, ar_est, y0, pmin, pmax, ic_type, dc, detr, ic_scale, h_rs, range);
          prog.increment();
        } else {
          break_flag = true;
        }
      }
    } else {
      for (int iB = 0; iB < B; iB++) {
        if ( ! Progress::check_abort() ) {
          for (int iN = 0; iN < N; iN++) {
            output.subcube(iB, 0, iN, iB, dclength * detrlength - 1, iN) = bootstrap_tests_cpp(u0.col(iN), e0.col(iN), boot_f, z.col(iB), i.col(iB), l, s, ar, ar_est.col(iN), y0.col(iN), pmin, pmax, ic_type, dc, detr, ic_scale, h_rs, range.col(iN));
          }
          prog.increment();
        }
      }
    }
  }
  if (break_flag) {
    Rcpp::stop("Bootstrap loop aborted.");
  }
  return(output);
}

arma::mat Quantile(const arma::mat& x, const arma::vec& prob, const bool& interp = false) {
  const arma::mat x_sort = sort(x, "ascend", 0);
  const arma::vec index = prob * x.n_rows - 1;
  arma::mat q = zeros(prob.n_elem, x.n_cols);
  const arma::uvec ceil_index = conv_to<uvec>::from(ceil(index));
  if (!interp) {
    q = x_sort.rows(ceil_index);
  } else {
    const arma::vec g = ceil_index - index;
    for (unsigned int i = 0; i < prob.n_elem; i++) {
      q.row(i) = g(i) * x_sort.row(ceil_index(i) - 1) + (1. - g(i)) * x_sort.row(ceil_index(i));
    }
  }
  return(q);
}

arma::rowvec Quantile(const arma::mat& x, const double& prob, const bool& interp = false) {
  const arma::mat x_sort = sort(x, "ascend", 0);
  const double index = prob * x.n_rows - 1;
  const int ceil_index = int(ceil(index));
  arma::rowvec q;
  if (!interp) {
    q = x_sort.row(ceil_index);
  } else {
    double g = ceil_index - index;
    q = g * x_sort.row(ceil_index - 1) + (1. - g) * x_sort.row(ceil_index);
  }
  return(q);
}

arma::vec Quantile(const arma::vec& x, const arma::vec& prob, const bool& interp = false) {
  const arma::vec x_sort = sort(x);
  const arma::vec index = prob * x.n_rows - 1;
  const arma::uvec ceil_index = conv_to<uvec>::from(ceil(index));
  arma::vec q;
  if (!interp) {
    q = x_sort.elem(ceil_index);
  } else {
    const arma::vec g = ceil_index - index;
    q = g % x_sort.elem(ceil_index - 1) + (1. - g) % x_sort.elem(ceil_index);
  }
  return(q);
}

double Quantile(const arma::vec& x, const double& prob, const bool& interp = false) {
  const arma::vec x_sort = sort(x);
  const double index = prob * x.n_rows - 1;
  const int ceil_index = int(ceil(index));
  double q;
  if (!interp) {
    q = x_sort(ceil_index);
  } else {
    const double g = ceil_index - index;
    q = g * x_sort(ceil_index - 1) + (1. - g) * x_sort(ceil_index);
  }
  return(q);
}

// [[Rcpp::export]]
arma::mat scaling_factors_cpp(const arma::cube& u, const double& prob){
  const int B = u.n_rows;
  const int D = u.n_cols;
  const int N = u.n_slices;
  arma::mat sc_f = zeros(D, N), oneD;

  for (int iD = 0; iD < D; iD++){
    oneD = u.subcube(0, iD, 0, B-1, iD, N-1);
    sc_f.row(iD) = Quantile(oneD, prob);
  }
  return(sc_f);
}

// [[Rcpp::export]]
arma::mat union_tests_cpp(const arma::cube& t, arma::mat& s){
  int B = t.n_rows;
  int D = t.n_cols;
  int N = t.n_slices;
  arma::mat un_tests = zeros(B, N), test;
  for (int iB = 0; iB < B; iB++) {
    test = t.subcube(iB, 0, 0, iB, D-1, N-1 );
    test = -test / s;
    un_tests.row(iB) = min(test, 0);
  }
  return(un_tests);
}

// [[Rcpp::export]]
arma::vec union_test_cpp(const arma::mat& t, arma::vec& s){
  int B = t.n_rows;
  // arma::mat s_rep = repelem(s.t(), B, 1)
  arma::vec un_tests = zeros(B);
  arma::rowvec test;
  for (int iB = 0; iB < B; iB++) {
    test = t.row(iB);
    test = -test / s.t();
    un_tests(iB) = min(test);
  }
  return(un_tests);
}

arma::rowvec BSQT_step_cpp(const int& p0, const int& p1, const arma::mat& test_i, const arma::uvec& ranks, const arma::mat& t_star){
  const int N = ranks.n_elem;
  const arma::uvec ind = ranks.subvec(p0, N - 1);
  const arma::mat t_star_sub = t_star.cols(sort(ind));

  const int index_j = ranks(p1-1);
  const double test_j = test_i(0, index_j);

  const arma::mat sort_tstar = arma::sort(t_star_sub, "ascend", 1);
  const arma::vec order_tstar = sort_tstar.col(p1 - p0 - 1);
  const double pval_j = double(sum(order_tstar < test_j)) / t_star_sub.n_rows;

  const arma::rowvec results = {double(p0), double(p1), double(index_j + 1), test_j, pval_j};
  return(results);
}

// [[Rcpp::export]]
Rcpp:: List BSQT_cpp(const arma::vec& pvec, const arma::mat& test_i, const arma::mat& t_star, const double& level){

  const int N = test_i.n_elem;
  const arma::uvec ranks = sort_index(test_i);
  const int K = pvec.n_elem - 1;
  int p_hat = N;
  arma::mat step_stats = zeros(K, 5);
  for (int i = 0; i < K; i++) {
    step_stats.row(i) = BSQT_step_cpp(pvec(i), pvec(i + 1), test_i, ranks, t_star);
    if (step_stats(i, 4) > level) {
      p_hat = step_stats(i, 0);
      step_stats = step_stats.head_rows(i + 1);
      break;
    }
  }

  arma::uvec rej_H0(N, arma::fill::zeros);
  if (p_hat > 0) {
    rej_H0.elem(ranks.head(p_hat)).ones();
  }

  return Rcpp::List::create(
    Rcpp::Named("rej_H0") = rej_H0,
    Rcpp::Named("BSQT_steps") =  step_stats,
    Rcpp::Named("ranks") = ranks
  );
}

// [[Rcpp::export]]
arma::mat iADF_cpp(const arma::mat& test_i, const arma::mat& t_star, const double& level){
  const int N = test_i.n_cols;
  arma::vec p_val_i = zeros(N);
  for (int i = 0; i < N; i++) {
    p_val_i(i) = double(sum(t_star.col(i) < test_i(i))) / t_star.n_rows;
  }
  const arma::mat Ind_Tests = join_rows(test_i.t(), p_val_i);
  return Ind_Tests;
}

// [[Rcpp::export]]
Rcpp:: List FDR_cpp(const arma::mat& test_i, const arma::mat& t_star, const double& level){
  const int N = test_i.n_elem;
  const int B = t_star.n_rows;
  const arma::uvec ranks = sort_index(test_i);

  arma::mat t_star_sub, sorted_t_star, noreject, noreject_jplus, cv_star;
  arma::mat cv_fdr = zeros(1, N);
  arma::vec noreject_prod, FDR_est;
  arma::uvec ranks_c1;
  double pseudo_inf = std::numeric_limits<double>::max();

  for (int j = 0; j < N; j++) {
    // Extract the j+1 "least" significant statistics
    t_star_sub = t_star.cols(ranks.subvec(N - 1 - j, N - 1));
    // Sort them within each bootstrap replication
    sorted_t_star = sort(t_star_sub, "ascend", 1);
    // Sort them along the first column
    ranks_c1 = sort_index(sorted_t_star.col(0));
    cv_star = join_cols(sorted_t_star.rows(ranks_c1), pseudo_inf * ones(1, j + 1));

    if (j == 0) {
      cv_fdr(0, N-1) = cv_star(std::min(B, int(N * level * B)) - (N * level <= 1), 0);
    } else if ((j > 0) & (j < (N-1))) {
      // Count number of non-rejections in the j "least-significant" statistics, add one supperflous column for easier coding
      noreject_jplus = ones(B + 1, j);
      noreject_jplus.elem( find(cv_star.cols(1, j) <= repelem(cv_fdr(0, span(N - j, N - 1) ), B + 1, 1)) ).zeros();
      noreject_jplus = join_rows(noreject_jplus, ones(B + 1));
      // First case: no rejection for the 1st among j statistics (set everything else to 1 as it doesn't count)
      noreject = ones(B + 1, j + 1);
      noreject.col(0) = noreject_jplus.col(0);
      // FDR estimate based on first case
      noreject_prod = prod(noreject, 1);
      FDR_est = noreject_prod/(N - (j + 1) + 1);
      // Loop over the cases of the first non-rejection in the last (j-1) statistics. The last case means all (j-1) rejections.
      for (int no_r = 1; no_r < (j + 1); no_r ++){
        // Set the (no.r-1)-th statistic to rejection
        noreject.col(no_r - 1) = 1 - noreject.col(no_r - 1);
        // Add the (r+1)-th non-rejection
        noreject.col(no_r) = noreject_jplus.col(no_r);
        // FDR estimate based on (no.r - 1) rejections, added to previous cases
        noreject_prod = prod(noreject, 1);
        FDR_est += (no_r + 1) * noreject_prod / (N - j + no_r);
      }
      cv_fdr(0, N-j-1) = cv_star(sum(cumsum(FDR_est / B) <= level) - 1, 0);
    } else if (j == (N - 1)) {
      cv_fdr(0, 0) = cv_star(int(level * B) - 1, 0);
    }
  }

  cv_fdr.reshape(N, 1);
  arma::uvec rej_H0 = cumprod((test_i.elem(ranks) < cv_fdr));
  rej_H0.elem(ranks) = rej_H0;
  const int p_hat = sum(rej_H0);

  arma::mat FDR_Tests = zeros(std::min(N, p_hat + 1), 3);
  FDR_Tests.col(0) = linspace(1, N, N).elem(ranks.subvec(0, std::min(N - 1, p_hat)));
  FDR_Tests.col(1) = test_i.elem(ranks.subvec(0, std::min(N - 1, p_hat)));
  FDR_Tests.col(2) = cv_fdr.submat(0, 0, std::min(N - 1, p_hat), 0);

  return Rcpp::List::create(
    Rcpp::Named("rej_H0") = rej_H0,
    Rcpp::Named("FDR_Tests") =  FDR_Tests,
    Rcpp::Named("ranks") = ranks
  );
}
