// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// workNat
Eigen::VectorXd workNat(const Eigen::Map<Eigen::VectorXd>& para_tilde, const bool& LEVIER, const int& Model_type, const Rcpp::Nullable<Rcpp::Function>& fixed_pars, const Rcpp::Nullable<Rcpp::Function>& fixed_values);
RcppExport SEXP _MDSV_workNat(SEXP para_tildeSEXP, SEXP LEVIERSEXP, SEXP Model_typeSEXP, SEXP fixed_parsSEXP, SEXP fixed_valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type para_tilde(para_tildeSEXP);
    Rcpp::traits::input_parameter< const bool& >::type LEVIER(LEVIERSEXP);
    Rcpp::traits::input_parameter< const int& >::type Model_type(Model_typeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::Function>& >::type fixed_pars(fixed_parsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::Function>& >::type fixed_values(fixed_valuesSEXP);
    rcpp_result_gen = Rcpp::wrap(workNat(para_tilde, LEVIER, Model_type, fixed_pars, fixed_values));
    return rcpp_result_gen;
END_RCPP
}
// natWork
Eigen::VectorXd natWork(const Eigen::Map<Eigen::VectorXd>& para, const bool& LEVIER, const int& Model_type);
RcppExport SEXP _MDSV_natWork(SEXP paraSEXP, SEXP LEVIERSEXP, SEXP Model_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type para(paraSEXP);
    Rcpp::traits::input_parameter< const bool& >::type LEVIER(LEVIERSEXP);
    Rcpp::traits::input_parameter< const int& >::type Model_type(Model_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(natWork(para, LEVIER, Model_type));
    return rcpp_result_gen;
END_RCPP
}
// volatilityVector
Eigen::VectorXd volatilityVector(const Eigen::VectorXd& para, const int& K, const int& N);
RcppExport SEXP _MDSV_volatilityVector(SEXP paraSEXP, SEXP KSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type para(paraSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(volatilityVector(para, K, N));
    return rcpp_result_gen;
END_RCPP
}
// probapi
Eigen::VectorXd probapi(const double& omega, const int& K, const int& N);
RcppExport SEXP _MDSV_probapi(SEXP omegaSEXP, SEXP KSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(probapi(omega, K, N));
    return rcpp_result_gen;
END_RCPP
}
// P
Eigen::MatrixXd P(const Eigen::VectorXd& para, const int& K, const int& N);
RcppExport SEXP _MDSV_P(SEXP paraSEXP, SEXP KSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type para(paraSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(P(para, K, N));
    return rcpp_result_gen;
END_RCPP
}
// levierVolatility1
Rcpp::List levierVolatility1(const NumericVector& ech, const Eigen::VectorXd& para, const int& Nl, const int& Model_type);
RcppExport SEXP _MDSV_levierVolatility1(SEXP echSEXP, SEXP paraSEXP, SEXP NlSEXP, SEXP Model_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type ech(echSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type para(paraSEXP);
    Rcpp::traits::input_parameter< const int& >::type Nl(NlSEXP);
    Rcpp::traits::input_parameter< const int& >::type Model_type(Model_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(levierVolatility1(ech, para, Nl, Model_type));
    return rcpp_result_gen;
END_RCPP
}
// levierVolatility
Rcpp::List levierVolatility(const NumericVector& ech, const Eigen::VectorXd& para, const int& Nl, const int& Model_type);
RcppExport SEXP _MDSV_levierVolatility(SEXP echSEXP, SEXP paraSEXP, SEXP NlSEXP, SEXP Model_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type ech(echSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type para(paraSEXP);
    Rcpp::traits::input_parameter< const int& >::type Nl(NlSEXP);
    Rcpp::traits::input_parameter< const int& >::type Model_type(Model_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(levierVolatility(ech, para, Nl, Model_type));
    return rcpp_result_gen;
END_RCPP
}
// logLik
double logLik(const Eigen::Map<Eigen::VectorXd>& para_tilde, const Eigen::MatrixXd& ech, const int& Model_type, const bool& LEVIER, const int& K, const int& N, const int& Nl, const Rcpp::Nullable<Rcpp::Function>& fixed_pars, const Rcpp::Nullable<Rcpp::Function>& fixed_values, const std::string& dis);
RcppExport SEXP _MDSV_logLik(SEXP para_tildeSEXP, SEXP echSEXP, SEXP Model_typeSEXP, SEXP LEVIERSEXP, SEXP KSEXP, SEXP NSEXP, SEXP NlSEXP, SEXP fixed_parsSEXP, SEXP fixed_valuesSEXP, SEXP disSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type para_tilde(para_tildeSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type ech(echSEXP);
    Rcpp::traits::input_parameter< const int& >::type Model_type(Model_typeSEXP);
    Rcpp::traits::input_parameter< const bool& >::type LEVIER(LEVIERSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type Nl(NlSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::Function>& >::type fixed_pars(fixed_parsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::Function>& >::type fixed_values(fixed_valuesSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type dis(disSEXP);
    rcpp_result_gen = Rcpp::wrap(logLik(para_tilde, ech, Model_type, LEVIER, K, N, Nl, fixed_pars, fixed_values, dis));
    return rcpp_result_gen;
END_RCPP
}
// logLik2
Rcpp::List logLik2(const Eigen::MatrixXd& ech, const Eigen::VectorXd& para, const int& Model_type, const bool& LEVIER, const int& K, const int& N, const double& r, const int& t, const int& Nl, const std::string& dis);
RcppExport SEXP _MDSV_logLik2(SEXP echSEXP, SEXP paraSEXP, SEXP Model_typeSEXP, SEXP LEVIERSEXP, SEXP KSEXP, SEXP NSEXP, SEXP rSEXP, SEXP tSEXP, SEXP NlSEXP, SEXP disSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type ech(echSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type para(paraSEXP);
    Rcpp::traits::input_parameter< const int& >::type Model_type(Model_typeSEXP);
    Rcpp::traits::input_parameter< const bool& >::type LEVIER(LEVIERSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const double& >::type r(rSEXP);
    Rcpp::traits::input_parameter< const int& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const int& >::type Nl(NlSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type dis(disSEXP);
    rcpp_result_gen = Rcpp::wrap(logLik2(ech, para, Model_type, LEVIER, K, N, r, t, Nl, dis));
    return rcpp_result_gen;
END_RCPP
}
// levierVolatilityMat
NumericMatrix levierVolatilityMat(const NumericMatrix& ech, const NumericMatrix Levier, const NumericVector& para, const int& Nl, const int& Model_type);
RcppExport SEXP _MDSV_levierVolatilityMat(SEXP echSEXP, SEXP LevierSEXP, SEXP paraSEXP, SEXP NlSEXP, SEXP Model_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type ech(echSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Levier(LevierSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type para(paraSEXP);
    Rcpp::traits::input_parameter< const int& >::type Nl(NlSEXP);
    Rcpp::traits::input_parameter< const int& >::type Model_type(Model_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(levierVolatilityMat(ech, Levier, para, Nl, Model_type));
    return rcpp_result_gen;
END_RCPP
}
// R_hat
Rcpp::List R_hat(const int& H, const Eigen::Map<Eigen::VectorXd>& ech, const Eigen::Map<Eigen::MatrixXi>& MC_sim, const Eigen::Map<Eigen::MatrixXd>& z_t, NumericMatrix Levier, const NumericVector& sig, const NumericVector& para, const int& Model_type, const int& N, const int& Nl);
RcppExport SEXP _MDSV_R_hat(SEXP HSEXP, SEXP echSEXP, SEXP MC_simSEXP, SEXP z_tSEXP, SEXP LevierSEXP, SEXP sigSEXP, SEXP paraSEXP, SEXP Model_typeSEXP, SEXP NSEXP, SEXP NlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type ech(echSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXi>& >::type MC_sim(MC_simSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type z_t(z_tSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Levier(LevierSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sig(sigSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type para(paraSEXP);
    Rcpp::traits::input_parameter< const int& >::type Model_type(Model_typeSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type Nl(NlSEXP);
    rcpp_result_gen = Rcpp::wrap(R_hat(H, ech, MC_sim, z_t, Levier, sig, para, Model_type, N, Nl));
    return rcpp_result_gen;
END_RCPP
}
// f_sim
Rcpp::List f_sim(const int& H, const Eigen::Map<Eigen::VectorXd>& sig, const Eigen::Map<Eigen::VectorXd>& pi_0, const Eigen::Map<Eigen::MatrixXd>& matP, const double& varphi, const double& xi, const double& shape, const double& delta1, const double& delta2);
RcppExport SEXP _MDSV_f_sim(SEXP HSEXP, SEXP sigSEXP, SEXP pi_0SEXP, SEXP matPSEXP, SEXP varphiSEXP, SEXP xiSEXP, SEXP shapeSEXP, SEXP delta1SEXP, SEXP delta2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type sig(sigSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type pi_0(pi_0SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type matP(matPSEXP);
    Rcpp::traits::input_parameter< const double& >::type varphi(varphiSEXP);
    Rcpp::traits::input_parameter< const double& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const double& >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< const double& >::type delta1(delta1SEXP);
    Rcpp::traits::input_parameter< const double& >::type delta2(delta2SEXP);
    rcpp_result_gen = Rcpp::wrap(f_sim(H, sig, pi_0, matP, varphi, xi, shape, delta1, delta2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MDSV_workNat", (DL_FUNC) &_MDSV_workNat, 5},
    {"_MDSV_natWork", (DL_FUNC) &_MDSV_natWork, 3},
    {"_MDSV_volatilityVector", (DL_FUNC) &_MDSV_volatilityVector, 3},
    {"_MDSV_probapi", (DL_FUNC) &_MDSV_probapi, 3},
    {"_MDSV_P", (DL_FUNC) &_MDSV_P, 3},
    {"_MDSV_levierVolatility1", (DL_FUNC) &_MDSV_levierVolatility1, 4},
    {"_MDSV_levierVolatility", (DL_FUNC) &_MDSV_levierVolatility, 4},
    {"_MDSV_logLik", (DL_FUNC) &_MDSV_logLik, 10},
    {"_MDSV_logLik2", (DL_FUNC) &_MDSV_logLik2, 10},
    {"_MDSV_levierVolatilityMat", (DL_FUNC) &_MDSV_levierVolatilityMat, 5},
    {"_MDSV_R_hat", (DL_FUNC) &_MDSV_R_hat, 10},
    {"_MDSV_f_sim", (DL_FUNC) &_MDSV_f_sim, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_MDSV(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
