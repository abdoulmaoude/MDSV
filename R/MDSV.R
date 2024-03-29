#' @title MDSV package
#' @description The \pkg{MDSV} package implements the Multifractal Discrete Stochastic Volatility developed in Augustyniak, et al. (2021). 
#' This MDSV model is proposed as a generalization of other high dimensional hidden markov models such as, MSM of Calvet and Fisher (2004), CDRS of Fleming and Kirby (2013), DSARV of Cordis and Kirby (2014), FHMV of Augustyniak et al. (2019).
#' To make the computations faster \pkg{MDSV} uses \code{C++} through the \pkg{Rcpp} package (Eddelbuettel et al., 2011). 
#'
#' @seealso For fitting \code{\link{MDSVfit}}, filtering \code{\link{MDSVfilter}}, bootstrap forecasting \code{\link{MDSVboot}} and rolling estimation and forecast \code{\link{MDSVroll}}.
#' 
#' @references 
#' Calvet, L. E. and Fisher, A. J. (2004). How to forecast long-run volatility: Regime switching and the estimation of multifractal processes. 
#' \emph{Journal of Financial Econometrics}, 2(1):49-83. \url{https://doi.org/10.1093/jjfinec/nbh003}
#' @references Eddelbuettel, D., Fran?ois, R., Allaire, J., Ushey, K., Kou, Q., Russel, N., ... & Bates, D., (2011).
#' \pkg{Rcpp}: Seamless \R and \code{C++} integration. \emph{Journal of Statistical Software}, 40(8), 1-18.
#' \url{https://www.jstatsoft.org/v40/i08/}
#' @references 
#' Fleming, J., & Kirby, C. (2013). Component-driven regime-switching volatility. \emph{Journal of Financial Econometrics}, 11(2), 263-301. \url{https://doi.org/10.1093/jjfinec/nbs023}
#' @references  
#' Cordis, A. S., & Kirby, C. (2014). Discrete stochastic autoregressive volatility. \emph{Journal of Banking & Finance}, 43, 160-178. \url{https://doi.org/10.1016/j.jbankfin.2014.03.020}
#' @references  
#' Augustyniak, M., Bauwens, L., & Dufays, A. (2019). A new approach to volatility modeling: the factorial hidden Markov volatility model. 
#' \emph{Journal of Business & Economic Statistics}, 37(4), 696-709. \url{https://doi.org/10.1080/07350015.2017.1415910}
#' @references 
#' Augustyniak, M., Dufays, A., & Maoude, K.H.A. (2021). Multifractal Discrete Stochastic Volatility.
#' 
#' @useDynLib MDSV, .registration = TRUE
"_PACKAGE"
NULL