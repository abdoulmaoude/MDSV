% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MDSVfilter.R
\name{MDSVfilter}
\alias{MDSVfilter}
\title{MDSV Filtering}
\usage{
MDSVfilter(
  N,
  K,
  data,
  para,
  ModelType = 0,
  LEVIER = FALSE,
  calculate.VaR = TRUE,
  VaR.alpha = c(0.01, 0.05),
  dis = "lognormal"
)
}
\arguments{
\item{N}{An integer designing the number of components for the MDSV process}

\item{K}{An integer designing the number of states of each MDSV process component}

\item{data}{A univariate or bivariate data matrix. Can only be a matrix of 1 or 2 columns. If data has 2 columns, the first one has to be the log-returns and the second the realized variances.}

\item{para}{A vector of parameters use for the MDSV filtering on \code{data}. For more informations see Details.}

\item{ModelType}{An integer designing the type of model to be fit. \eqn{0} for univariate log-returns, \eqn{1} for univariate realized variances and \eqn{2} for joint log-return and realized variances.}

\item{LEVIER}{if \code{TRUE}, estime the MDSV model with leverage.}

\item{calculate.VaR}{Whether to calculate forecast Value at Risk during the estimation.}

\item{VaR.alpha}{The Value at Risk tail level to calculate.}
}
\value{
A list consisting of: 
\itemize{
    \item ModelType : type of model to be filtered.
    \item LEVIER : wheter the filter take the leverage effect into account or not.
    \item N : number of components for the MDSV process.
    \item K : number of states of each MDSV process component.
    \item data : data use for the filtering.
    \item dates : vector or names of data designing the dates.
    \item estimates : input parameters.
    \item LogLikelihood : log-likelihood of the model on the data.
    \item AIC : Akaike Information Criteria of the model on the data.
    \item BIC : Bayesian Information Criteria of the model on the data.
    \item Levier : numeric vector representing the leverage effect at each date. \code{Levier} is 1 when no leverage is detected
    \item filtred_proba : matrix containing the filtred probabilities \eqn{P(C_t=c_i | x_1,\dots, x_t)} of the Markov Chain.
    \item smoothed_proba : matrix containing the smoothed probabilities \eqn{P(C_t=c_i | x_1,\dots, x_T)} of the Markov Chain.
    \item Marg_loglik : marginal log-likelihood corresponding to the log-likelihood of log-returns. This is only return when \eqn{ModelType = 2}.
    \item VaR : Value-at-Risk compute empirically.
}
}
\description{
Method for filtering the MDSV model on log-retruns and realized variances (uniquely or jointly).
}
\details{
The MDSVfilter help to simply filter a set of data with a predefined set of parameters. Like in \code{\link{MDSVfit}}, the 
likelihood calculation is performed in \code{C++} through the \pkg{Rcpp} package. 
According to Augustyniak et al., (2021), the para consist on a vector of :
\itemize{
    \item{\eqn{\omega}}{ probability of success of the stationnary distribution of the Markov chain. (\eqn{0 < \omega < 1}).}
    \item{\eqn{a}}{ highest persistance of the component of the MDSV process. (\eqn{0 < a < 1}).}
    \item{\eqn{b}}{ decroissance rate of the persistances. (\eqn{1 < b}).}
    \item{\eqn{\sigma^2}}{ unconditionnal variance of the MDSV process.}
    \item{\eqn{\nu_0}}{ states defined parameter. (\eqn{0 < \nu_0 < 1}).}
    \item{\eqn{\gamma}}{ parameter of the realized variances innovation. This parameter is required only if \eqn{ModelType = 1} or \eqn{ModelType = 2}. (\eqn{0 < \gamma}).}
    \item{\eqn{\xi}}{ parameter of realized variances equation in joint estimation. This parameter is required only if \eqn{ModelType = 2}.}
    \item{\eqn{\phi}}{ parameter of realized variances equation in joint estimation. This parameter is required only if \eqn{ModelType = 2}.}
    \item{\eqn{\delta1}}{ parameter of realized variances equation in joint estimation. This parameter is required only if \eqn{ModelType = 2}.}
    \item{\eqn{\delta2}}{ parameter of realized variances equation in joint estimation. This parameter is required only if \eqn{ModelType = 2}.}
    \item{\eqn{l}}{ leverage effect parameter. This parameter is required only if \eqn{LEVIER = TRUE}. (\eqn{0 < l}).}
    \item{\eqn{\theta}}{ leverage effect parameter. This parameter is required only if \eqn{LEVIER = TRUE}. (\eqn{0 < \theta < 1}).}
}
The leverage effect is taken into account according to the FHMV model (see Augustyniak et al., 2019). While filtering an
univariate realized variances data, log-returns are required to add leverage effect. 
AIC and BIC are computed using the formulas : 
\itemize{
\item{AIC : }{\eqn{L - k}}
\item{BIC : }{\eqn{L - (k/2)*log(n)}}
}
where \eqn{L} is the log-likelihood, \eqn{k} is the number of parameters and \eqn{n} the number of observations in the dataset.
The \link[base]{class} of the output of this function is \code{\link{MDSVfilter}}. This class has a \link[base]{summary}, 
\link[base]{print} and \link[base]{plot} \link[utils]{methods} to summarize, print and plot the results. See 
\code{\link{summary.MDSVfilter}}, \code{\link{print.MDSVfilter}} and \code{\link{plot.MDSVfilter}} for more details.
}
\examples{
\dontrun{
# MDSV(N=2,K=3) without leverage on univariate log-returns S&P500
data(sp500)        # Data loading
N         <- 2     # Number of components
K         <- 3     # Number of states
ModelType <- 0     # Univariate log-returns
LEVIER    <- FALSE # No leverage effect

# Model estimation
out_fit   <- MDSVfit(K = K, N = N, data = sp500, ModelType = ModelType, LEVIER = LEVIER)
# Model filtering
para      <-out_fit$estimates # parameter
out_filter<- MDSVfilter(K = K, N = N, data = sp500, para = para, ModelType = ModelType, LEVIER = LEVIER)
# Summary
summary(out_filter)
# Plot
plot(out_filter)


# MDSV(N=3,K=3) with leverage on joint log-returns and realized variances NASDAQ
data(nasdaq)      # Data loading
N         <- 3    # Number of components
K         <- 3    # Number of states
ModelType <- 2    # Joint log-returns and realized variances
LEVIER    <- TRUE # No leverage effect

para      <- c(omega = 0.52, a = 0.99, b = 2.77, sigma = 1.95, v0 = 0.72, 
              xi = -0.5, varphi = 0.93, delta1 = 0.93, delta2 = 0.04, shape = 2.10,
              l = 0.78, theta = 0.876)
# Model filtering
out       <- MDSVfilter(K = K, N = N,data = nasdaq, para = para, ModelType = ModelType, LEVIER = LEVIER)
# Summary
summary(out)
# Plot
plot(out)
}

}
\references{
Augustyniak, M., Bauwens, L., & Dufays, A. (2019). A new approach to volatility modeling: the factorial hidden Markov volatility model. 
\emph{Journal of Business & Economic Statistics}, 37(4), 696-709. \url{https://doi.org/10.1080/07350015.2017.1415910}

Augustyniak, M., Dufays, A., & Maoude, K.H.A. (2021). Multifractal Discrete Stochastic Volatility.
}
\seealso{
For fitting \code{\link{MDSVfit}}, bootstrap forecasting \code{\link{MDSVboot}}, simulation \code{\link{MDSVsim}} and rolling estimation and forecast \code{\link{MDSVroll}}.
}
