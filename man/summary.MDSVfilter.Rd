% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MDSVfilter.R
\name{summary.MDSVfilter}
\alias{summary.MDSVfilter}
\alias{print.summary.MDSVfilter}
\alias{print.MDSVfilter}
\title{Summarize and print MDSV Filtering}
\usage{
\method{summary}{MDSVfilter}(object, ...)

\method{print}{summary.MDSVfilter}(x, ...)

\method{print}{MDSVfilter}(x, ...)
}
\arguments{
\item{object}{An object of class \link{MDSVfilter}, output of the function \code{\link{MDSVfilter}}.}

\item{...}{Further arguments passed to or from other methods.}

\item{x}{An object of class \link{summary.MDSVfilter}, output of the function \code{\link{summary.MDSVfilter}}
or class \link{MDSVfilter} of the function \code{\link{MDSVfilter}}.}
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
    \item filtred_proba : matrix containing the filtred probabilities \eqn{\mathbb{P}(C_t=c_i\mid x_1,\dots,\x_t)} of the Markov Chain.
    \item smoothed_proba : matrix containing the smoothed probabilities \eqn{\mathbb{P}(C_t=c_i\mid x_1,\dots,\x_T)} of the Markov Chain.
    \item Marg_loglik : marginal log-likelihood corresponding to the log-likelihood of log-returns. This is only return when \eqn{ModelType = 2}.
    \item VaR : Value-at-Risk compute empirically.
}
}
\description{
Summary and print methods for the class \link{MDSVfilter} as returned by the function \code{\link{MDSVfilter}}.
}
\seealso{
For fitting \code{\link{MDSVfit}}, filtering \code{\link{MDSVfilter}}, bootstrap forecasting \code{\link{MDSVboot}} and rolling estimation and forecast \code{\link{MDSVroll}}.
}
