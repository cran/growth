\name{potthoff}
\title{Potthoff and Roy growth curve model}
\usage{
potthoff(response, x, ccov=NULL, times=NULL, torder=0, orthogonal=T)
}
\alias{potthoff}
\alias{coefficients.potthoff}
\alias{deviance.potthoff}
\alias{print.potthoff}
\alias{residuals.potthoff}
\arguments{
\item{response}{A matrix or dataframe of response values.}
\item{x}{A matrix defining the complete intersubject differences.}
\item{ccov}{A matrix of columns of the baseline covariates to be
actually fitted, with one row per individual.}
\item{times}{A vector of unequally spaced times when they are the same
for all individuals. Not necessary if equally spaced.}
\item{torder}{Order of the polynomial in time to be fitted. If
non-numeric, the full model in time is fitted.}
\item{orthogonal}{If TRUE, uses orthogonal polynomials for time,
otherwise only centres times at their mean.}
}
\description{
\code{potthoff} fits the Potthoff and Roy repeated measurements growth
curve model with unstructured covariance matrix to completely balanced
data.
}
\value{
A list of class potthoff is returned. 
}
\author{J.K. Lindsey}
\seealso{
\code{\link{carma}}, \code{\link{elliptic}}, \code{\link{lm}}.
}
\examples{
y <- matrix(rnorm(40),ncol=5)
x <- gl(2,4)
# full model with treatment effect
potthoff(y, ~x, torder="f", ccov=~x)
# no time trend with treatment effect
potthoff(y, ~x, torder=0, ccov=~x)
# quadratic time with treatment effect
potthoff(y, ~x, torder=2, ccov=~x)
# full model without treatment effect
potthoff(y, ~x, torder="f")
# linear time without treatment effect
potthoff(y, ~x, torder=1)
}