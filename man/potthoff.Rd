% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/potthoff.r
\name{potthoff}
\alias{potthoff}
\title{Potthoff and Roy Growth Curve Model}
\usage{
potthoff(response, x = NULL, ccov = NULL, times = NULL, torder = 0,
  orthogonal = TRUE)
}
\arguments{
\item{response}{A matrix or dataframe of response values.}

\item{x}{A matrix defining the complete intersubject differences or a
Wilkinson and Rogers formula that will create one.}

\item{ccov}{A matrix of columns of the baseline covariates to be actually
fitted, with one row per individual or a W&R formula that will create one.}

\item{times}{A vector of unequally spaced times when they are the same for
all individuals. Not necessary if equally spaced.}

\item{torder}{Order of the polynomial in time to be fitted. If non-numeric,
the full model in time is fitted.}

\item{orthogonal}{If TRUE, uses orthogonal polynomials for time, otherwise
only centres times at their mean.}
}
\value{
A list of class potthoff is returned.
}
\description{
\code{potthoff} fits the Potthoff and Roy repeated measurements growth curve
model with unstructured covariance matrix to completely balanced data.
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
\seealso{
\code{\link[growth]{carma}}, \code{\link[growth]{elliptic}},
\code{\link{lm}}.
}
\author{
J.K. Lindsey
}
\keyword{models}
