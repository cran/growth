\name{carma}
\title{Continuous ARMA for Unequally Spaced Repeated
Measurements}
\usage{
carma(response, ccov=NULL, times=NULL, torder=0, interaction,
	transform="identity", arma=c(0,0,0), parma=NULL, pre=NULL,
	position=NULL, iopt=T, resid=T, delta=NULL, print.level=0,
	iterlim=100, typsiz=abs(p), ndigit=10, gradtol=0.00001,
	fscale=1, stepmax=10*sqrt(p\%*\%p), steptol=0.00001)
}
\alias{carma}
\alias{coefficients.carma}
\alias{deviance.carma}
\alias{fitted.carma}
\alias{residuals.carma}
\alias{print.carma}
\arguments{
\item{response}{A list of two column matrices with response values and
times for each individual, one matrix or dataframe of response
values, or an object of either class, response (created by
\code{\link{restovec}}) or repeated (created by \code{\link{rmna}}).}
\item{ccov}{A matrix of columns of baseline covariates with one row per
individual, a model formula using vectors of the same size, or an
object of class, tccov (created by \code{\link{tcctomat}}). Ignored if
response has class, repeated.}
\item{times}{When response is a matrix, a vector of possibly unequally
spaced times when they are the same for all individuals or a matrix of
times. Not necessary if equally spaced. Ignored if response has class,
response or repeated.}
\item{torder}{Order of the polynomial in time to be fitted.}
\item{interaction}{Vector indicating order of interactions of covariates
with time.}
\item{transform}{Transformation of the response variable: \code{identity},
\code{exp}, \code{square}, \code{sqrt}, or \code{log}.}
\item{arma}{Vector of three values: order of AR, order of MA, binary
indicator for presence of measurement error. Not required for an AR(1)
if an initial estimate is supplied. If only one value is supplied, it
is assumed to be the order of the AR.}
\item{parma}{Initial estimates of ARMA parameters. For example, with
\code{arma=c(1,0,0)}, an AR(1), the parameter is
\code{parma[1]=log(theta)}, where \code{theta} is the positive,
continuous time autoregressive coefficient.  The finite step
autoregression coefficient for a step of length \code{delta} is then
\code{alpha=exp(-delta*theta)}
i.e. \code{alpha=exp(-delta*exp(parma[1]))}.}
\item{pre}{Initial estimates of random effect parameters.}
\item{position}{Two column matrix with rows giving index positions of
random effects in the covariance matrix.}
\item{iopt}{TRUE if optimization should be performed.}
\item{resid}{TRUE if residuals to be calculated.}
\item{delta}{Scalar or vector giving the unit of measurement
for each response value, set to unity by default. For example, if a
response is measured to two decimals, delta=0.01. Ignored if response
has class, response or repeated.}
\item{others}{Arguments controlling \code{\link{nlm}}.}
}
\description{
\code{carma} is designed to handle a polynomial within subject
design matrix with unequally spaced observations which can be at
different times for different subjects. The origin of time is taken
as the mean time of all the subjects. The within subject errors are
assumed to be independent Gaussian or have a continuous time ARMA(p,q)
Gaussian structure with the option to include measurement error.
The between subject random coefficients are assumed to have an
arbitrary covariance matrix. The fixed effect design matrix is a
polynomial of equal or higher order than the within subject design
matrix. This matrix can be augmented by covariates multiplied by
polynomial design matrices of any order up to the order of the first
partition of the design matrix. The method is based on exact maximum
likelihood using the Kalman filter to calculate the likelihood.

For clustered (non-longitudinal) data, where only random effects will
be fitted, the `times' may be any strictly increasing sequence
distinguishing the responses on an individual.

Marginal and individual profiles can be plotted using
\code{\link{plot.profile}} and \code{\link{plot.iprofile}} and
residuals with \code{\link{plot.residuals}}.
}
\value{
A list of class \code{carma} is returned. For any ARMA of order
superior to an AR(1), the (complex) roots of the characteristic
equation are printed out; see Jones and Ackerson (1991) for their use
in calculation of the covariance function.
}
\references{
Jones, R. H. and Ackerson, L. M. (1991) Serial correlation in
unequally spaced longitudinal data. Biometrika, 77, 721-731.

Jones, R.H. (1993) Longitudinal Data Analysis with Serial Correlation:
A State-space Approach. Chapman and Hall
}
\author{R.H. Jones and J.K. Lindsey}
\seealso{
\code{\link{elliptic}}, \code{\link{gar}}, \code{\link{glmm}},
\code{\link{gnlmm}}, \code{\link{kalseries}}, \code{\link{plot.profile}},
\code{\link{plot.iprofile}}, \code{\link{plot.residuals}},
\code{\link{potthoff}}, \code{\link{read.list}}, \code{\link{restovec}},
\code{\link{rmna}}, \code{\link{tcctomat}}, \code{\link{tvctomat}}.
}
\examples{
y <- matrix(rnorm(40),ncol=5)
x1 <- gl(2,4)
x2 <- gl(2,1,8)
# independence with time trend
carma(y, ccov=~x1, torder=2)
# AR(1)
carma(y, ccov=~x1, torder=2, arma=c(1,0,0), parma=-0.5)
carma(y, ccov=~x1, torder=3, interact=3, arma=c(1,0,0), parma=-1)
# ARMA(2,1)
carma(y, ccov=~x1+x2, interact=c(2,0), torder=3,arma=c(2,1,0),
	parma=c(0.3,2,0.7))
# random intercept
carma(y, ccov=~x1+x2, interact=c(2,0), torder=3, pre=-0.4, 
	position=c(1,1))
# random coefficients
carma(y, ccov=~x1+x2, interact=c(2,0), torder=3, pre=c(-0.4,0.1), 
	position=rbind(c(1,1),c(2,2)))
}
