\name{pergram}
\title{Calculate and plot a periodogram}
\usage{
pergram(y)
plot.pergram(y, add=F, lty=1, xlab="Frequency", ylab="Periodogram",
	main="Periodogram", ylim=c(0,max(po[,2])))
plot.cum.pergram(y, xlab="Frequency", ylab="Periodogram",
	main="Cumulative periodogram",
	ylim=c(0,max(cpo+1.358/(a+0.12+0.11/a))))
}
\alias{plot.pergram}
\alias{plot.cum.pergram}
\arguments{
\item{y}{A time series vector.}
\item{add}{If TRUE, adds a new periodogram to an existing plot.}
\item{others}{Graphics parameters}
}
\description{
\code{pergram} calculates the values of a periodogram,
\code{plot.pergram} plots it, and \code{plot.cum.pergram} plots the
corresponding cumulative periodogram.
}
\value{
\code{pergram} prints and returns a two-column matrix of class,
\code{pergram}, containing the periodogram.
}
\author{J.K. Lindsey}
\examples{
y <- rnorm(100)
print(z <- pergram(y))
plot(z)
plot.cum(z)
}