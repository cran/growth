\name{corgram}
\title{Calculate and plot a correlogram}
\usage{
corgram(y, wt=1, add=F, lty=1, ylim=NULL, xlab=NULL, ylab=NULL, main=NULL)
}
\arguments{
\item{y}{A time series vector.}
\item{wt}{Indicator vector with zeros for values to be ignored.}
\item{add}{If TRUE, adds a new correlogram to an existing plot.}
\item{others}{Graphics parameters}
}
\description{
\code{corgram} calculates the values of a correlogram and plots it.
}
\value{
\code{corgram} returns a two-column matrix containing the correlogram.
}
\author{J.K. Lindsey}
\examples{
y <- rnorm(100)
corgram(y)
}

