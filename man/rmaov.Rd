\name{rmaov}
\title{Split-plot ANOVA Model}
\usage{
rmaov(response, tvcov=NULL, ccov=NULL, analysis=T)
}
\arguments{
\item{response}{A matrix or dataframe of response values with units as
rows and repeated measures as columns.}
\item{tvcov}{A numeric vector or factor variable defining the
clusters. If there are several levels of nesting, a matrix or
dataframe with columns of such variables defining the nested clusters
starting with the highest level (that is, from slowest to fastest
varying). If not provided, each response value of a unit is assumed to
belong to a different cluster (that is, one factor with
\code{ncol(response)} levels is assumed).}
\item{ccov}{A vector or factor variable for one inter-subject
covariate or a matrix, dataframe, or list of several such variables.}
\item{analysis}{If FALSE, the design matrix is set up, but the
analysis is not performed.}
}
\description{
\code{rmaov} performs the classical balanced split-plot ANOVA, with
\code{summary} providing the table. This is the so-called repeated
measures ANOVA.

For unbalanced data, \code{\link{elliptic}} will
perform the analysis for one or two levels of nesting.
}
\value{
The fitted model is returned. 
}
\author{Ralf Goertz (ralf.goertz@uni-jena.de)}
\seealso{
\code{\link{carma}}, \code{\link{elliptic}}, \code{\link{lm}},
\code{\link{potthoff}}.
}
\examples{
# vision data for 7 individuals, with response a 7x8 matrix
# two levels of nesting: 4 levels of power for each eye
y <- matrix(rnorm(56),ncol=8)
tvc <- data.frame(eye=c(rep(1,4),rep(2,4)),power=c(1:4,1:4))
summary(rmaov(y, tvc))
}