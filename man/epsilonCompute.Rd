\name{epsilonCompute}
\alias{epsilonCompute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute default diffusion map epsilon. }
\description{
  Uses the pair-wise distances to estimate a diffusion map epsilon value by
  the median p*n-th nearest neighbor
}
\usage{
epsilonCompute(D, p = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{n-by-n pairwise distance matrix for a data set with n points,
  or alternatively output from the dist() function }
  \item{p}{ distances to p*n-th nearest neighbor are used.  Default
    value is .01}
}
\details{
  Function is used as the default value in diffuse().  For 
  inference problems, it is advised that the results be optimized over epsilon.
}
\value{
  \item{epsilon }{value of epsilon to be used in diffusion map}
}
\author{ Joseph Richards \email{joeyrichar@gmail.com} }
\seealso{ \code{\link{diffuse}} }
\examples{
data(annulus)
D = dist(annulus) # use Euclidean distance
epsilonCompute(D,.005)
epsilonCompute(D,.01)
epsilonCompute(D,.05)
epsilonCompute(D,.1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line

