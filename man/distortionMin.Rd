\name{distortionMin}
\alias{distortionMin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distortion Minimization via K-means}
\description{
  Runs one K-means loop based on the diffusion coordinates of a data
  set, beginning from an initial set of cluster centers.  
}
\usage{
distortionMin(X, phi0, K, c0, epsilon = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ diffusion coordinates, each row corresponds to a data point }
  \item{phi0}{ trivial left eigenvector of Markov matrix (stationary
    distribution of Markov random walk) in diffusion map construction}
  \item{K}{number of clusters }
  \item{c0}{ initial cluster centers}
  \item{epsilon}{ stopping criterion for relative change in distortion}
}
\details{
  Used by diffusionKmeans().
}
\value{
  The returned value is a list with components

  \item{S}{ labelling from K-means loop. n-dimensional vector with
    integers between 1 and K}
  \item{c}{ K geometric centroids found by K-means}
  \item{D}{minimum of total distortion (loss function of K-means) found
    in K-means run}
  \item{DK}{n by k matrix of squared (Euclidean) distances from each
    point to every centroid} 
}
\references{
    Lafon, S., \& Lee, A., (2006), IEEE Trans. Pattern Anal. and
  Mach. Intel., 28, 1393
}
\author{ Joseph Richards \email{jwrichar@stat.cmu.edu} }

\seealso{ \code{\link{diffusionKmeans}} }
\examples{
data(annulus)
n = dim(annulus)[1]
D = dist(annulus) # use Euclidean distance
dmap = diffuse(D,0.03) # compute diffusion map 
km = distortionMin(dmap$X,dmap$phi0,2,dmap$X[sample(n,2),])
plot(annulus,col=km$S,pch=20)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line

