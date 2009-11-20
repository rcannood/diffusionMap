\name{nystrom}
\alias{nystrom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Perform Nystrom Extension to estimate diffusion coordinates of data.}
\description{
  Given the diffusion map coordinates of a training data set, estimates the
  diffusion map coordinates of a new set of data using the
  pairwise distance matrix from the new data to the original data.
}
\usage{
nystrom(dmap,Dnew,sigma=dmap$epsilon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dmap}{a '"dmap"' object from the original data set, computed by diffuse()}
  \item{Dnew}{distance matrix between each new data point and every
    point in the training data set.  Matrix is m-by-n, where m is the
    number of data points in the new set and n is the number of training
    data points}
  \item{sigma}{scalar giving the size of the Nystrom extension kernel.
    Default uses the tuning parameter of the original diffusion map }
}
\details{
  Often, it is computationally infeasible to compute the exact diffusion
  map coordinates for large data sets.  In this case, one may use the
  exact diffusion coordinates of a training data set to extend to a new
  data set using the Nystrom approximation.

  A Gaussian kernel is used:  exp(-D(x,y)^2/sigma).
  The default value of sigma is the epsilon value used in the
  construction of the original diffusion map.  Other methods to
  select sigma, such as Algorithm 2 in Lafon, Keller, and Coifman (2006)
  have been proposed.

  The dimensionality of the diffusion map representation of the new data
  set will be the same as the dimensionality of the diffusion map
  constructed on the original data.
  
}
\value{
  The estimated diffusion coordinates for the new data, a matrix of
  dimensions m by p, where p is the dimensionality of the input
  diffusion map 
}
\references{
  Freeman, P. E., Newman, J. A., Lee, A. B., Richards, J. W., and Schafer,
  C. M. (2009), MNRAS, Volume 398, Issue 4, pp. 2012-2021

  Lafon, S., Keller, Y., and Coifman, R. R. (2006), IEEE Trans. Pattern Anal. and
  Mach. Intel., 28, 1784
}
\author{ Joseph Richards \email{jwrichar@stat.cmu.edu} }
\seealso{ \code{\link{diffuse}} }
\examples{
Norig = 1000
Next = 4000
t=runif(Norig+Next)^.7*10
al=.15;bet=.5;
x1=bet*exp(al*t)*cos(t)+rnorm(length(t),0,.1)
y1=bet*exp(al*t)*sin(t)+rnorm(length(t),0,.1)

D = as.matrix(dist(cbind(x1,y1)))
Dorig = D[1:Norig,1:Norig] # training distance matrix
DExt = D[(Norig+1):(Norig+Next),1:Norig] # new data distance matrix
# compute original diffusion map
dmap = diffuse(Dorig,neigen=2) 
 # use Nystrom extension
dmapExt = nystrom(dmap,DExt)
plot(dmapExt[,1:2],pch=8,col=2,main="Diffusion map, black = original, red = new data",xlab="1st diffusion coefficient",ylab="2nd diffusion coefficient")
points(dmap$X[,1:2],pch=19,cex=.5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line

