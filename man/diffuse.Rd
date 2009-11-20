\name{diffuse}
\alias{diffuse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute diffusion map coordinates from pair-wise
  distances.}
\description{
  Uses the pair-wise distance matrix for a data set to compute the
  diffusion map coefficients.  Computes the Markov transition
  probability matrix, and its
  eigenvalues 
  and left & right eigenvectors.  Returns a 'dmap' object.
}
\usage{
diffuse(D, eps.val = epsilonCompute(D), neigen = NULL, t = 0, maxdim = 50, delta=10^-5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{n-by-n pairwise distance matrix for a data set with n points,
  or alternatively output from the dist() function}.
  \item{eps.val}{epsilon parameter for the diffusion weight matrix,
    exp(-D$^2$/(eps.val)).  Default is to use the epsilon
    corresponding to the median distance to the 0.01*n nearest neighbor}
  \item{neigen}{ number of dimensions of final diffusion map representation.
    Default uses number of dimensions corresponding to a 95\% drop-off
    in eigenvalue multiplier.}
  \item{t}{optional time-scale parameter in the diffusion map.
    The (recommended) default uses multiscale geometry. }
  \item{maxdim}{the maximum number of diffusion map dimensions returned
    if 95\% drop-off is not attained.}
  \item{delta}{sparsity cut-off for the symmetric graph Laplacian.
    Default of 10^-5 is used.  Higher value induces more sparsity in
    Laplacian (and faster computations)}
}
\details{
  Diffusion map is a powerful tool for data parametrization that
  exploits the natural geometry of a data set.  Diffusion map uses local
  interactions between data points, propogated to larger
  scales, to construct a global representation of the data.

  The parameter eps.val controls the degree of localness in the diffusion
  weight matrix.  For most statisitical inference problems using
  diffusion map, results should be optimized over eps.val.  Generally
  a good starting point is to pick eps.val as $2*$med.knn$^2$, where
  med.knn is 
  the median distance to the kth nearest neighbor, and k is chosen 1-2\% of
  n.  The default uses 1\% of n.

  Computation of the diffusion map coordinates requires singular value
  decomposition of the normalized graph Laplacian.  This operation is
  optimized for speed by exploiting the sparseness of the graph
  Laplacian and by using ARPACK for fast matrix decomposition.
  Increasing the sparseness parameter, delta, will speed up the algorithm.
}
\value{
  The returned value is an object of 'class' 'dmap'.

 The function 'plot' is used to plot the diffusion coordinates in 1, 2,
or 3 dimensions.  The function 'print' displays the computed eigen-multipliers
and the value of epsilon used.

   An object of class 'dmap' is a list containing the
   following components:

  \item{X}{matrix of n diffusion map coordinates, entered column-wise
    (does not include the trivial coordinate)}
  \item{phi0}{ trivial left eigenvector of Markov matrix (stationary
    distribution of Markov random walk) in diffusion map construction}
  \item{eigenvals}{eigen-values of the svd of the symmetric graph Laplacian}
  \item{eigenmult}{eigen-multipliers of the diffusion map}
  \item{psi}{right eigenvectors of the Markov matrix (first row is the
    trivial right eigenvector)}
  \item{phi}{left eigenvectors of the Markov matrix (first row is the
    trivial left eigenvector)}
  \item{neigen}{number of diffusion map dimensions used}
  \item{epsilon}{the value of epsilon used}
}
\references{
  Coifman, R. R., \& Lafon, S., (2006), Appl. Comput. Harmon. Anal., 21, 5
  
  Lafon, S., \& Lee, A., (2006), IEEE Trans. Pattern Anal. and
  Mach. Intel., 28, 1393

  Richards, J. W., Freeman, P. E., Lee, A. B., Schafer, C. M., (2009),
  ApJ, 691, 32
   }
\author{ Joseph Richards \email{jwrichar@stat.cmu.edu} }
\seealso{ \code{\link{plot.dmap}},\code{\link{print.dmap}} }
\examples{

## example with noisy spiral
n=2000
t=runif(n)^.7*10
al=.15;bet=.5;
x1=bet*exp(al*t)*cos(t)+rnorm(n,0,.1)
y1=bet*exp(al*t)*sin(t)+rnorm(n,0,.1)
plot(x1,y1,pch=20,main="Noisy spiral")
D = dist(cbind(x1,y1))
dmap = diffuse(D,neigen=10) # compute diffusion map
par(mfrow=c(2,1))
plot(t,dmap$X[,1],pch=20,axes=FALSE,xlab="spiral parameter",ylab="1st diffusion coefficient")
box()
plot(1:10,dmap$eigenmult,typ='h',xlab="diffusion map dimension",ylab="eigen-multipliers")

## example with annulus data set
data(annulus)
plot(annulus,main="Annulus Data",pch=20,cex=.7)
D = dist(annulus) # use Euclidean distance
dmap = diffuse(D,eps.val=.1) # compute diffusion map & plot
print(dmap)
plot(dmap)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
