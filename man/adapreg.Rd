\name{adapreg}
\alias{adapreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Adaptive Regression}
\description{
  Non-parametric adaptive regression method for diffusion map basis.
}
\usage{
adapreg(D,y,mmax=min(50,length(y)),fold=NULL,nfolds=10,nrep=5)
}
\arguments{
  \item{D}{n-by-n pairwise distance matrix for a data set with n points,
    or alternatively output from the dist() function}.
  \item{y}{vector of responses to model}
  \item{mmax}{maximum model size to consider}
  \item{fold}{vector of integers of size n specifying the k-fold
    cross-validation allocation.  Default does nfolds-fold CV by
    sample(1:nfolds,length(y),replace=T)}
  \item{nfolds}{number of folds to do CV.  If fold is supplied, nfolds
    is ignored}
  \item{nrep}{number of times optimization algorithm is run (with random
    initializations).  Higher nrep allows algorithm to avoid getting stuck
    in local minima}
}
\details{
  Fits an adaptive regression model leaving as free parameters both
  the diffusion map localness parameter, epsilon, and the size
  of the regression model, m.  The adaptive regression
  model is the expansion of the response function on the first m
  diffusion map basis functions.

  This routine searches for the optimal (epsilon,m) by
  minimizing the cross-validation risk (CV MSE) of the regression
  estimate.  The function uses \code{\link{optimize}} to search over an
  appropriate range of epsilon and calls the function
  \code{\link{adapreg.m}} to find the optimal m for each epsilon.

  Default uses 10-fold cross-validation to choose the optimal model
  size.  User may also supply a vector 
  of fold allocations.  For instance,  sample(1:10,length(y),replace=T)
  does 10-fold CV while 1:length(y) performs leave-one-out CV.
}
\value{
  The returned value is a list with components

  \item{mincvrisk}{minimum cross-validation risk for the adaptive
    regression model for the given epsilon}
  \item{mopt}{size of the optimal regression model.  If mopt == mmax, it
    is advised to increase mmax.}  
  \item{epsopt}{optimal value of epsilon used in diffusion map construction}
  \item{y.hat}{predictions of the response, y-hat, for the optimal
    model}
  \item{coeff}{coefficients of the optimal model}
}
\references{
  Richards, J. W., Freeman, P. E., Lee, A. B., and Schafer, C. M., (2009),
  ApJ, 691, 32
  }
\author{ Joseph Richards \email{jwrichar@stat.cmu.edu} }

\seealso{ \code{\link{diffuse}},\code{\link{adapreg.m}}}

\examples{

## trig function on circle
t=seq(-pi,pi,.01)
x=cbind(cos(t),sin(t))
y = cos(3*t) + rnorm(length(t),0,.1)
tcol = topo.colors(32)
colvec = floor((y-min(y))/(max(y)-min(y))*32); colvec[colvec==0] = 1
scatterplot3d(x[,1],x[,2],y,color=tcol[colvec],pch=20,main="Cosine function supported on circle",angle=55,cex.main=2,col.axis="gray",cex.symbols=2,cex.lab=2,xlab=expression("x"[1]),ylab=expression("x"[2]),zlab="y")

D = as.matrix(dist(x))
# do 10-fold cross-validation to optimize (epsilon, m):
AR = adapreg(D,y)
print(paste("optimal model size:",AR$mopt,"; optimal epsilon:",round(AR$epsopt,4),"; min. CV risk:",round(AR$mincvrisk,5)))
plot(y,AR$y.hat,ylab=expression(hat("y")),cex.lab=1.5,cex.main=1.5,main="Predictions")
abline(0,1,col=2,lwd=2)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line


