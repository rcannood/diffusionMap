\name{adapreg.m}
\alias{adapreg.m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Adaptive Regression}
\description{
  Non-parametric adaptive regression method for diffusion map basis.
}
\usage{
adapreg.m(epsilon,D,y,mmax=min(50,length(y)),fold=NULL,nfolds=10,objfun=FALSE)
}
\arguments{
  \item{epsilon}{diffusion map kernel parameter}
  \item{D}{n-by-n pairwise distance matrix for a data set with n points,
    or alternatively output from the dist() function}.
  \item{y}{vector of responses to model}
  \item{mmax}{maximum model size to consider}
  \item{fold}{vector of integers of size n specifying the k-fold
    cross-validation allocation.  Default does nfolds-fold CV by
    sample(1:nfolds,length(y),replace=T)}
  \item{nfolds}{number of folds to do CV.  If fold is supplied, nfolds
    is ignored}
  \item{objfun}{if the function is to be passed into an optimization routine
    (such as minimize()), then this needs to be set to TRUE, so that
    only the minimal CV risk is returned}
}
\details{
  Fits an adaptive regression model using the estimated diffusion map
  coordinates of a data set, while holding epsilon fixed and optimizing
  over m.  The adaptive regression
  model is the expansion of the response function on the first m
  diffusion map basis functions.

  For a given epsilon value, this routine finds the optimal m by
  minimizing the cross-validation risk (CV MSE) of the regression
  estimate.  To optimize over (epsilon,m), use the function
  \code{\link{adapreg}}.

  Default uses 10-fold cross-validation to choose the optimal model
  size.  User may also supply a vector 
  of fold allocations.  For instance,  sample(1:10,length(y),replace=T)
  does 10-fold CV while 1:length(y) does leave-one-out CV.
}
\value{
  The returned value is a list with components

  \item{mincvrisk}{minimum cross-validation risk for the adaptive
    regression model for the given epsilon}
  \item{mopt}{size of the optimal regression model.  If mopt equals mmax, it
    is advised to increase mmax.}
  \item{cvrisk}{vector of CV risk estimates for model sizes from 1:mmax}
  \item{epsilon}{value of epsilon used in diffusion map construction}
  \item{y.hat}{predictions of the response, y-hat, for
    the optimal model}
  \item{coeff}{coefficients of the optimal model}
    
  If objfun is set to TRUE, then the returned value is the minimum
  cross-validation risk for the adaptive regression model for the given
  epsilon. 
}
\references{
  Richards, J. W., Freeman, P. E., Lee, A. B., and Schafer, C. M., (2009),
  ApJ, 691, 32
  }
\author{ Joseph Richards \email{joeyrichar@gmail.com} }

\seealso{ \code{\link{diffuse}},\code{\link{adapreg}} }

\examples{

library(scatterplot3d)
## trig function on circle
t=seq(-pi,pi,.01)
x=cbind(cos(t),sin(t))
y = cos(3*t) + rnorm(length(t),0,.1)
tcol = topo.colors(32)
colvec = floor((y-min(y))/(max(y)-min(y))*32); colvec[colvec==0] = 1
scatterplot3d(x[,1],x[,2],y,color=tcol[colvec],pch=20,
  main="Cosine function supported on circle",angle=55,
  cex.main=2,col.axis="gray",cex.symbols=2,cex.lab=2,
  xlab=expression("x"[1]),ylab=expression("x"[2]),zlab="y")

D = as.matrix(dist(x))
# leave-one-out cross-validation:
AR = adapreg.m(.01,D,y,fold=1:length(y))
print(paste("optimal model size:",AR$mopt,"; min. CV risk:",
  round(AR$mincvrisk,4)))
par(mfrow=c(2,1),mar=c(5,5,4,1))
plot(AR$cvrisks,typ='b',xlab="Model size",ylab="CV risk",
  cex.lab=1.5,cex.main=1.5,main="CV risk estimates")
plot(y,AR$y.hat,ylab=expression(hat("y")),cex.lab=1.5,cex.main=1.5,
  main="Predictions")
abline(0,1,col=2,lwd=2)



## swiss roll data
N=2000
t = (3*pi/2)*(1+2*runif(N));  height = runif(N);
X = cbind(t*cos(t), height, t*sin(t))
X = scale(X) + matrix(rnorm(N*3,0,0.05),N,3)
tcol = topo.colors(32)
colvec = floor((t-min(t))/(max(t)-min(t))*32); colvec[colvec==0] = 1
scatterplot3d(X,pch=18,color=tcol[colvec],xlab=expression("x"[1]),
  ylab=expression("x"[2]),zlab=expression("x"[3]),cex.lab=1.5,
  main="Swiss Roll, Noise = 0.05",cex.main=1.5,xlim=c(-2,2),
  ylim=c(-2,2),zlim=c(-2,2),col.axis="gray")

D = as.matrix(dist(X))
# 10-fold cross-validation:
AR = adapreg.m(.2,D,t,mmax=25,nfolds=5)
print(paste("optimal model size:",AR$mopt,"; min. CV risk:",
  round(AR$mincvrisk,4)))
par(mfrow=c(2,1),mar=c(5,5,4,1))
plot(AR$cvrisks,typ='b',xlab="Model size",ylab="CV risk",
  cex.lab=1.5,cex.main=1.5,main="CV risk estimates")
plot(t,AR$y.hat,ylab=expression(hat("t")),cex.lab=1.5,cex.main=1.5,
  main="Predictions")
abline(0,1,col=2,lwd=2)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line


