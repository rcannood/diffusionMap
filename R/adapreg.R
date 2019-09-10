#' Adaptive Regression
#'
#' Non-parametric adaptive regression method for diffusion map basis.
#'
#' Fits an adaptive regression model leaving as free parameters both the
#' diffusion map localness parameter, epsilon, and the size of the regression
#' model, m.  The adaptive regression model is the expansion of the response
#' function on the first m diffusion map basis functions.
#'
#' This routine searches for the optimal (epsilon,m) by minimizing the
#' cross-validation risk (CV MSE) of the regression estimate.  The function
#' uses [optimize()] to search over an appropriate range of epsilon
#' and calls the function [adapreg.m()] to find the optimal m for
#' each epsilon.
#'
#' Default uses 10-fold cross-validation to choose the optimal model size.
#' User may also supply a vector of fold allocations.  For instance,
#' sample(1:10,length(y),replace=T) does 10-fold CV while 1:length(y) performs
#' leave-one-out CV.
#'
#' @param D n-by-n pairwise distance matrix for a data set with n points, or
#' alternatively output from the dist() function
#' @param y vector of responses to model
#' @param mmax maximum model size to consider
#' @param fold vector of integers of size n specifying the k-fold
#' cross-validation allocation.  Default does nfolds-fold CV by
#' sample(1:nfolds,length(y),replace=T)
#' @param nfolds number of folds to do CV.  If fold is supplied, nfolds is
#' ignored
#' @param nrep number of times optimization algorithm is run (with random
#' initializations).  Higher nrep allows algorithm to avoid getting stuck in
#' local minima
#' @return The returned value is a list with components
#'
#' \item{mincvrisk}{minimum cross-validation risk for the adaptive regression
#' model for the given epsilon} \item{mopt}{size of the optimal regression
#' model.  If mopt == mmax, it is advised to increase mmax.}
#' \item{epsopt}{optimal value of epsilon used in diffusion map construction}
#' \item{y.hat}{predictions of the response, y-hat, for the optimal model}
#' \item{coeff}{coefficients of the optimal model}
#' @seealso [diffuse()],[adapreg.m()]
#' @references Richards, J. W., Freeman, P. E., Lee, A. B., and Schafer, C. M., (2009), ApJ, 691, 32
#' @keywords multivariate nonparametric
#'
#' @importFrom stats optimize runif
#' @export
#' @examples
#' library(scatterplot3d)
#' ## trig function on circle
#' t=seq(-pi,pi,.01)
#' x=cbind(cos(t),sin(t))
#' y = cos(3*t) + rnorm(length(t),0,.1)
#' tcol = topo.colors(32)
#' colvec = floor((y-min(y))/(max(y)-min(y))*32); colvec[colvec==0] = 1
#' scatterplot3d(x[,1],x[,2],y,color=tcol[colvec],pch=20,
#'   main="Cosine function supported on circle",angle=55,
#'   cex.main=2,col.axis="gray",cex.symbols=2,cex.lab=2,
#'   xlab=expression("x"[1]),ylab=expression("x"[2]),zlab="y")
#'
#' D = as.matrix(dist(x))
#' # do 10-fold cross-validation to optimize (epsilon, m):
#' AR = adapreg(D,y, mmax=5,nfolds=2,nrep=2)
#' print(paste("optimal model size:",AR$mopt,"; optimal epsilon:",
#'   round(AR$epsopt,4),"; min. CV risk:",round(AR$mincvrisk,5)))
#' plot(y,AR$y.hat,ylab=expression(hat("y")),cex.lab=1.5,cex.main=1.5,
#'   main="Predictions")
#' abline(0,1,col=2,lwd=2)
adapreg = function(D,y,mmax=min(50,length(y)),fold=NULL,nfolds=10,nrep=5){

  # set up folds for CV
  if(length(fold)!=length(y)){ #defaults to nfolds-fold CV
    fold = sample(1:nfolds,length(y),replace=T)
  }
  nfolds = length(table(fold))

  print(paste("Doing ",nfolds,"-fold cross-validation",sep=""))
  objtmp = Inf
  for(kk in 1:nrep){
    popt = optimize(adapreg.m,lower=0,upper=epsilonCompute(D,runif(1,.05,.2)),D=D,y=y,mmax=mmax,fold=fold, objfun=TRUE) # optimize CV risk over epsilon

    if(popt$objective < objtmp){
      objtmp = popt$objective
      mintmp = popt$minimum
    }
  }
  ARopt = adapreg.m(mintmp,D,y,mmax=mmax,fold=fold)

  return(list(mincvrisk = objtmp, mopt = ARopt$mopt, epsopt = mintmp, y.hat= ARopt$y.hat,  coeff = ARopt$coeff))

}
