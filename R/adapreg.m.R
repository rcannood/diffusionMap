#' Adaptive Regression
#'
#' Non-parametric adaptive regression method for diffusion map basis.
#'
#' Fits an adaptive regression model using the estimated diffusion map
#' coordinates of a data set, while holding epsilon fixed and optimizing over
#' m.  The adaptive regression model is the expansion of the response function
#' on the first m diffusion map basis functions.
#'
#' For a given epsilon value, this routine finds the optimal m by minimizing
#' the cross-validation risk (CV MSE) of the regression estimate.  To optimize
#' over (epsilon,m), use the function [adapreg()].
#'
#' Default uses 10-fold cross-validation to choose the optimal model size.
#' User may also supply a vector of fold allocations.  For instance,
#' sample(1:10,length(y),replace=T) does 10-fold CV while 1:length(y) does
#' leave-one-out CV.
#'
#' @param epsilon diffusion map kernel parameter
#' @param D n-by-n pairwise distance matrix for a data set with n points, or
#' alternatively output from the dist() function
#' @param y vector of responses to model
#' @param mmax maximum model size to consider
#' @param fold vector of integers of size n specifying the k-fold
#' cross-validation allocation.  Default does nfolds-fold CV by
#' sample(1:nfolds,length(y),replace=T)
#' @param nfolds number of folds to do CV.  If fold is supplied, nfolds is
#' ignored
#' @param objfun if the function is to be passed into an optimization routine
#' (such as minimize()), then this needs to be set to TRUE, so that only the
#' minimal CV risk is returned
#' @return The returned value is a list with components
#'
#' \item{mincvrisk}{minimum cross-validation risk for the adaptive regression
#' model for the given epsilon} \item{mopt}{size of the optimal regression
#' model.  If mopt equals mmax, it is advised to increase mmax.}
#' \item{cvrisk}{vector of CV risk estimates for model sizes from 1:mmax}
#' \item{epsilon}{value of epsilon used in diffusion map construction}
#' \item{y.hat}{predictions of the response, y-hat, for the optimal model}
#' \item{coeff}{coefficients of the optimal model}
#'
#' If objfun is set to TRUE, then the returned value is the minimum
#' cross-validation risk for the adaptive regression model for the given
#' epsilon.
#' @seealso [diffuse()],[adapreg()]
#' @references Richards, J. W., Freeman, P. E., Lee, A. B., and Schafer, C. M.,
#' (2009), ApJ, 691, 32
#' @keywords multivariate nonparametric
#'
#' @importFrom stats lm
#'
#' @export
#' @examples
#' library(stats)
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
#' # leave-one-out cross-validation:
#' AR = adapreg.m(.01,D,y,fold=1:length(y))
#' print(paste("optimal model size:",AR$mopt,"; min. CV risk:",
#'   round(AR$mincvrisk,4)))
#' par(mfrow=c(2,1),mar=c(5,5,4,1))
#' plot(AR$cvrisks,typ='b',xlab="Model size",ylab="CV risk",
#'   cex.lab=1.5,cex.main=1.5,main="CV risk estimates")
#' plot(y,AR$y.hat,ylab=expression(hat("y")),cex.lab=1.5,cex.main=1.5,
#'   main="Predictions")
#' abline(0,1,col=2,lwd=2)
#'
#' ## swiss roll data
#' N=2000
#' t = (3*pi/2)*(1+2*runif(N));  height = runif(N);
#' X = cbind(t*cos(t), height, t*sin(t))
#' X = scale(X) + matrix(rnorm(N*3,0,0.05),N,3)
#' tcol = topo.colors(32)
#' colvec = floor((t-min(t))/(max(t)-min(t))*32); colvec[colvec==0] = 1
#' scatterplot3d(X,pch=18,color=tcol[colvec],xlab=expression("x"[1]),
#'   ylab=expression("x"[2]),zlab=expression("x"[3]),cex.lab=1.5,
#'   main="Swiss Roll, Noise = 0.05",cex.main=1.5,xlim=c(-2,2),
#'   ylim=c(-2,2),zlim=c(-2,2),col.axis="gray")
#'
#' D = as.matrix(dist(X))
#' # 10-fold cross-validation:
#' AR = adapreg.m(.2,D,t,mmax=25,nfolds=5)
#' print(paste("optimal model size:",AR$mopt,"; min. CV risk:",
#'   round(AR$mincvrisk,4)))
#' par(mfrow=c(2,1),mar=c(5,5,4,1))
#' plot(AR$cvrisks,typ='b',xlab="Model size",ylab="CV risk",
#'   cex.lab=1.5,cex.main=1.5,main="CV risk estimates")
#' plot(t,AR$y.hat,ylab=expression(hat("t")),cex.lab=1.5,cex.main=1.5,
#'   main="Predictions")
#' abline(0,1,col=2,lwd=2)
adapreg.m = function(epsilon,D,y,mmax=min(50,length(y)),fold=NULL,nfolds=10,objfun=FALSE){

  print(paste("Finding optimal model for epsilon=",round(epsilon,4)))

  dmap = diffuse(as.matrix(D),eps.val=epsilon,neigen=mmax,t=1)

  pred.y = matrix(0,length(y),mmax)

  if(length(fold)!=length(y)){ #defaults to nfolds-fold CV
    fold = sample(1:nfolds,length(y),replace=T)
  }
  nfolds = length(table(fold))

  for(ii in 1:nfolds){
    if(sum(fold==ii)>1){
        # regress y on diffusion map coordinates
      AR = lm(y[fold!=ii]~dmap$X[fold!=ii,])
      for(jj in 1:mmax){ # predictions for each size model
        pred.y[fold==ii,jj] = cbind(rep(1,sum(fold==ii)),dmap$X[fold==ii,1:jj])%*%AR$coeff[1:(jj+1)]
      }
    }
    if(sum(fold==ii)==1){ # if only 1 data point in the fold
      # regress y on diffusion map coordinates
      AR = lm(y[fold!=ii]~dmap$X[fold!=ii,])
      for(jj in 1:mmax){ # predictions for each size model
        pred.y[fold==ii,jj] = c(1,dmap$X[fold==ii,1:jj])%*%AR$coeff[1:(jj+1)]
      }
    }
  }

  # compute MSE
  risk.dmap = apply((pred.y-y)^2,2,mean)
  mopt = which.min(risk.dmap)

  ARopt = lm(y~dmap$X[,1:mopt])

  if(objfun==TRUE){
    return(min(risk.dmap))
  } else{
    return(list(mincvrisk = min(risk.dmap), mopt = mopt, cvrisks = risk.dmap, epsilon = epsilon, y.hat = ARopt$fitted.values, coeff = ARopt$coeff))
  }
}
