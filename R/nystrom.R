#' Perform Nystrom Extension to estimate diffusion coordinates of data.
#'
#' Given the diffusion map coordinates of a training data set, estimates the
#' diffusion map coordinates of a new set of data using the pairwise distance
#' matrix from the new data to the original data.
#'
#' Often, it is computationally infeasible to compute the exact diffusion map
#' coordinates for large data sets.  In this case, one may use the exact
#' diffusion coordinates of a training data set to extend to a new data set
#' using the Nystrom approximation.
#'
#' A Gaussian kernel is used: exp(-D(x,y)^2/sigma).  The default value of sigma
#' is the epsilon value used in the construction of the original diffusion map.
#' Other methods to select sigma, such as Algorithm 2 in Lafon, Keller, and
#' Coifman (2006) have been proposed.
#'
#' The dimensionality of the diffusion map representation of the new data set
#' will be the same as the dimensionality of the diffusion map constructed on
#' the original data.
#'
#' @param dmap a '"dmap"' object from the original data set, computed by
#' diffuse()
#' @param Dnew distance matrix between each new data point and every point in
#' the training data set.  Matrix is m-by-n, where m is the number of data
#' points in the new set and n is the number of training data points
#' @param sigma scalar giving the size of the Nystrom extension kernel.
#' Default uses the tuning parameter of the original diffusion map
#' @return The estimated diffusion coordinates for the new data, a matrix of
#' dimensions m by p, where p is the dimensionality of the input diffusion map
#' @seealso [diffuse()]
#' @references Freeman, P. E., Newman, J. A., Lee, A. B., Richards, J. W., and
#' Schafer, C. M. (2009), MNRAS, Volume 398, Issue 4, pp. 2012-2021
#'
#' Lafon, S., Keller, Y., and Coifman, R. R. (2006), IEEE Trans. Pattern Anal.
#' and Mach. Intel., 28, 1784
#' @keywords multivariate nonparametric
#'
#' @export
#' @examples
#' library(stats)
#' Norig = 1000
#' Next = 4000
#' t=runif(Norig+Next)^.7*10
#' al=.15;bet=.5;
#' x1=bet*exp(al*t)*cos(t)+rnorm(length(t),0,.1)
#' y1=bet*exp(al*t)*sin(t)+rnorm(length(t),0,.1)
#'
#' D = as.matrix(dist(cbind(x1,y1)))
#' Dorig = D[1:Norig,1:Norig] # training distance matrix
#' DExt = D[(Norig+1):(Norig+Next),1:Norig] # new data distance matrix
#' # compute original diffusion map
#' dmap = diffuse(Dorig,neigen=2)
#'  # use Nystrom extension
#' dmapExt = nystrom(dmap,DExt)
#' plot(dmapExt[,1:2],pch=8,col=2,
#'   main="Diffusion map, black = original, red = new data",
#'   xlab="1st diffusion coefficient",ylab="2nd diffusion coefficient")
#' points(dmap$X[,1:2],pch=19,cex=.5)
nystrom <- function(dmap,Dnew,sigma=dmap$epsilon){

  Nnew = dim(Dnew)[1]
  Nold = dim(Dnew)[2]

  if(Nold != dim(dmap$X)[1]){
    stop("dimensions don't match")
  }

  Xnew = exp(-Dnew^2/(sigma))
  v = apply(Xnew, 1, sum)
  Xnew = Xnew/matrix(v,Nnew ,Nold)
  #nystrom extension:
  Xnew = Xnew %*% dmap$X %*% diag(1/dmap$eigenvals)

  return(Xnew)

}
