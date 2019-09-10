#' Compute default diffusion map epsilon.
#'
#' Uses the pair-wise distances to estimate a diffusion map epsilon value by
#' the median p*n-th nearest neighbor
#'
#' Function is used as the default value in diffuse().  For inference problems,
#' it is advised that the results be optimized over epsilon.
#'
#' @param D n-by-n pairwise distance matrix for a data set with n points, or
#' alternatively output from the dist() function
#' @param p distances to p*n-th nearest neighbor are used.  Default value is
#' .01
#' @return \item{epsilon }{value of epsilon to be used in diffusion map}
#' @seealso [diffuse()]
#' @keywords multivariate nonparametric
#' @export
#' @importFrom stats median
#' @examples
#' data(annulus)
#' D = dist(annulus) # use Euclidean distance
#' epsilonCompute(D,.005)
#' epsilonCompute(D,.01)
#' epsilonCompute(D,.05)
#' epsilonCompute(D,.1)
epsilonCompute <- function(D,p=.01){

  D = as.matrix(D)
  n = dim(D)[1]
  k = ceiling(p*n)
  k = ifelse(k<2,2,k) # use k of at least 2
  D.sort = apply(D,1,sort)
  dist.knn = D.sort[(k+1),] # find dists. to kth nearest neighbor
  epsilon = 2*median(dist.knn)^2

  return(epsilon)
}

