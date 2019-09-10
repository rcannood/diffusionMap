#' Compute diffusion map coordinates from pair-wise distances.
#'
#' Uses the pair-wise distance matrix for a data set to compute the diffusion
#' map coefficients.  Computes the Markov transition probability matrix, and
#' its eigenvalues and left & right eigenvectors.  Returns a 'dmap' object.
#'
#' Diffusion map is a powerful tool for data parametrization that exploits the
#' natural geometry of a data set.  Diffusion map uses local interactions
#' between data points, propogated to larger scales, to construct a global
#' representation of the data.
#'
#' The parameter eps.val controls the degree of localness in the diffusion
#' weight matrix.  For most statisitical inference problems using diffusion
#' map, results should be optimized over eps.val.  Generally a good starting
#' point is to pick eps.val as $2*$med.knn$^2$, where med.knn is the median
#' distance to the kth nearest neighbor, and k is chosen 1-2\% of n.  The
#' default uses 1\% of n.
#'
#' Computation of the diffusion map coordinates requires singular value
#' decomposition of the normalized graph Laplacian.  This operation is
#' optimized for speed by exploiting the sparseness of the graph Laplacian and
#' by using ARPACK for fast matrix decomposition.  Increasing the sparseness
#' parameter, delta, will speed up the algorithm.
#'
#' @param D n-by-n pairwise distance matrix for a data set with n points, or
#' alternatively output from the dist() function
#' @param eps.val epsilon parameter for the diffusion weight matrix,
#' exp(-D$^2$/(eps.val)).  Default is to use the epsilon corresponding to the
#' median distance to the 0.01*n nearest neighbor
#' @param neigen number of dimensions of final diffusion map representation.
#' Default uses number of dimensions corresponding to a 95\% drop-off in
#' eigenvalue multiplier.
#' @param t optional time-scale parameter in the diffusion map.  The
#' (recommended) default uses multiscale geometry.
#' @param maxdim the maximum number of diffusion map dimensions returned if
#' 95\% drop-off is not attained.
#' @param delta sparsity cut-off for the symmetric graph Laplacian.  Default of
#' 10^-5 is used.  Higher value induces more sparsity in Laplacian (and faster
#' computations)
#' @return The returned value is an object of 'class' 'diffuse'.
#'
#' The function 'plot' is used to plot the diffusion coordinates in 1, 2, or 3
#' dimensions.  The function 'print' displays the computed eigen-multipliers
#' and the value of epsilon used.
#'
#' An object of class 'dmap' is a list containing the following components:
#'
#' \item{X}{matrix of n diffusion map coordinates, entered column-wise (does
#' not include the trivial coordinate)} \item{phi0}{ trivial left eigenvector
#' of Markov matrix (stationary distribution of Markov random walk) in
#' diffusion map construction} \item{eigenvals}{eigen-values of the svd of the
#' symmetric graph Laplacian} \item{eigenmult}{eigen-multipliers of the
#' diffusion map} \item{psi}{right eigenvectors of the Markov matrix (first row
#' is the trivial right eigenvector)} \item{phi}{left eigenvectors of the
#' Markov matrix (first row is the trivial left eigenvector)}
#' \item{neigen}{number of diffusion map dimensions used} \item{epsilon}{the
#' value of epsilon used}
#' @references Coifman, R. R., \& Lafon, S., (2006), Appl. Comput. Harmon.
#' Anal., 21, 5
#'
#' Lafon, S., \& Lee, A., (2006), IEEE Trans. Pattern Anal. and Mach. Intel.,
#' 28, 1393
#'
#' Richards, J. W., Freeman, P. E., Lee, A. B., Schafer, C. M., (2009), ApJ,
#' 691, 32
#' @keywords multivariate nonparametric
#'
#' @export
#' @examples
#' library(stats)
#' ## example with noisy spiral
#' n=2000
#' t=runif(n)^.7*10
#' al=.15;bet=.5;
#' x1=bet*exp(al*t)*cos(t)+rnorm(n,0,.1)
#' y1=bet*exp(al*t)*sin(t)+rnorm(n,0,.1)
#' plot(x1,y1,pch=20,main="Noisy spiral")
#' D = dist(cbind(x1,y1))
#' dmap = diffuse(D,neigen=10) # compute diffusion map
#' par(mfrow=c(2,1))
#' plot(t,dmap$X[,1],pch=20,axes=FALSE,xlab="spiral parameter",ylab="1st diffusion coefficient")
#' box()
#' plot(1:10,dmap$eigenmult,typ='h',xlab="diffusion map dimension",ylab="eigen-multipliers")
#'
#' ## example with annulus data set
#' data(annulus)
#' plot(annulus,main="Annulus Data",pch=20,cex=.7)
#' D = dist(annulus) # use Euclidean distance
#' dmap = diffuse(D,eps.val=.1) # compute diffusion map & plot
#' print(dmap)
#' plot(dmap)
diffuse <- function(D, eps.val=epsilonCompute(D), neigen=NULL,
                    t=0, maxdim=50, delta=10^-5) {

  start = proc.time()[3]

  D = as.matrix(D)
  n=dim(D)[1]
  K = exp(-D^2/(eps.val))
  v=sqrt(apply(K,1,sum))
  A=K/(v%*%t(v));   # symmetric graph Laplacian

  # make A matrix sparse
  ind = which(A>delta, arr.ind=TRUE)
  Asp = sparseMatrix(i = ind[,1], j = ind[,2], x = A[ind], dims = c(n,n))

  f = function(x, A = NULL){ # matrix multiplication for ARPACK
    as.matrix(A %*% x)
  }

  cat('Performing eigendecomposition\n') # eigendecomposition
  if(is.null(neigen)){
    neff = min(maxdim+1,n)
  }else{
    neff =  min(neigen+1, n)
  }

  # eigendecomposition using ARPACK
  decomp = arpack(f,extra=Asp,sym=TRUE,
    options=list(which='LA',nev=neff,n=n,ncv=max(min(c(n,4*neff)))))
  psi = decomp$vectors/(decomp$vectors[,1]%*%matrix(1,1,neff))#right ev
  phi = decomp$vectors * (decomp$vectors[,1]%*%matrix(1,1,neff))#left ev
  eigenvals = decomp$values #eigenvalues


  cat('Computing Diffusion Coordinates\n')
  if(t<=0){# use multi-scale geometry
    lambda=eigenvals[-1]/(1-eigenvals[-1])
    lambda=rep(1,n)%*%t(lambda)
    if(is.null(neigen)){#use no. of dimensions corresponding to 95% dropoff
      lam = lambda[1,]/lambda[1,1]
      neigen = min(which(lam<.05)) # default number of eigenvalues
      neigen = min(neigen,maxdim)
      eigenvals = eigenvals[1:(neigen+1)]
      cat('Used default value:',neigen,'dimensions\n')
    }
    X = psi[,2:(neigen+1)]*lambda[,1:neigen] #diffusion coords. X
  }
  else{# use fixed scale t
    lambda=eigenvals[-1]^t
    lambda=rep(1,n)%*%t(lambda)

    if(is.null(neigen)){#use no. of dimensions corresponding to 95% dropoff
      lam = lambda[1,]/lambda[1,1]
      neigen = min(which(lam<.05)) # default number of eigenvalues
      neigen = min(neigen,maxdim)
      eigenvals = eigenvals[1:(neigen+1)]
      cat('Used default value:',neigen,'dimensions\n')
    }

    X = psi[,2:(neigen+1)]*lambda[,1:neigen] #diffusion coords. X
  }
  cat('Elapsed time:',signif(proc.time()[3]-start,digits=4),'seconds\n')

  y = list(X=X,phi0=phi[,1],eigenvals=eigenvals[-1],eigenmult=lambda[1,1:neigen],
    psi=psi,phi=phi,neigen=neigen,epsilon=eps.val)
  class(y) = "diffuse"
  return(y)
}


#' @export
#' @importFrom graphics plot stripchart
#' @importFrom scatterplot3d scatterplot3d
plot.diffuse <- function(x, y, ...){
  if(x$neigen>=3){
    scatterplot3d(x$X[,1],x$X[,2],x$X[,3],pch=20,cex.symbols=.6,xlab="1st Diffusion Coord.",ylab="2nd Diffusion Coord.",zlab="3rd Diffusion Coord.",main="3-Dimensional Diffusion Map")
  }else if(x$neigen==2){
    plot(x$X[,1],x$X[,2],pch=20,cex=.6,xlab="1st Diffusion Coord.",ylab="2nd Diffusion Coord.",main="2-Dimensional Diffusion Map")
  }else {
    stripchart(x$X,method='jitter',pch=20,cex=.6,xlab="1st Diffusion Coord.",main="1-Dimensional Diffusion Map")
  }
}
#' @method print diffuse
#' @export
print.diffuse <- function(x, ...) {
  lis <- list(
    eigenmult=x$eigenmult,
    epsilon=x$epsilon
  )
  print(
    lis,
    digits = NULL,
    quote = TRUE,
    na.print = NULL,
    print.gap = NULL,
    right = FALSE,
    max = NULL,
    ...
  )
}
