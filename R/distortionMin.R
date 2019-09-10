#' Distortion Minimization via K-means
#'
#' Runs one K-means loop based on the diffusion coordinates of a data set,
#' beginning from an initial set of cluster centers.
#'
#' Used by diffusionKmeans().
#'
#' @param X diffusion coordinates, each row corresponds to a data point
#' @param phi0 trivial left eigenvector of Markov matrix (stationary
#' distribution of Markov random walk) in diffusion map construction
#' @param K number of clusters
#' @param c0 initial cluster centers
#' @param epsilon stopping criterion for relative change in distortion
#' @return The returned value is a list with components
#'
#' \item{S}{ labelling from K-means loop. n-dimensional vector with integers
#' between 1 and K} \item{c}{ K geometric centroids found by K-means}
#' \item{D}{minimum of total distortion (loss function of K-means) found in
#' K-means run} \item{DK}{n by k matrix of squared (Euclidean) distances from
#' each point to every centroid}
#' @seealso [diffusionKmeans()]
#' @references Lafon, S., \& Lee, A., (2006), IEEE Trans. Pattern Anal. and
#' Mach. Intel., 28, 1393
#' @keywords multivariate nonparametric
#' @export
#' @examples
#'
#' data(annulus)
#' n = dim(annulus)[1]
#' D = dist(annulus) # use Euclidean distance
#' dmap = diffuse(D,0.03) # compute diffusion map
#' km = distortionMin(dmap$X,dmap$phi0,2,dmap$X[sample(n,2),])
#' plot(annulus,col=km$S,pch=20)
#' table(km$S,c(rep(1,500),rep(2,500)))
distortionMin <- function(X,phi0,K,c0,epsilon=0.001){

  n=dim(X)[1]
  c=c0
  oldD=Inf

  MaxIter=1000

  for(ii in 1:MaxIter){ #K-means loop
    DX=c()
    for(jj in 1:K){
      dX = X-matrix(1,n,1)%*%c[jj,] # n by p
      DX=cbind(DX,apply(dX*dX,1,sum))# n by K
    }
    S = apply(DX,1,which.min) # new labels
    Dtmp = apply(DX,1,min) # new dists. to centroids

    for(jj in 1:K){ # check for empty clusters
      ind=which(S==jj)
      if(length(ind)==0){# if cluster jj empty
        S[which.max(Dtmp)]=jj # give it the pt. furthest from its centroid
        Dtmp[which.max(Dtmp)]=0 # and set its distance to 0
      }
    }

    for(jj in 1:K){# update diffusion centroids
      ind=which(S==jj)
      c[jj,]=t(phi0[ind])%*%X[ind,]/sum(phi0[ind]) #centroid of cluster jj
    }
    D = Dtmp%*%phi0 # distortion

    if((oldD-D)/D < epsilon){ #stopping criterion
      break
    }
    oldD=D
  }
  if(ii==MaxIter) print('Maximum # of iterations reached')

  return(list(S=S,c=c,D=D,DX=DX))
}
