#' Diffusion K-means
#'
#' Clusters a data set based on its diffusion coordinates.
#'
#' A '"dmap"' object computed by diffuse() is the input, so diffuse() must be
#' performed first.  Function is written this way so the K-means parameters may
#' be varied without having to recompute the diffusion map coordinates in each
#' run.
#'
#' Diffusion K-means is a special form of spectral clustering.  It is a unique
#' algorithm because the eigenvectors of the symmetric Laplacian are weighted
#' in such a way to guarantee that Euclidean distance in diffusion space will
#' be approximately equal to the diffusion distance between objects.
#' Clustering by Euclidean distance in diffusion space exploits this fact.
#'
#' @param dmap a '"dmap"' object, computed by diffuse()
#' @param K number of clusters
#' @param params optional parameters for each data point.  Entry can be a
#' vector of length n, or a matrix with n rows. If this argument is given,
#' cluster centroid parameters are returned.
#' @param Niter number of K-means iterations performed.
#' @param epsilon stopping criterion for relative change in distortion for each
#' K-means iteration
#' @return The returned value is a list with components
#'
#' \item{part}{final labelling of data from K-means. n-dimensional vector with
#' integers between 1 and K} \item{cent}{ K geometric centroids found by
#' K-means} \item{D}{minimum of total distortion (loss function of K-means)
#' found across K-means runs} \item{DK}{n by k matrix of squared (Euclidean)
#' distances from each point to every centroid for the optimal K-means run}
#' \item{centparams}{optional parameters for each centroid.  Only returned if
#' params is specified in the function call.  Is a matrix with k rows.}
#' @seealso [diffuse()]
#' @references Lafon, S., \& Lee, A., (2006), IEEE Trans. Pattern Anal. and
#' Mach. Intel., 28, 1393
#'
#' Richards, J. W., Freeman, P. E., Lee, A. B., and Schafer, C. M., (2009),
#' ApJ, 691, 32
#'
#' Richards, J. W., Freeman, P. E., Lee, A. B., Schafer, C. M., (2009), MNRAS,
#' Volume 399, Issue 2, pp. 1044-1057
#' @keywords multivariate nonparametric
#'
#' @export
#' @examples
#' library(scatterplot3d)
#'
#' ## example with annulus data set
#' data(annulus)
#' par(mfrow=c(2,1))
#' plot(annulus,main="Annulus Data",pch=20,cex=.7)
#' D = dist(annulus) # use Euclidean distance
#' dmap = diffuse(D,eps.val=0.05) # compute diffusion map
#' k=2  # number of clusters
#' dkmeans = diffusionKmeans(dmap, k)
#' plot(annulus,main="Colored by diffusion K-means clustering",pch=20,
#'    cex=.7,col=dkmeans$part)
#' table(dkmeans$part,c(rep(1,500),rep(2,500)))
#'
#'
#' ## example with Chainlink data set
#' data(Chainlink)
#' lab.col = c(rep("red",500),rep("blue",500)); n=1000
#' scatterplot3d(Chainlink$C1,Chainlink$C2,Chainlink$C3,color=lab.col,
#'    main="Chainlink Data") # plot Chainlink data
#' D = dist(Chainlink) # use Euclidean distance
#' dmap = diffuse(D,neigen=3,eps.val=.01) # compute diffusion map & plot
#' plot(dmap)
#' dkmeans = diffusionKmeans(dmap, K=2)
#' col.dkmeans=ifelse(dkmeans$part==1,"red","blue")
#' scatterplot3d(Chainlink,color=col.dkmeans,
#'    main="Chainlink Data, colored by diff. K-means class")
#' table(dkmeans$part,lab.col)
diffusionKmeans <- function(dmap, K, params=c(), Niter=10, epsilon=0.001){

  n=dim(dmap$X)[1]
  D=Inf # max. distortion

  for(ii in 1:Niter){# run Niter K-means loops
    print(paste('Iteration',ii,'of',Niter))
    c0 = dmap$X[sample(1:n,K),] # choose random initial centroids
    kmeans = distortionMin(dmap$X,dmap$phi0,K,c0,epsilon) # K-means loop
    if(kmeans$D<D){ # keep best result
      D = kmeans$D
      DX = kmeans$DX
      part = kmeans$S
      cent = kmeans$c
    }
  }

  if(!is.null(params)){# if parameters are given for each data point, compute the centroid parameters
    npar = dim(params)[2]
    npar = ifelse(is.null(npar),1,npar)
    params = matrix(params,n,npar)
    centparams=matrix(0,K,npar)
    for(jj in 1:K){
      ind=which(part==jj)
      centparams[jj,] = t(dmap$phi0[ind])%*%params[ind,]/sum(dmap$phi0[ind])
    }
    return(list(part=part,cent=cent,D=D,DX=DX,centparams=centparams))
  }
  else{
    return(list(part=part,cent=cent,D=D,DX=DX))
  }
}
