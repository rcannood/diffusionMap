diffuse <- function(D,eps.val=epsilonCompute(D),neigen=0,t=0,maxdim=50) {

  D = as.matrix(D)
  n=dim(D)[1]
  K = exp(-D^2/(4*eps.val)) #or equivalently, K=exp(-(D/sigmaK)^2) 
  v=sqrt(apply(K,1,sum))
  A=K/(v%*%t(v));   # symmetric graph Laplacian

  print('Performing eigendecomposition')#Eigendecomposition
  if(neigen<=0){ # eigendecomposition of symmetric matrix
    decomp = svd(A, nu = min(maxdim+1, n), nv = min(maxdim+1, n))
    psi = decomp$u/(decomp$u[,1]%*%matrix(1,1,min(maxdim+1,n)))#right ev
    phi = decomp$u * (decomp$u[,1]%*%matrix(1,1,min(maxdim+1,n)))#left ev
    eigenvals = decomp$d[1:min(maxdim+1, n)]#eigenvalues
  }else{
    decomp = svd(A, nu = min(neigen+1, n), nv = min(neigen+1, n))
    psi = decomp$u/(decomp$u[,1]%*%matrix(1,1,min(neigen+1,n)))#right ev
    phi = decomp$u * (decomp$u[,1]%*%matrix(1,1,min(neigen+1,n)))#left ev
    eigenvals = decomp$d[1:min(neigen+1, n)]#eigenvalues
  }  

  print('Computing Diffusion Coordinates') # Compute Diffusion Coordinates
  if(t<=0){# use multi-scale geometry
    lambda=eigenvals[-1]/(1-eigenvals[-1])
    lambda=rep(1,n)%*%t(lambda)
    if(neigen<=0){#use no. of dimensions corresponding to 95% dropoff
      lam = lambda[1,]/lambda[1,1]
      neigen = min(which(lam<.05)) # default number of eigenvalues
      neigen = min(neigen,maxdim)
      eigenvals = eigenvals[1:(neigen+1)]  
      print(paste('Used default value:',neigen,'dimensions'))
    }
    X = psi[,2:(neigen+1)]*lambda[,1:neigen] #diffusion coords. X
  }
  else{# use fixed scale t
    lambda=eigenvals[-1]^t
    lambda=rep(1,n)%*%t(lambda)
    if(neigen<=0){#use no. of dimensions corresponding to 95% dropoff
      lam = lambda[1,]/lambda[1,1]
      neigen = min(which(lam<.05)) # default number of eigenvalues
      neigen = min(neigen,maxdim)
      eigenvals = eigenvals[1:(neigen+1)]  
      print(paste('Used default value:',neigen,'dimensions'))
    }
    X = psi[,2:(neigen+1)]*lambda[,1:neigen] #diffusion coords. X
  }

  y = list(X=X,phi0=phi[,1],eigenmult=lambda[1,1:neigen],psi=psi,phi=phi,neigen=neigen,epsilon=eps.val)
  class(y) = "dmap"
  return(y)

}
