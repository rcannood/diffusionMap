diffuse <- function(D,eps.val=epsilonCompute(D),neigen=NULL,t=0,maxdim=50,delta=10^-5) {

  start = proc.time()[3]

  D = as.matrix(D)
  n=dim(D)[1]
  K = exp(-D^2/(eps.val)) 
  v=sqrt(apply(K,1,sum))
  A=K/(v%*%t(v));   # symmetric graph Laplacian

  # make A matrix sparse
  ind = which(A>delta, arr.ind=TRUE)
  Asp = sparseMatrix(i = ind[,1], j = ind[,2],x=A[ind], dims=c(n,n))

  f = function(x,A=NULL){ # matrix multiplication for ARPACK
    as.matrix(A%*%x)
  }

  print('Performing eigendecomposition') # eigendecomposition
  if(is.null(neigen)){ 
    neff = min(maxdim+1,n)  
  }else{
    neff =  min(neigen+1, n)
  }
  
  # eigendecomposition using ARPACK
  decomp = arpack(f,extra=Asp,sym=TRUE,options=list(which='LA',nev=neff,n=n,ncv=max(30,2*neff)))
  psi = decomp$vectors[,neff:1]/(decomp$vectors[,neff]%*%matrix(1,1,neff))#right ev
  phi = decomp$vectors[,neff:1] * (decomp$vectors[,neff]%*%matrix(1,1,neff))#left ev
  eigenvals = decomp$values[neff:1]#eigenvalues


  print('Computing Diffusion Coordinates') # Compute Diffusion Coordinates
  if(t<=0){# use multi-scale geometry
    lambda=eigenvals[-1]/(1-eigenvals[-1])
    lambda=rep(1,n)%*%t(lambda)
    if(is.null(neigen)){#use no. of dimensions corresponding to 95% dropoff
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
    if(is.null(neigen)){#use no. of dimensions corresponding to 95% dropoff
      lam = lambda[1,]/lambda[1,1]
      neigen = min(which(lam<.05)) # default number of eigenvalues
      neigen = min(neigen,maxdim)
      eigenvals = eigenvals[1:(neigen+1)]  
      print(paste('Used default value:',neigen,'dimensions'))
    }
    X = psi[,2:(neigen+1)]*lambda[,1:neigen] #diffusion coords. X
  }
  print(paste('Elapsed time:',signif(proc.time()[3]-start,digits=4),'seconds'))

  y = list(X=X,phi0=phi[,1],eigenvals=eigenvals[-1],eigenmult=lambda[1,1:neigen],psi=psi,phi=phi,neigen=neigen,epsilon=eps.val)
  class(y) = "dmap"
  return(y)

}


