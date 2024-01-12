CV_KRR <- function(X,Y,maxY=NULL,minY=NULL){
  require(Rfast)
  Y=as.matrix(Y)
  d=ncol(Y)
  mY = colMeans(Y)
  sY = apply(Y,2,sd)
  for (i in 1:d){
    Y[,i] = (Y[,i,drop=FALSE] - mY[i])/sY[i]
  }
  
  
  X=normalizeData(X)
  
  n = nrow(Y);
  # shuf = sample(1:n,n,replace=FALSE)
  # X[1:n,]=X[shuf,]
  # Y[1:n,] = Y[shuf,]
  
  dotX = as.matrix(dist(X))
  med_dist = median(dotX);
  if (med_dist == 0){ med_dist = 1;}
  # med_dist = 1
  dotX = dotX^2;
  
  lambdas = c(1,1E-1,1E-2,1E-3,1E-4,1E-5,1E-6);
  sigmas = med_dist*seq(0.5,4,0.25);
  
  ## estimate conditional expectations
  err = Inf; #maximum error allowed
  
  for (s in seq_len(length(sigmas))){
    K = RBF_kernel(dotX,sigmas[s])
    # Ke = eigen(K,symmetric=TRUE)
    for (l in seq_len(length(lambdas))){
      # K_inv = Ke$vector %*% diag(1/(Ke$values+n*lambdas[l])) %*% t(Ke$vector)
      K_inv = spdinv( K  + n*diag(n)*lambdas[l]);
      H = K_inv %*% K;
      te_y = t(t(Y) %*% H)
      res = LOO_KRR(Y, Y - te_y, H)
      if (res$err<err){
        err = res$err
        pre = res$pre
      }
    }
  }
  
  if (!is.null(maxY)){
    pre = pmax(pmin(pre*sY[i] + mY[i],maxY),minY)
  } else{
    pre = pre*sY[i] + mY[i]
  }
  
  return(pre)
  
}