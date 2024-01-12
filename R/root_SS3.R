root_SS3 <- function(tr,out,p,err=F){
  
  require(fastcluster)
  require(pdfCluster)
 
  cl_RE <- hclust.vector(abs(out$REs), method="ward")
  MSE_REs_out = matrix(0,nrow(tr$REs),p-1)
  MSE_REs_out[,out$interv] = abs(out$REs)
  MSE_REs_tr = matrix(0,nrow(tr$REs),p-1)
  
  if (err){
    cl_tRE <- hclust.vector(abs(tr$tREs), method="ward")
    MSE_REs_tr[,tr$interv] = tr$tREs
  } else{
    cl_tRE <- hclust.vector(abs(tr$REs), method="ward")
    MSE_REs_tr[,tr$interv] = tr$REs
  }
  
  RMSE = rep(0,10)
  for (c in 1:10){
    ix_RE = cutree(cl_RE, k = c)
    ix_tRE = cutree(cl_tRE, k = c)
    
    for (cc in 1:c){
      RMSE[c] = RMSE[c] + sqrt(mean(  sweep(MSE_REs_out[ix_RE,],2,colMeans(MSE_REs_tr[ix_tRE,]),"-")^2 ))
    }
  }
  
  return( list(RMSE_RE = RMSE/c(1:10))  )
  
}