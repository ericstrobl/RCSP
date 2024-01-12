compute_RMSE3 <- function(tr,out,p,err=F){
  
  MSE_REs_out = matrix(0,nrow(tr$REs),p-1)
  MSE_REs_out[,out$interv] = abs(out$REs)
  MSE_REs_tr = matrix(0,nrow(tr$REs),p-1)
  if (err){
    MSE_REs_tr[,tr$interv] = abs(tr$tREs)
  } else{
    MSE_REs_tr[,tr$interv] = abs(tr$REs) 
  }
  
  RMSE = sqrt(mean(  (MSE_REs_out - MSE_REs_tr)^2 ))
  
  return(list( RMSE_REs = RMSE ))
  
}