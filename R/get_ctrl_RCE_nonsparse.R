get_ctrl_RCE_nonsparse <- function(samps,ctrl_ind){
  
  Yi = ncol(samps$data)
  
  ind = to_indicators(samps$batches)
  
  psi1 = CV_KRR(cbind(samps$controls[,ctrl_ind],ind),samps$data[,Yi],maxY=max(samps$data[,Yi]),minY=min(samps$data[,Yi]))
  
  if (ncol(ind)==0){
    psi2 = mean(samps$data[,Yi])
  } else{
    psi2 = CV_KRR(ind,samps$data[,Yi],maxY=max(samps$data[,Yi]),minY=min(samps$data[,Yi]))
  }
  
  return(psi1 - psi2)
  
}