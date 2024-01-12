reg_norm <- function(samps,desL,stable){# needs bulk samples
  
  ## multivariate regression residuals
  
  require(earth)#
  
  # desL = find_des2(samps,ctrl_ind)
  interv = desL$interv
  desL = desL$desL
  
  ind = to_indicators(samps$batches,unique(samps$batches))#
  data = samps$data####
  batches = samps$batches
  
  p = ncol(data)#
  Yi = p
  
  REs = matrix(0,nrow(samps$data),length(interv))#
  
  
  for (i in seq_len(length(interv)) ){#
    
    REs[,i] = earth(cbind(data[,-c(interv[i],Yi)],ind),data[,Yi])$fitted.values - data[,Yi]
    
  }
  
  ix = which(colMeans(abs(REs))==0)
  if (length(ix)>0){
    REs = REs[,-ix]
    interv = interv[-ix]
  }
  genes = colnames(samps$data)[interv]
  
  return(list(REs = as.matrix(REs), interv=interv, genes = genes))
  
}