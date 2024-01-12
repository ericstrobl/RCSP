lin_norm <- function(samps,desL,stable){# needs bulk samples
  
  ## LiNGAM
  
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
  resid = REs
  
  for (i in seq_len(length(interv)) ){#
    # print(i)
    
    pa = c()
    for (k in setdiff(interv,c(desL[[i]],stable))){# variables that are not descendants of i or stable
      if (interv[i] %in% desL[[which(interv==k)]]){ # variables that are ancestors of i
        pa = c(pa, k)#
      }
    }
    
    # pa = which(DAGs$graph[,interv[i]]!=0) ###
    pa = setdiff(pa,stable) ###
    # print(c(interv[i],pa))
    
    if (length(pa)>0){
      resid[,i]  = data[,interv[i]]-lm.fit(cbind(data[,pa],1),data[,interv[i]])$fitted.values
    } else{
      resid[,i]  = data[,interv[i]]-mean(data[,interv[i]])
    }
    
  }
  
  # if (!is.null(stable)){
  #   L = rowSums(data[,stable,drop=FALSE]); #L = rep(1,nrow(data)) ###
  # } else{
  #   L = rowSums(data)
  # }
  
  psi = earth(cbind(samps$data[,-Yi],ind),data[,Yi])$fitted.values
  for (i in seq_len(length(interv)) ){#
    REs[,i] = psi - earth( resid[,-i],data[,Yi])$fitted.values
  }
  
  # REs[,i] = earth( resid[,i],data[,Yi] )$fitted.values - mean(data[,Yi])
  
  ix = which(colMeans(abs(REs))==0)
  if (length(ix)>0){
    REs = REs[,-ix]
    interv = interv[-ix]
  }
  genes = colnames(samps$data)[interv]
  
  return(list(REs = as.matrix(REs), interv=interv, genes = genes))
  
}