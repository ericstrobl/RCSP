cor_norm <- function(samps,desL,stable){#
  
  ## univariate regression residuals
  
  require(earth)#
  
  load("samps_bulk_synth.RData")#
  
  interv = desL$interv
  desL = desL$desL
  
  Yi = ncol(samps$data)
  
  ind = to_indicators(samps$batches,unique(samps$batches))#
  
  REs = matrix(0,nrow(samps$data),length(interv))#
  
  for (i in seq_len(length(interv))){#
    # 
    # pa = c()
    # for (k in setdiff(interv,c(desL[[i]],stable))){# variables that are not descendants of i or stable
    #   if (interv[i] %in% desL[[which(interv==k)]]){ # variables that are ancestors of i
    #     pa = c(pa, k)#
    #   }
    # }
    # 
    # # pa = which(DAGs$graph[,interv[i]]!=0) ###
    # pa = setdiff(pa,stable) ###
    
    if (!is.null(stable)){
      L = rowSums(samps$data[,setdiff(stable,interv[i]),drop=FALSE])
    } else{
      L = rowSums(samps$data)
    }

    psi = earth(cbind(samps$data[,interv[i]],L,ind),samps$data[,Yi])$fitted.values#
    
    REs[,i] = psi  - samps$data[,Yi]
  }
  
  return(list(REs=REs,interv=interv))#
  
}