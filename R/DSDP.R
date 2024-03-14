DSDP <- function(bulk_samps,desL=NULL,stable=NULL,reg="KRR",verbose=FALSE){# needs bulk samples
  # Similar to the RCSP algorithm but computes Deviation of Statistical Dependence (D-SD) by removing the conditioning on the surrogate ancestors
  # 
  # bulk_samps denotes a list of bulk RNAseq samples, where:
  #     (1) bulk_samps$data contains RNA counts
  #     (2) bulk_samps$controls contains control varaibles (e.g., age, gender)
  #     (3) bulk_samps$batches contains batch indicators
  #
  # reg is set to either "KRR" for kernel ridge regression, or "MARS" for multivariate adaptive regression splines
  #     KRR is good for dense output, MARS is good for speed
  #
  # desL corresponds to the list of descendants inferred from Perturb-seq. This is computed inside the algorithm by default.
  
  require(earth)#
  
  if (is.null(desL)){
    desL = find_des_final(bulk_samps)
  }
  
  interv = desL$interv
  desL = desL$desL
  
  ind = to_indicators(bulk_samps$batches,unique(bulk_samps$batches))#
  data = bulk_samps$data####
  batches = bulk_samps$batches
  
  p = ncol(data)#
  Yi = p
  
  REs = matrix(0,nrow(bulk_samps$data),length(interv))#
  
  maxY = max(data[,Yi])
  minY = min(data[,Yi])
  
  for (i in seq_len(length(interv)) ){#
     if (verbose) print(i)
    
    pa = c()
    for (k in setdiff(interv,desL[[i]])){# variables that are not descendants of i or stable
      if (interv[i] %in% desL[[which(interv==k)]]){ # variables that are ancestors of i
        pa = c(pa, k)#
      }
    }
    
    if (length(stable)==0){
      L = c()
    } else{
      L = rowSums(data)
    }
    
    if (reg == "KRR"){
      psi1_1 = CV_KRR(cbind(data[,c(pa,interv[i])],L,ind,bulk_samps$controls),data[,Yi],maxY,minY)
      psi1_2 = CV_KRR(cbind(data[,pa],L,ind,bulk_samps$controls),data[,Yi],maxY,minY)
    } else if (reg == "KRR"){
      psi1_1 = earth(cbind(data[,c(pa,interv[i])],L,ind,bulk_samps$controls),data[,Yi])$fitted.values
      psi1_2 = earth(cbind(data[,pa],L,ind,bulk_samps$controls),data[,Yi])$fitted.values
    }
    
    REs[,i] = psi1_1 - psi1_2#
    
  }
  
  ix = which(colMeans(abs(REs))==0)
  if (length(ix)>0){
    REs = REs[,-ix]
    interv = interv[-ix]
  }
  genes = colnames(bulk_samps$data)[interv]
  
  return(list(REs = as.matrix(REs), interv=interv, genes = genes))
  
}
