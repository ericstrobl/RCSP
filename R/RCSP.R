RCSP <- function(bulk_samps,reg="KRR"){# needs bulk samples
  #
  # bulk_samps denotes a list of bulk RNAseq samples, where:
  #     (1) bulk_samps$data contains RNA counts
  #     (2) bulk_samps$controls contains control varaibles (e.g., age, gender)
  #     (3) bulk_samps$batches contains batch indicators
  #
  # reg is set to either "KRR" for kernel ridge regRCSsion, or "MARS" for multivariate adaptive regRCSsion splines
  #     KRR is good for dense output, MARS is good for speed
  #
  
  require(earth)#
  
  desL = find_des_final(bulk_samps)
  
  interv = desL$interv
  desL = desL$desL
  
  ind = to_indicators(bulk_samps$batches,unique(bulk_samps$batches))#
  data = bulk_samps$data####
  batches = bulk_samps$batches
  
  p = ncol(data)#
  Yi = p
  
  RCS = matrix(0,nrow(bulk_samps$data),length(interv))#
  
  maxY = max(data[,Yi])
  minY = min(data[,Yi])
  
  for (i in seq_len(length(interv)) ){#
    # print(i)
    
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
    
    RCS[,i] = psi1_1 - psi1_2#
    
  }
  
  ix = which(colMeans(abs(RCS))==0)
  if (length(ix)>0){
    RCS = RCS[,-ix]
    interv = interv[-ix]
  }
  genes = colnames(bulk_samps$data)[interv]
  
  return(list(RCS = as.matrix(RCS), interv=interv, genes = genes))
  
}