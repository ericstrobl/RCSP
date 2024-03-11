find_des_AMD_final <- function(samps,ctrl_ind = NULL,alpha=0.01){
  
  ### find variables that are d-connected to Y in BULK
  load("samps_bulk_AMD.RData")#
  p = ncol(samps$data)-1 # exclude target
  
  interv = 1:p

  ### use single cell interventions to find descendants
  
  load("samps_interv0_AMD.RData")#
  data0 = samps$data
  
  desL = vector("list",length(interv))
  for (i in 1:length(interv)){
    print(i)
    load(paste("samps_interv",as.character(interv[i]),"_AMD.RData",sep=""))
    desL[[i]] = interv[i]
    for (j in setdiff(interv,interv[i])){
      if( (mean(data0[,j]) != 0) | (mean(samps$data[,j]) != 0) ){ # check if t.test can be run
        if (t.test(data0[,j],samps$data[,j])$p.value < alpha){
          desL[[i]] = c(desL[[i]],j)
        }
      }
    }
  }
  
  return(list(desL=desL,interv=interv))
  
}