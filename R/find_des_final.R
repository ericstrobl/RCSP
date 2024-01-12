find_des_final <- function(samps,alpha=0.01){
  
  p = ncol(samps$data)-1 # exclude target
  interv = 1:p # list of variables with interventions/perturbations
  
  ### use Perturb-seq to find descendants
  
  load("samps_interv0.RData") # control samples
  data0 = samps$data
  
  desL = vector("list",length(interv)) # list of descendants of each variable
  for (i in 1:length(interv)){
    load(paste("samps_interv",as.character(interv[i]),".RData",sep="")) # data with X_i perturbed
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
