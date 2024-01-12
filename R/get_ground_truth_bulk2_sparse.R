get_ground_truth_bulk2_sparse <- function(DAGs,interv,stable){
  require(earth)
  
  load("samps_bulk_synth.RData")
  # plot(as(as.matrix(DAGs$graph),"graphNEL"))
  p = nrow(DAGs$graph)
  
  dataT0 = samps$dataT#
  n = nrow(dataT0)
  
  AncY = isAncAll(DAGs$graph,1:p,DAGs$Y) # ancestors of Y and intervened on
  
  intervf = c()
  for (i in interv){
    if(isDec(DAGs$graph,DAGs$Y,i)){
      intervf = c(intervf,i)
    } 
  }
  interv = intervf
  
  REs = matrix(0,nrow(samps$dataT),length(interv))#
  # REs_test = REs
  
  weightsT = DAGs$weights[c(interv,p),c(interv,p)]
  total = solve(diag(length(interv)+1) - weightsT)
  
  tREs = matrix(0,nrow(samps$dataT),length(interv))#
  for (i in seq_len(length(interv)) ){
    tREs[,i] = samps$dataE[,interv[i]] * total[i,length(interv)+1]
  }
  
  ind = to_indicators(samps$batches,unique(samps$batches))#
  
  for (i in seq_len(length(interv)) ){
    
    if( !(interv[i] %in% AncY) ) next
    
    pa = which(DAGs$graph[,interv[i]]!=0)
    
    # print(c(interv[i],pa))
    psi1_1 = earth(dataT0[,c(pa,interv[i])],dataT0[,DAGs$Y])$fitted.values
    if (length(pa)>0){
      psi1_2 = earth(dataT0[,pa],dataT0[,DAGs$Y])$fitted.values
    } else{
      psi1_2 = mean(dataT0[,DAGs$Y])
    }
    
    REs[,i] = psi1_1 - psi1_2#
    
    ####
    L = rowSums(samps$data[,DAGs$HK,drop=FALSE])
    
    # psi1_1 = earth( cbind(samps$data[,c(pa,interv[i])],L,ind),samps$data[,DAGs$Y])$fitted.values
    # psi1_2 = earth( cbind(samps$data[,pa],L,ind),samps$data[,DAGs$Y])$fitted.values
    # REs_test[,i] = psi1_1 - psi1_2#
    
    
  }
  
  return(list(REs = as.matrix(REs), tREs = as.matrix(tREs), interv=interv))
  
}

