sample_DAG_NB_linear <- function(nsamps,DAG){
  
  require(MASS)
  
  G = DAG$graph
  p = nrow(G)
  
  dataT = matrix(0,nsamps,p)
  dataE = dataT
  off = DAG$off
  
  L = length(DAG$ord)
  # for (i in setdiff(1:L,DAG$HK)){
  #   dataE[,i] = 0.5*rnorm(nsamps)
  # }
  # dataT = dataE %*% ginv(diag(L)-as.matrix(DAG$weights))
  for (i in 1:L){
    c = DAG$ord[i]
    h = which(G[DAG$ord[1:(i-1)],c]==1)
    h = DAG$ord[1:(i-1)][h]

    if (length(h)>0){
      A = dataT[,h,drop=FALSE]%*%as.matrix(DAG$weights[h,c,drop=FALSE]) # do not need to include offset because the gamma error terms introduce the offset
    } else{
      A = rep(0,nsamps)
    }

    if (c == DAG$Y){
      dataT[,c]= A + 0.5*rnorm(nsamps)
    } else{
      if (c %in% DAG$HK){
        dataT[,c]= off[c]
      } else{
        dataE[,c] = 0.5*rnorm(nsamps)
        dataT[,c]= A + dataE[,c]
      }
      # dataT[,c] = rgamma(nsamps,shape = r[c], rate = r[c]/softplus(A + off[c]))
    }

  }
  
  dataT[,1:(L-1)] = softplus(dataT[,1:(L-1)] + off[-L])
  
  # batch modification
  batches = sample(DAG$n_batch,nsamps,replace=TRUE) # each sample is from the same batch
  datap = dataT
  for (b in 1:DAG$n_batch){
    ib = which(batches == b)
    datap[ib,-DAG$Y] = t(t(dataT[ib,-DAG$Y])* DAG$batch[b,-DAG$Y])  #####
  }
  # # introduce Poisson measurement error from multinomial
  total = rowSums(datap[,-DAG$Y])###
  
  datap[,-DAG$Y] = datap[,-DAG$Y]/total # compute proportions ###
  
  N = rep(0,nsamps)
  data = dataT
  for (n in 1:nsamps){ # sample from multinomial with proportions
    N[n] = rpois(1,total[n])
    data[n,-DAG$Y] = rmultinom(1,N[n],prob=datap[n,-DAG$Y])
  } # measurement error marginal distributions then follow a Poisson
  
  
  return(list(data=data,dataT=dataT,dataE=dataE,batches=batches))
}

