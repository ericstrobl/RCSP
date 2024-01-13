save_samps_by_file_mult3 <- function(DAG2,nsamps=200){
  
  cnt = -1
  sampsA = list(); sampsA$data = c(); sampsA$dataT = c(); sampsA$batches = c(); sampsA$interv = c()
  
  sc_batches = 1:DAG2$DAGs$n_batch
  
  for (j in 0:(ncol(DAG2$DAGs$graph)-1)){
    cnt = cnt + 1
    # print(j)
    DAGs = DAG2$DAGs
    if (j > 0){
      DAGs$off[j] = -2  # knockdown
    }
    samps = sample_DAG_NB_linear(nsamps,DAGs) #######
    
    sampsA$data = rbind(sampsA$data, samps$data[1:nsamps,])
    sampsA$batches = c(sampsA$batches, samps$batches[1:nsamps])
    sampsA$interv = c(sampsA$interv, rep(j,nsamps))
    
    dataT = samps$dataT
    samps$dataT = c(); samps$data = samps$data[1:nsamps,]; samps$batches = samps$batches[1:nsamps]
    save(file=paste("samps_interv",as.character(j),".RData",sep=""),samps) #data, batches
    
    samps$dataT = dataT
    samps$data = c(); samps$batches = c()
    save(file=paste("samps_intervT",as.character(j),".RData",sep=""),samps) #dataT
    
    
    if ( (cnt == 250) | j== (ncol(DAG2$DAGs$graph)-1) ){
      
      for (b in sc_batches){
        ix = which(sampsA$batches==b)
        if (j > cnt ){
          load(file=paste("samps_batch",as.character(b),".RData",sep=""))
        } else{
          samps$data = c(); samps$batches = c(); samps$dataT = c();
          samps$interv = c(); samps$ix = c()
        }
        samps$data = rbind(samps$data, sampsA$data[ix,])
        samps$batches = c(samps$batches, sampsA$batches[ix])
        samps$interv = c(samps$interv, sampsA$interv[ix])
        
        for (jj in 0:(ncol(DAG2$DAGs$graph)-1)){
          samps$ix = c(samps$ix, which(sampsA$batches[sampsA$interv==jj]==b))
        }
        
        save(file=paste("samps_batch",as.character(b),".RData",sep=""),samps)
      }
      
      cnt = 0
      sampsA = list(); sampsA$data = c(); sampsA$dataT = c(); sampsA$batches = c(); sampsA$interv = c()
    }
    
    
  }
}
