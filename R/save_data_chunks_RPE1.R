save_data_chunks_RPE1 <- function(gene_symbols,RPE1_directory){
  # gene symbols = gene symbols of bulk RNA-seq data
  # RPE1 directory = directory of ReplogleWeissman2022_rpe1.h5ad, e.g., RPE1_directory = "C:/Users/ericv/Documents/ReplogleWeissman2022_rpe1.h5ad"
  
  require(anndata)
  
  bulk_vars = gene_symbols
  
  ad <- read_h5ad(RPE1_directory,backed='r') # load single cell data
  var_names = ad$var_names # get gene names
  ad = ad$obs # overwrite ad with meta data, so that the code below goes faster
  
  perturbs = unique(ad[,17]) # get perturbations
  perturbs = intersect(intersect(perturbs,var_names),bulk_vars) # perturbed, measured in single cell, measured in bulk
  ip = match(perturbs,var_names) # only perturb things that are measured
  
  perturbs = c("control",perturbs) # add in control perturbation
  
  #### BATCH DATA
  
  sc_batches = sort(unique(ad$obs[,1])) # get all single cell batches and sort them
  for (b in sc_batches){ 
    
    print(b)
    
    is = which(ad$obs[,1] == b) # samps with same batch
    is = is[which(ad$obs[is,17] %in% perturbs)]  # samps with overlapping perturbs
    
    samps = list()
    samps$data = ad$chunk_X(as.integer(is-1),replace=FALSE)
    samps$data = as.matrix(samps$data[,ip]) # do NOT include NaN for Y
    
    samps$interv = ad$obs[is,17]
    
    samps$is = is # sample labels / indices
    
    save(file=paste("samps_batch",as.character(b),"_AMD.RData",sep=""),samps)
  }
  
  
  #### ALL DATA
  
  sampsA = c(); sampsA$data = c(); sampsA$interv = c(); sampsA$is = c(); sampsA$batches = c()
  for (b in sc_batches){ 
    
    print(b)
    load(paste("samps_batch",as.character(b),"_AMD.RData",sep=""))
    
    sampsA$data = rbind(sampsA$data, samps$data)
    sampsA$interv = c(sampsA$interv, as.vector(samps$interv))
    sampsA$is = c(sampsA$is, samps$is)
    sampsA$batches = c(sampsA$batches, rep(b,length(samps$is)))
  }  
  samps = sampsA
  save(file="samps_all_AMD.RData",samps)
  
  
  ### INTERVENTIONAL DATA
  
  load("samps_all.RData")
  sampsA = samps
  for (j in 0:(length(perturbs)-1)){ 
    
    ix = which(sampsA$interv == perturbs[j+1])
    
    samps = list(); samps$data = sampsA$data[ix,]; 
    samps$batches = sampsA$batches[ix]
    samps$interv = sampsA$interv[ix]
    samps$is = sampsA$is[ix]
    
    save(file=paste("samps_interv",as.character(j),"_AMD.RData",sep=""),samps) #data, batches
    
  }
  
  
}