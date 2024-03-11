save_data_chunks_K562 <- function(gene_symbols,K562_directory){
  # gene symbols = gene symbols of bulk RNA-seq data
  # K562 directory = directory of ReplogleWeissman2022_K562_gwps.h5ad, e.g., K562_directory = "C:/Users/ericv/Documents/ReplogleWeissman2022_K562_gwps.h5ad"
  
  require(anndata)
  
  bulk_vars = gene_symbols
  
  ad <- read_h5ad(K562_directory,backed='r') # load single cell data
  var_names = ad$var_names # get gene names
  ad = ad$obs # overwrite ad with meta data, so that the code below goes faster
  
  perturbs = unique(ad$perturbation) # get perturbations
  perturbs = intersect(intersect(perturbs,var_names),bulk_vars) # perturbed, measured in single cell, measured in bulk
  perturbsf = c()
  for (p in perturbs){
    if(sum(ad[,17]==p)>=50){ # each perturbation must have at least 50 samples
      perturbsf = c(perturbsf,p)
    }
  }
  perturbs = perturbsf
  ip = match(perturbs,var_names) # only perturb things that are measured
  perturbs = c("control",perturbs) # add in control perturbation
  
  #### BATCH DATA
  
  ad <- read_h5ad(K562_directory,backed='r')
  sc_batches = sort(unique(ad$obs$batch)) # get all single cell batches and sort them
  for (b in sc_batches){ 
    
    print(b)
    
    is = which(ad$obs$batch == b) # samps with same batch
    is = is[which(ad$obs$perturbation[is] %in% perturbs)]  # samps with overlapping perturbs
    
    samps = list()
    samps$data = ad$chunk_X(as.integer(is-1),replace=FALSE)
    samps$data = cbind(as.matrix(samps$data[,ip]),NaN); # add NaN for labels to be filled in
    
    samps$interv = ad$obs$perturbation[is]
    
    samps$is = is # sample labels / indices
    
    save(file=paste("samps_batch",as.character(b),"_MS.RData",sep=""),samps)
  }
  
  ### INTERVENTIONAL DATA
  
  interv_new = rep(TRUE, length(perturbs)) ## USE THIS ONE
  for (b in sc_batches){ 
    print(b)
    load(paste("samps_batch",as.character(b),"_MS.RData",sep=""))
    sampsb = samps
    for (j in 0:(length(perturbs)-1)){ 
      
      if(interv_new[j+1]){
        samps = list(); samps$batches = c(); samps$interv = c(); samps$is = c(); samps$data = c()
      } else{
        load(paste("samps_interv",as.character(j),"_MS.RData",sep=""))
      }
      interv_new[j+1] = FALSE
      
      ix = which(sampsb$interv == perturbs[j+1])
      
      samps$data = rbind(samps$data, sampsb$data[ix,])
      samps$batches = c(samps$batches, rep(b,length(ix)))
      samps$interv = c(samps$interv, sampsb$interv[ix])
      samps$is = c(samps$is, sampsb$is[ix])
      
      save(file=paste("samps_interv",as.character(j),"_MS.RData",sep=""),samps) #data, batches
    }
  }  
  
  
}