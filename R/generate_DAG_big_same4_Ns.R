generate_DAG_big_same4_Ns <- function(p,en,nHK=5,sc_nbatch=10,bulk_nbatch=3,ns=990){
  
  require(Matrix)
  
  N = p*p - p;
  
  DAGs = list()
  DAGs$Y = 0
  DAGs$HK = 0
  
  graph = Matrix(0,p,p,sparse=TRUE)
  while( (DAGs$Y != p) | (length(DAGs$HK)<nHK)  ){
    
    ### SINGLE CELL
    
    samplesB = rbinom(N/2,1, en/(p-1) ); # sample edges
    graph[upper.tri(graph, diag=FALSE)] <- samplesB; # put in edges in upper triangular
    
    ord = sample(1:p,p,replace=FALSE) # permute order
    DAGs$graph = graph[ord,ord]
    # ord2 = rep(0,p)
    # ord2[ord] = 1:p
    DAGs$ord = ord
    
    Ys = which( (rowSums(DAGs$graph)==0) & (colSums(DAGs$graph)>0) ) # no children, some observed parents
    
    HK = which( (rowSums(DAGs$graph)==0) & (colSums(DAGs$graph)==0) ) # no children and no parents
    
    if (p %in% Ys) { # ensure that DAGs$Y is the last variable
      DAGs$Y = p
      
      if ( length(HK)>=nHK ) {
        DAGs$HK = HK
      }
    }
  }
  
  # partial order
  ord = rep(0,p)
  ord[DAGs$ord] = 1:p
  DAGs$ord = ord
  
  ix = which(DAGs$graph!=0)
  weights = matrix(0.75*runif(length(ix))+0.25)*sample(c(-1,1),length(ix),replace=TRUE)
  DAGs$weights = DAGs$graph
  DAGs$weights[ix] = weights
  
  DAGs$off = runif(p)*0.4 + 0.1 #
  
  
  DAGs$n_batch = sc_nbatch
  DAGs$batch = runif(p*DAGs$n_batch)*0.9 + 0.1#0.1 to 2, volume times pi_{ij}
  DAGs$batch = matrix(DAGs$batch,DAGs$n_batch,p) #batches are rows, variables are columns
  
  ### BULK
  
  DAGb = DAGs
  
  DAGb$n_batch = bulk_nbatch
  DAGb$batch = runif(p*DAGb$n_batch)*ns + 10 #10 to 200, volume times pi_{ij}
  DAGb$batch = matrix(DAGb$batch,DAGb$n_batch,p) #batches are rows, variables are columns
  
  return(list(DAGs=DAGs,DAGb=DAGb))
}