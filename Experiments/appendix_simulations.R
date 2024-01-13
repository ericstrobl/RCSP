require(earth)


## housekeeping genes

nvar = 250
nsamps = 200

nreps = 100

Ns = c(10,20,40,80,160,320,640,1280)
nHK = c(0,10,20)

RMSE = array(0,c(nreps,length(Ns),length(nHK)+1))

for (n in 1:length(Ns)){
  for (s in 1:(length(nHK)+1)){
    for (i in 1:nreps){
      
      print(i)
      
      if (s <= length(nHK)){
        DAG2 = generate_DAG_big_same4_Ns(p=nvar,en=2,nHK=nHK[s],sc_nbatch=1,bulk_nbatch=1,ns=Ns[n]) #nHK = 100
      } else{
        DAG2 = generate_DAG_big_same4_Ns(p=nvar,en=2,nHK=0,sc_nbatch=1,bulk_nbatch=1,ns=Ns[n])
      }
      
      Yi = DAG2$DAGb$Y
      # plot(as(as.matrix(DAG2$DAGb$graph),"graphNEL"))
      
      samps = sample_DAG_NB_linear(nsamps*100,DAG2$DAGb) # bulk data only
      
      AncY = isAncAll(DAG2$DAGs$graph,1:(nvar-1),DAG2$DAGs$Y) # ancestors of Y and intervened on
      if (length(AncY)==1){
        id = AncY
      } else{
        id = sample(AncY,1)
      }
      
      pa = which(as.matrix(DAG2$DAGb$graph)[,id]!=0)
      
      
      if (s==1){
        L = c()
      } else if (s<=length(nHK)){
        stable = DAG2$DAGb$HK
        L = rowSums(samps$data[1:nsamps,stable,drop=FALSE])
      } else{
        L = rowSums(samps$data[1:nsamps,,drop=FALSE])
      }
      
      ### ground truth
      psi1_1 = earth(cbind(samps$dataT[,c(pa,id)]),samps$data[,Yi])$fitted.values
      if (length(pa)>0){
        psi1_2 = earth(cbind(samps$dataT[,pa]),samps$data[,Yi])$fitted.values
      } else{
        psi1_2 = rep(mean(samps$data[,Yi]),nrow(samps$data))
      }
      psi_tr = psi1_1 - psi1_2
      psi_tr = psi_tr[1:nsamps]
      
      ### estimate
      psi1_1 = earth(cbind(samps$data[1:nsamps,c(pa,id)],L),samps$data[1:nsamps,Yi])$fitted.values
      if (length(pa)>0){
        psi1_2 = earth(cbind(samps$data[1:nsamps,pa],L),samps$data[1:nsamps,Yi])$fitted.values
      } else{
        psi1_2 = rep(mean(samps$data[1:nsamps,Yi]),nsamps)
      }
      psi_es = psi1_1 - psi1_2
      
      ###
      RMSE[i,n,s] = sqrt(mean( (psi_es - psi_tr)^2 ))
      
    }
  }
}

print(colMeans(RMSE))
apply(RMSE, c(2,3), mean)-1.96*apply(RMSE, c(2,3), sd)/sqrt(100)
apply(RMSE, c(2,3), mean)+1.96*apply(RMSE, c(2,3), sd)/sqrt(100)




### RCE vs signed RCS

nvar = 250
nsamps = c(100,200,400,800,Inf)

nreps = 100

RMSE = matrix(0,nreps,length(nsamps))
Signs = matrix(0,nreps,length(nsamps))

for (s in length(nsamps)){
  for (i in 1:nreps){
    
    print(i)
    
    DAG2 = generate_DAG_big_same4(p=nvar,en=2,nHK=0,sc_nbatch=1,bulk_nbatch=1)
    Yi = DAG2$DAGb$Y
    
    samps = sample_DAG_NB_linear(20000,DAG2$DAGb) # bulk data only
    
    AncY = isAncAll(DAG2$DAGs$graph,1:(nvar-1),DAG2$DAGs$Y) # ancestors of Y and intervened on
    if (length(AncY)==1){
      id = AncY
    } else{
      id = sample(AncY,1)
    }
    
    pa = which(as.matrix(DAG2$DAGb$graph)[,id]!=0)
    
    ### ground truth
    # RCE
    psi1_1 = earth(cbind(samps$dataT[,id]),samps$data[,Yi])$fitted.values
    psi1_2 = rep(mean(samps$data[,Yi]),nrow(samps$data))
    psi_tr = psi1_1 - psi1_2
    if (nsamps[s]!=Inf){
      psi_tr = psi_tr[1:nsamps[s]]
    }
    
    # RCS
    if (nsamps[s]!=Inf){
      psi1_1 = earth(cbind(samps$data[1:nsamps[s],c(pa,id)]),samps$data[1:nsamps[s],Yi])$fitted.values
      if (length(pa)>0){
        psi1_2 = earth(cbind(samps$data[1:nsamps[s],pa]),samps$data[1:nsamps[s],Yi])$fitted.values
      } else{
        psi1_2 = rep(mean(samps$data[1:nsamps[s],Yi]),nsamps[s])
      }
      psi_es = psi1_1 - psi1_2
      psi_es = psi_es[1:nsamps[s]]
    } else{
      psi1_1 = earth(cbind(samps$dataT[,c(pa,id)]),samps$data[,Yi])$fitted.values
      if (length(pa)>0){
        psi1_2 = earth(cbind(samps$dataT[,pa]),samps$data[,Yi])$fitted.values
      } else{
        psi1_2 = rep(mean(samps$data[,Yi]),nrow(samps$data))
      }
      psi_es = psi1_1 - psi1_2
    }
    
    ###
    RMSE[i,s] = sqrt(mean( (psi_es - psi_tr)^2 ))
    Signs[i,s] = sum(sign(psi_es) == sign(psi_tr))/length(psi_tr)
    
  }
}

print(colMeans(RMSE))
print(colMeans(Signs))


### LiNGAM and ANM

nvar = 100
nHK = 5

nreps = 100
ANM_res = matrix(0,200,nreps)

nsamps = c(100,200,400,800,1600,3200,6400)

RMSE_ANM = matrix(0,nreps,length(nsamps))
RMSE_LiNGAM = matrix(0,nreps,length(nsamps))

for (s in 1:length(nsamps)){
  for (i in 1:nreps){
    
    print(i)
    DAG2 = generate_DAG_big_same4(p=nvar,en=2,nHK=nHK,sc_nbatch=1,bulk_nbatch=1) #nHK = 100
    
    # plot(as(as.matrix(DAG2$DAGb$graph),"graphNEL"))
    
    samps = sample_DAG_NB_linear(nsamps[s],DAG2$DAGb) # bulk data only
    
    rel = which(colSums(as.matrix(DAG2$DAGb$graph))>=1)
    rel = setdiff(rel,DAG2$DAGb$Y)
    ind = sample(rel,1)
    pa = which(as.matrix(DAG2$DAGb$graph[,ind])!=0)
    
    stable = DAG2$DAGb$HK
    
    ### ANM
    resid = samps$data[,ind]-earth(samps$data[,pa],samps$data[,ind])$fitted.values
    RMSE_ANM[i,s] = sqrt(mean( (normalizeData(resid) - normalizeData(samps$dataE[,ind]))^2   ))
    
    ### LiNGAM
    resid = samps$data[,ind]-lm.fit( cbind(samps$data[,pa],1),samps$data[,ind])$fitted.values
    RMSE_LiNGAM[i,s] = sqrt(mean( (normalizeData(resid) - normalizeData(samps$dataE[,ind]))^2   ))
    
  }
}

print(colMeans(RMSE_ANM))
print(colMeans(RMSE_LiNGAM))
