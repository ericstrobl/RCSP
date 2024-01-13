library(earth)
library(pcalg)
library(fastcluster)

reps = 100

Gs = vector("list",reps)
RCSP_res = Gs
reg_res = Gs
correl_res = Gs
ANM_res = Gs
lin_res = Gs

nvar = 2500
nHK = 0
nsamps = 200 #nsamps per interv
sc_nbatch = 100
for (i in 1:reps){
  print(i)
  DAG = generate_DAG_big_same4(p=nvar,en=2,nHK=nHK,sc_nbatch=sc_nbatch) #nHK = 100
  save_samps_by_file_mult3(DAG, nsamps=nsamps) # single-cell data
  
  samps = sample_DAG_NB_linear(nsamps*100,DAG$DAGb) # bulk data only
  samps$data = samps$data[1:nsamps,]
  samps$batches = samps$batches[1:nsamps]
  save(file="samps_bulk_synth.RData",samps) # bulk data
  
  Gs[[i]] = DAG
  
  print("ground truth")
  tr = get_ground_truth_bulk2_sparse(DAG$DAGs,1:(ncol(DAG$DAGs$graph)-1))
  tr$REs = tr$REs[1:nsamps,,drop=FALSE]
  tr$tREs = tr$tREs[1:nsamps,,drop=FALSE]
  
  print("find descendants")
  ptm <- proc.time()
  desL = find_des_final(samps)
  Gs[[i]]$time_desL = (proc.time() - ptm)[3]
  
  ## find stable variables
  stable = c()
  
  load("samps_bulk_synth.RData")
  
  
  print("RCSP")
  ptm <- proc.time()
  shaps = RCSP(samps,desL,reg="MARS")  ##
  RCSP_res[[i]]DAG$time = (proc.time() - ptm)[3]
  RCSP_res[[i]]DAG$RMSE_RE = compute_RMSE3(tr,shaps,ncol(samps$data))$RMSE_REs
  RCSP_res[[i]]DAG$clusters = root_SS3(tr,shaps,ncol(samps$data))
  
  print("univariate regression")
  ptm <- proc.time()
  shaps = cor_norm(samps,desL,stable) ##
  correl_res[[i]]DAG$time = (proc.time() - ptm)[3]
  correl_res[[i]]DAG$RMSE_RE = compute_RMSE3(tr,shaps,ncol(samps$data))$RMSE_REs
  correl_res[[i]]DAG$clusters = root_SS3(tr,shaps,ncol(samps$data))
  
  print("multivariate regression")
  ptm <- proc.time()
  shaps = reg_norm(samps,desL,stable) ##
  reg_res[[i]]DAG$time = (proc.time() - ptm)[3]
  reg_res[[i]]DAG$RMSE_RE = compute_RMSE3(tr,shaps,ncol(samps$data))$RMSE_REs
  reg_res[[i]]DAG$clusters = root_SS3(tr,shaps,ncol(samps$data))
  
  print("ANM")
  ptm <- proc.time()
  shaps = ANM_norm(samps,desL,stable) ##
  ANM_res[[i]]DAG$time = (proc.time() - ptm)[3]
  ANM_res[[i]]DAG$RMSE_RE = compute_RMSE3(tr,shaps,ncol(samps$data))$RMSE_REs
  ANM_res[[i]]DAG$clusters = root_SS3(tr,shaps,ncol(samps$data))
  
  print("LiNGAM")
  ptm <- proc.time()
  shaps = lin_norm(samps,desL,stable) ##
  lin_res[[i]]DAG$time = (proc.time() - ptm)[3]
  lin_res[[i]]DAG$RMSE_RE = compute_RMSE3(tr,shaps,ncol(samps$data))$RMSE_REs
  lin_res[[i]]DAG$clusters = root_SS3(tr,shaps,ncol(samps$data))
  
  save(file="Results_synth.RData", Gs, RCSP_res, reg_res, correl_res, ANM_res, lin_res)
  
}

#RMSE
imax = 30
scores = matrix(0,imax,5)
for (i in 1:imax){
  scores[i,1] = RCSP_res[[i]]DAG$RMSE_RE
  scores[i,2] = reg_res[[i]]DAG$RMSE_RE
  scores[i,3] = correl_res[[i]]DAG$RMSE_RE
  scores[i,4] = ANM_res[[i]]DAG$RMSE_RE
  scores[i,5] = lin_res[[i]]DAG$RMSE_RE
  
}
colMeans(scores)
colMeans(scores)+apply(scores,2,sd)*1.96/sqrt(imax)
colMeans(scores)-apply(scores,2,sd)*1.96/sqrt(imax)


# time
imax = 30
scores = matrix(0,imax,5)
for (i in 1:imax){
  scores[i,1] = RCSP_res[[i]]DAG$time + Gs[[i]]$time_desL
  scores[i,2] = reg_res[[i]]DAG$time
  scores[i,3] = correl_res[[i]]DAG$time
  scores[i,4] = ANM_res[[i]]DAG$time + Gs[[i]]$time_desL
  scores[i,5] = lin_res[[i]]DAG$time + Gs[[i]]$time_desL
  
}
colMeans(scores)
colMeans(scores)+apply(scores,2,sd)*1.96/sqrt(imax)
colMeans(scores)-apply(scores,2,sd)*1.96/sqrt(imax)

write.csv(file="time_synth.csv",scores)

# RMSE clusters
scores = array(0,dim=c(imax,5,10))
for (i in 1:imax){
  scores[i,1,] = RCSP_res[[i]]DAG$clusters$RMSE_RE
  scores[i,2,] = reg_res[[i]]DAG$clusters$RMSE_RE
  scores[i,3,] = correl_res[[i]]DAG$clusters$RMSE_RE
  scores[i,4,] = ANM_res[[i]]DAG$clusters$RMSE_RE
  scores[i,5,] = lin_res[[i]]DAG$clusters$RMSE_RE
  
}
apply(scores,c(2,3),mean)
write.csv(file="RCSP_clusters_synth.csv",scores[,1,])
