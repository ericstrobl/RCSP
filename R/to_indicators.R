to_indicators <- function(discrete,uniq = sort(unique(discrete)) ){
  
  ind = matrix(0,length(discrete),length(uniq)-1)
  for (u in seq_len(length(uniq)-1)){
    ind[which(discrete == uniq[u]),u]=1
  }
  
  return(ind)
}
