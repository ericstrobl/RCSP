isDec <- function(graph,a,b,visited=rep(FALSE,nrow(graph)))
{
  # a is a descendant of b
  
  if (a %in% b){
    return(TRUE)
  }
  
  visited[a] = TRUE;
  
  adj = which(graph[,a] & !visited);
  
  out=FALSE;
  for (j in adj){
    out=isDec(graph, j, b, visited);
    if(out==TRUE){
      break;
    }
  }
  
  return(out)
  
}