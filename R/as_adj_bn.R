#' @useDynLib laplacemc3
#' @export
#' @description Lorem ipsum dolor sit amet, consectetur adipiscing elit,
#' sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.
#' Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.
#' @title estim_gbn
#' @title Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna...
#' @param dag ...
#' @return ...
#' @examples
#'
#' #estim_gbn(ss,dag)
#'

as.adj.bn<-function(dag){
  p<-length(dag)
  adj<-matrix(0,nrow=p,ncol=p)
  for (i in 1:p){
    adj[dag[[i]],i]=1
  }
  return(adj=adj)
}
