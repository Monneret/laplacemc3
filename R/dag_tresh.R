#' @useDynLib laplacemc3
#' @export
#' @description Lorem ipsum dolor sit amet, consectetur adipiscing elit,
#' sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.
#' Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.
#' @title dag_tresh
#' @title Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna...
#' @param p ...
#' @param res ...
#' @param tresh ...
#' @return ...
#' @examples
#'
#' #dag_tresh(ss,res,tresh)
#'

dag_tresh<-function(p,res,tresh=0.85){
  ep=matrix(apply(sapply(res$dag[1000:length(res$score)],function(z) c(as.adjacency(z))),1,mean),p,p)
  consensus.dag=as.bn(ep>threshold)
  return(ep=ep,consensus.dag=consensus.dag)
}
