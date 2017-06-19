#' @useDynLib laplacemc3
#' @export
#' @description Lorem ipsum dolor sit amet, consectetur adipiscing elit,
#' sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.
#' Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.
#' @title summary_stat
#' @title Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna...
#' @param x ...
#' @param int ...
#' @return ...
#' @examples
#'
#' #summary_stat(x,int)
#'
summary_stat=function(x,int) {
  p=ncol(x)
  n=nrow(x)
  N=apply(!int,2,sum);
  K=lapply(as.list(data.frame(!int)),which)

  # centered data
  y=array(NA,dim=c(dim(x),p));
  for (j in 1:p) y[,,j]=t(t(x)-apply(x[K[[j]],],2,mean));


  z=array(NA,dim=c(p,p,p))
  for (j in 1:p) z[,,j]=crossprod(y[K[[j]],,j])

  return(list(n=n,p=p,N=N,z=z))
}
