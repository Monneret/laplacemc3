#' @useDynLib laplacemc3
#' @export
#' @description Lorem ipsum dolor sit amet, consectetur adipiscing elit,
#' sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.
#' Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.
#' @title estim_gbn
#' @title Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna...
#' @param ss ...
#' @param dag ...
#' @return ...
#' @examples
#'
#' #estim_gbn(ss,dag)
#'
estim_gbn=function(ss,dag) {
  p=ss$p
  d=sapply(dag,length) # number of edges towards j

  # informative sample size by variable
  N=ss$N

  # solve w and sigma and compute detH
  w=vector(p,mode="list")
  sigma=rep(NA,p)
  detH=rep(1.0,p)
  for (j in 1:p) {
    pa=dag[[j]]
    if (d[j]>0) {
      b=rep(0,d[j])
      A=matrix(0,nrow=d[j],ncol=d[j])
      for (i in 1:d[j]) {
        b[i]=ss$z[pa[i],j,j];
        for (ii in 1:d[j]) A[i,ii]=ss$z[pa[i],pa[ii],j];
      }
      w[[j]]=solve(A,b)
    }
    # now sigma
    aux=ss$z[j,j,j]-2*sum(w[[j]]*ss$z[pa,j,j])
    if (d[j]>0) aux=aux+sum(tcrossprod(w[[j]])*ss$z[pa,pa,j])
    sigma[j]=sqrt(aux/N[j])
    # detH
    if (d[j]>0) detH[j]=det(A/sigma[j]^2)
  }

  # finally the log-likelihood
  loglik=-log(2*3.141593)/2*N-N*log(sigma)-N/2

  laplace=1/2*((1+d)*log(2*3.141593)-log(2*N)-log(detH))


  return(list(dag=dag,sigma=sigma,s=log(sigma),w=w,d=d,detH=detH,
              loglik=loglik,laplace=laplace))

}
