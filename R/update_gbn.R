#' @useDynLib laplacemc3
#' @export
#' @description Lorem ipsum dolor sit amet, consectetur adipiscing elit,
#' sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.
#' Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.
#' @title update_gbn
#' @title Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna...
#' @param ss ...
#' @param dag ...
#' @param fit ...
#' @param J ...
#' @return ...
#' @examples
#'
#' #update_gbn(ss,dag,fit,J)
#'

update_gbn=function(ss,dag,fit,J) {
  p=ss$p
  d=sapply(dag,length) # number of edges towards j
  # sum(sapply(dag[J],length)^3) update_gbn cost

  # informative sample size by variable
  N=ss$N

  # solve w and sigma and compute detH
  w=fit$w
  sigma=fit$sigma
  detH=fit$detH
  for (j in J) {
    pa=dag[[j]]
    if (d[j]>0) {
      b=rep(0,d[j])
      A=matrix(0,nrow=d[j],ncol=d[j])
      for (i in 1:d[j]) {
        b[i]=ss$z[pa[i],j,j];
        for (ii in 1:d[j]) A[i,ii]=ss$z[pa[i],pa[ii],j];
      }
      w[[j]]=solve(A,b)
    } else {
      w[[j]]=numeric(0)
    }
    # now sigma
    aux=ss$z[j,j,j]-2*sum(w[[j]]*ss$z[pa,j,j])
    if (d[j]>0) aux=aux+sum(tcrossprod(w[[j]])*ss$z[pa,pa,j])
    sigma[j]=sqrt(aux/N[j])
    # detH
    if (d[j]>0) {
      detH[j]=det(A/sigma[j]^2)
    } else {
      detH[j]=1.0
    }
  }

  # finally the log-likelihood
  loglik=fit$loglik
  loglik[J]=-log(2*3.141593)/2*N[J]-N[J]*log(sigma[J])-N[J]/2
  laplace=fit$laplace
  laplace[J]=1/2*((1+d[J])*log(2*3.141593)-log(2*N[J])-log(detH[J]))


  return(list(dag=dag,sigma=sigma,s=log(sigma),w=w,d=d,detH=detH,
              loglik=loglik,laplace=laplace))
}
