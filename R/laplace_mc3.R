#' @useDynLib laplacemc3
#' @useDynLib structmcmc
#' @export
#' @description Lorem ipsum dolor sit amet, consectetur adipiscing elit,
#' sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.
#' Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.
#' @title laplace_mc3
#' @title Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna...
#' @param ss ...
#' @param dag ...
#' @return ...
#' @examples
#'
#' #laplace_mc3(ss)
#'

laplace_mc3<-function(ss,itermax=5000,burnin=1000,maxParents=ss$p-1,constraintT=t(matrix(0,ss$p,ss$p)),verbose=TRUE){
  res=list(dag=list(itermax,mode="list"),
           score=rep(NA,itermax),
           accepted=rep(NA,itermax))
  rfit<-list(itermax,mode="list")
  p=ss$p
  if (verbose==TRUE){
    start=proc.time()
    cat("Running",itermax,"MC3 iterations ...")
    pb=txtProgressBar(max=itermax,style=3)
    setTxtProgressBar(pb,0)
  }
  initial=sampleBN(n=p,maxNumberParents=3)
  dag0=initial
  fit0=estim_gbn(ss,dag0)
  routes0=routes(dag0)
  adj0=as.adjacency(dag0)
  logMoves0=logNumMHNeighbours(routes0,adj0,constraintT,maxParents)
  for (iter in 1:itermax) {
    # proposal move
    prop=proposal(routes0,adj0,constraintT,maxParents)
    # update proposal dag and aux data
    dag1=prop$dag
    routes1=prop$routes
    adj1=prop$adj
    logMoves1=logNumMHNeighbours(routes1,adj1,constraintT,maxParents)
    J=prop$changed
    # update fitting using the fast update function
    fit1=update_gbn(ss,dag1,fit0,J)
    # acceptance rate
    logAR=sum(fit1$loglik[J]+fit1$laplace[J]-fit0$loglik[J]-fit0$laplace[J])-logMoves1+logMoves0

    if (!is.nan(logAR)){
      accepted=runif(1)<min(1,exp(logAR))
    } else {
      accepted=0
    }

    if (accepted) {
      dag0=dag1
      fit0=fit1
      routes0=routes1
      adj0=adj1
      logMoves0=logMoves1
    }

    res$accepted[iter]=accepted
    res$dag[[iter]]=dag0
    res$score[iter]=sum(fit0$loglik+fit0$laplace)
    rfit[[iter]]=list(s=fit0$s,w=fit0$w,detH=fit0$detH)
  if (verbose==TRUE) setTxtProgressBar(pb,value=iter);
  }
  if (verbose==TRUE){
    elapsed=(proc.time()-start)[[3]]
    cat("Running time",elapsed, "seconds\n")
  }
  ratio=mean(res$accepted)
  return(list(res=res,ratio=ratio,fit=rfit))
}
