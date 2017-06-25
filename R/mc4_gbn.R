#' @useDynLib laplacemc3
#' @useDynLib structmcmc
#' @export
#' @description Lorem ipsum dolor sit amet, consectetur adipiscing elit,
#' sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.
#' Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.
#' @title mc4_gbn
#' @title Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna...
#' @param ss ...
#' @param beta ...
#' @param nswap ...
#' @param swap.interval ...
#' @param initial ...
#' @param maxParents ...
#' @param constraintT ...
#' @param parallel ...
#' @param cores ...
#' @return ...
#' @examples
#'
#' #laplace_mc3(ss)
#'

mc4_gbn=function(ss,beta=1.0,nswap=100,swap.interval=50,initial=NULL,maxParents=ssp$p-1,constraintT=matrix(0,ss$p,ss$p),parallel=FALSE,cores=NULL) {

  K=length(beta)

  start=proc.time()
  cat("Running MC4 with",K,"chains of",nswap*swap.interval,"iterations ...")
  pb=txtProgressBar(max=nswap,style=3)
  setTxtProgressBar(pb,0)

  res=vector(nswap,mode="list")
  nproposed.swap=0
  naccepted.swap=0

  if (is.null(initial)) initial=lapply(vector(K,mode="list"),function(z) sampleBN(n=ss$p,maxNumberParents=3))
  current=vector(K,mode="list")

  for (k in 1:K) current[[k]]=mc3_init(ss,initial[[k]],beta[k])

  if (parallel) {
    require(doParallel)
    registerDoParallel(cores=cores)
  }

  for (iter in 1:nswap) {

    if (!parallel) {
      current=lapply(current,mc3_update,ss,swap.interval)
    } else {
      current=foreach(k=1:K,
                      .export=c('as.bn'),
                      .packages='parental') %dopar% {
                        #list(res=sqrt(k))
                        mc3_update(current[[k]],ss,swap.interval)
                      }
    }
    # propose swap
    pos=sample.int(K-1,size=1)
    score0=sum(current[[pos]]$fit0$loglik)+sum(current[[pos]]$fit0$laplace)
    score1=sum(current[[pos+1]]$fit0$loglik)+sum(current[[pos+1]]$fit0$laplace)

    logAR=-(score1-score0)*(beta[pos+1]-beta[pos])

    accepted=runif(1)<min(1,exp(logAR))

    nproposed.swap=nproposed.swap+1
    naccepted.swap=naccepted.swap+accepted

    if (accepted) current[pos+0:1]=current[pos+1:0]

    res[[iter]]=current

    setTxtProgressBar(pb,value=iter)
  }

  if (parallel) stopImplicitCluster()

  elapsed=(proc.time()-start)[[3]]
  cat("Running time",elapsed, "seconds\n")


  return(list(res=res,beta=beta,elapsed=elapsed,nproposed.swap=nproposed.swap,naccepted.swap=naccepted.swap))
}
