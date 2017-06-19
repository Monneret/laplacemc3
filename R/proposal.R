#' @useDynLib laplacemc3
#' @useDynLib structmcmc
#' @useDynLib parental
#' @export
#' @description Lorem ipsum dolor sit amet, consectetur adipiscing elit,
#' sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.
#' Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.
#' @title proposal
#' @title Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna...
#' @param routes0 ...
#' @param adj0 ...
#' @param constraintT ...
#' @param maxParents ...
#' @return ...
#' @examples
#'
#' #proposal(routes0,adj0,constraintT,maxParents)
#'

proposal=function(routes0,adj0,constraintT,maxParents) {
  canAddOrRemove=transposeEdgeIsTogglable(routes0,adj0,constraintT,maxParents)
  canFlip=edgeIsFlippable(routes0,adj0,constraintT,maxParents)

  nonCycleInducing=which(canAddOrRemove, arr.ind = T)
  nonCycleInducingFlips=which(canFlip, arr.ind = T)

  nNonCycleInducing=nrow(nonCycleInducing)
  nNonCycleInducingFlips=nrow(nonCycleInducingFlips)
  numberOfMoves=nNonCycleInducing + nNonCycleInducingFlips

  select=sample.int(numberOfMoves,size=1)
  adj1=adj0
  routes1=routes0
  if (select <= nNonCycleInducing){
    # changing edge status (careful with the transposition)
    j=nonCycleInducing[[select, 1]]
    i=nonCycleInducing[[select, 2]]
    # only pa[[j]] are changed
    changed=j
    if (adj1[i,j]==1) {
      adj1[i,j]=0
      routes1=routesRemoveEdge(routes1,i,j)
    } else {
      adj1[i,j]=1
      routes1=routesAddEdge(routes1,i,j)
    }
  } else {
    # flipping edge
    i=nonCycleInducingFlips[[select-nNonCycleInducing, 1]]
    j=nonCycleInducingFlips[[select-nNonCycleInducing, 2]]
    # both pa[[i]] and pa[[j]] are changed
    changed=c(i,j)
    adj1[i,j]=0
    adj1[j,i]=1
    routes1=routesRemoveEdge(routes1,i,j)
    routes1=routesAddEdge(routes1,j,i)
  }
  return(list(dag=as.bn(adj1),adj=adj1,routes=routes1,changed=changed))
}
