## #################################################################
##
## setFinding etc:
##
## Must live to ensure backward compatibility (it is used in the GMwR
## book, in the JSS paper and in many places on the web
##
## The functions call the corresponding evidence functions,
## setEvidence etc.
##
## #################################################################

#' @title Set, retrieve, and retract finding in Bayesian network.
#' 
#' @description Set, retrieve, and retract finding in Bayesian
#'     network.  NOTICE: The functions described here are kept only
#'     for backward compatibility; please use the corresponding
#'     evidence-functions in the future.
#' 
#' @name finding
#' 
#' @aliases setFinding retractFinding getFinding pFinding
#' @param object A "grain" object
#' @param nodes A vector of nodes
#' @param states A vector of states (of the nodes given by 'nodes')
#' @param flist An alternative way of specifying findings, see
#'     examples below.
#' @param propagate Should the network be propagated?
#' @note NOTICE: The functions described here are kept only for
#'     backward compatibility; please use the corresponding
#'     evidence-functions in the future:
#' 
#' \code{setEvidence()} is an improvement of \code{setFinding()} (and as such
#' \code{setFinding} is obsolete). Users are recommended to use
#' \code{setEvidence()} in the future.
#' 
#' \code{setEvidence()} allows to specification of "hard evidence" (specific
#' values for variables) and likelihood evidence (also known as virtual
#' evidence) for variables.
#' 
#' The syntax of \code{setEvidence()} may change in the future.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{setEvidence}} \code{\link{getEvidence}}
#' \code{\link{retractEvidence}} \code{\link{pEvidence}}
#' \code{\link{querygrain}}
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#' with the gRain Package for R. Journal of Statistical Software, 46(10), 1-26.
#' \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords models utilities
#' @examples
#' 
#' 
#' ## setFindings
#' yn <- c("yes","no")
#' a    <- cptable(~asia, values=c(1,99),levels=yn)
#' t.a  <- cptable(~tub+asia, values=c(5,95,1,99),levels=yn)
#' s    <- cptable(~smoke, values=c(5,5), levels=yn)
#' l.s  <- cptable(~lung+smoke, values=c(1,9,1,99), levels=yn)
#' b.s  <- cptable(~bronc+smoke, values=c(6,4,3,7), levels=yn)
#' e.lt <- cptable(~either+lung+tub,values=c(1,0,1,0,1,0,0,1),levels=yn)
#' x.e  <- cptable(~xray+either, values=c(98,2,5,95), levels=yn)
#' d.be <- cptable(~dysp+bronc+either, values=c(9,1,7,3,8,2,1,9), levels=yn)
#' plist <- compileCPT(list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be))
#' chest <- grain(plist)
#' 
#' ## These two forms are equivalent
#' bn1 <- setFinding(chest, nodes=c("asia","xray"), states=c("yes", "yes"))
#' bn2 <- setFinding(chest, flist=list(c("asia","yes"), c("xray", "yes")))
#' 
#' getFinding(bn1)
#' getFinding(bn2)
#' 
#' pFinding(bn1)
#' pFinding(bn2)
#' 
#' bn1 <- retractFinding(bn1, nodes="asia")
#' bn2 <- retractFinding(bn2, nodes="asia")
#' 
#' getFinding(bn1)
#' getFinding(bn2)
#' 
#' pFinding(bn1)
#' pFinding(bn2)
#' 
#' 
#' @export setFinding
setFinding <- function(object, nodes=NULL, states=NULL, flist=NULL, propagate=TRUE){
    if (!is.null(flist)){
        flist2 <- do.call("rbind",flist)
        nodes   <- flist2[,1]
        states  <- flist2[,2]
    }
    setEvidence(object, nodes=nodes, states=states, propagate=propagate)
}


retractFinding <- retractEvidence

pFinding <- pEvidence

getFinding <- getEvidence








## print.grainFinding <- function(x, ...){
##   cat("Finding: \n")
##   v<-do.call("cbind",x)
##   colnames(v) <- c("variable", "state")
##   print(v, quote=FALSE)
##   if (!is.null(attr(x,"pFinding")))
##     cat("Pr(Finding)=",attr(x,"pFinding"),"\n")
##   return(x)
## }



## setFinding <- function(object, nodes=NULL, states=NULL, flist=NULL, propagate=TRUE){

## ###cat("setFinding\n")
##     if (!object$isCompiled){
##         ##cat("setFinding: Compiling model ...\n")
##         object <- compile(object)
##     }

##     if (!is.null(flist)){
##         flist2 <- do.call("rbind",flist)
##         nodes   <- flist2[,1]
##         states  <- flist2[,2]
##     }

##     len            <- length(nodes)
##     if (len>0){
##         netNodes       <- nodeNames(object) # object$universe$nodes
##         currFinding    <- getFinding(object)
##         for (i in 1:len){
##             ev1   <- nodes[i];
##             if (!(ev1 %in% netNodes)){
##                 nodes[i] <- states[i] <- NA
##                 ##warning("Node ", ev1, " is not in network, skipping it...\n",call.=FALSE)
##             }
##         }

##         ## Ignore NA's in the states
##         ##
##         nodes  <- nodes[!is.na(states)]
##         states <- states[!is.na(states)]
## ### print(nodes); print(states)

##         ## Drop nodes which are already given evidence in the network
##         ##
##         if (!is.na(currFinding)){
##             idx <- match(intersect(currFinding$nodes, nodes), nodes)
##             if (length(idx)>0){
##                 nodes  <- nodes[-idx]
##                 states <- states[-idx]
##             }
##         }
##         ##  Now insert the findings
##         ##
##         if (length(nodes)>0){
##             t0 <- proc.time()
##             ## setFinding: findings are inserted to temppot
##             object$temppot <- .insertFinding( nodes, states, object$temppot, object$rip )
##             object$isInitialized  <- FALSE

##             if (!is.na(currFinding)){
##                 ev <- list(nodes=c(currFinding$nodes,nodes), states=c(currFinding$states,states))
##             } else {
##                 ev <- list(nodes=nodes, states=states)
##             }

##             ## Set finding slot
##             class(ev)<-"grainFinding"
##             object$finding <- ev
##             if (object$control$timing)
##                 cat("Time: enter finding", proc.time()-t0, "\n")

##             if (propagate){
##                 object<-propagate(object)
##             } else {
##                 object$isPropagated <- FALSE
##             }
##         }
##     }
##     return(object)
## }


## .insertFinding <- function(nodes, states, APlist, rip, details=0){

##     .infoPrint(details, 1, cat(".insertFinding\n"))
##     cli <- rip$cliques

##     ## Note: perhaps create amat globally
##     amat <- glist2setMAT(cli,vn=rip$nodes)

##     for (i in 1:length(nodes)){
##         currn <- nodes[i]
##         currs <- states[i]
##         ##cat("Node:", currn, "State:", toString(currs), "\n")
##         idx <- which(amat[,currn]==1)#[1]

##         ##cat("Host cliques:",paste(idx,sep=' '),"\n");
##         for (j in idx){
##             cpot <- APlist[[j]]
##             ## cat("Current clique:", paste(varNames(cpot), sep=' '),"\n")
##             ## lev    <- valueLabels.array(cpot)[[currn]] ## BRIS
##             lev    <- dimnames(cpot)[[currn]]
##             evTab  <- .findingTable(currn, currs, lev)
##             ##print(evTab)
##             APlist[[j]]  <- tableOp(cpot, evTab, "*")
##         }
##     }
##     APlist
## }

## .findingTable <- function(node, state, levels){
##     pot   <- rep.int(0,length(levels))
##     pot[match(state, levels)] <- 1
##     t2  <- parray(node, list(levels), pot)
##     t2
## }



#' ## FIXME: retractFinding virker naeppe med evidence
#' retractFinding <- function(object, nodes=NULL, propagate=TRUE){

#'     .resetgrain <- function(x){
#'         x$temppot       <- x$origpot
#'         x$finding       <- NULL
#'         x$isInitialized <- TRUE  ## FIXME Do we need this? I doubt!
#'         x$isPropagated  <- FALSE
#'         x
#'     }

#'     if (is.null(nodes)){
#'         object <- .resetgrain( object )
#'     } else {

#'         ev <- getFinding(object)
#'         evnodes   <- ev$nodes
#'         evstates  <- ev$states
#'         idx <- match(intersect(evnodes,nodes), evnodes)

#'         if (length(idx)>0){
#'             newevnodes  <- evnodes[-idx]
#'             newevstates <- evstates[-idx]
#'             object <- .resetgrain(object)
#'             if (length(newevnodes)>0){
#'                 object <- setFinding(object, nodes=newevnodes,
#'                                      states=newevstates, propagate=FALSE)
#'             }
#'         }
#'     }

#'     if (propagate){
#'         propagate(object)
#'     } else {
#'         object
#'     }
#' }



#' setEvidence <- function(object, nodes=NULL, states=NULL, nslist=NULL, propagate=TRUE){

#'     if (!object$isCompiled){
#'         object <- compile(object)
#'     }
#'     object$isInitialized <- FALSE
#'     object$isPropagated  <- FALSE

#'     if (!is.null(nslist)){
#'         nodes   <- names(nslist)
#'         states  <- nslist
#'     }

#'     if (!is.list(states))
#'         states <- as.list( states )

#'     n.nodes            <- length(nodes)
#'     if (n.nodes>0){
#'         ## Skip nodes not in network
#'         netNodes       <- nodeNames(object)
#'         for (i in 1:n.nodes){
#'             v   <- nodes[ i ];
#'             if (!(v %in% netNodes)){
#'                 nodes[ i ] <- states[ i ] <- NA
#'                 ##warning("Node ", v, " is not in network, skipping it...\n",call.=FALSE)
#'             }
#'         }

#'         ## Ignore NA's in the states
#'         nodes  <- nodes[  !is.na(states) ]
#'         states <- states[ !is.na(states) ]

#'         ## Drop nodes which are already given evidence in the network
#'         currFinding    <- getFinding(object)
#'         if (!is.null(currFinding)){
#'             idx <- match(intersect(currFinding$nodes, nodes), nodes)
#'             if ( length(idx) > 0 ){
#'                 nodes  <- nodes[-idx]
#'                 states <- states[-idx]
#'             }
#'         }

#'         ##  Insert the findings to temppot
#'         if (length(nodes)>0){
#'             t0 <- proc.time()

#'             ## Maximum value of soft evidence is 1
#'             states <- lapply(states, function(x) if (is.numeric(x)) pmin(x,1) else x)
#'             object$temppot <- .insertEvidence( nodes, states, object$temppot, object$rip )
#'             object$isInitialized  <- FALSE

#'             if (!is.null(currFinding)){
#'                 ev <- list(nodes=c(currFinding$nodes, nodes), states=c(currFinding$states, states))
#'             } else {
#'                 ev <- list(nodes=nodes, states=states)
#'             }

#'             ## Set finding slot
#'             ##class(ev)      <- "grainFinding"
#'             object$finding <- ev

#'             if (object$control$timing)
#'                 cat("Time: enter finding", proc.time()-t0, "\n")

#'             if (propagate){
#'                 object<-propagate(object)
#'             } else {
#'                 object$isPropagated <- FALSE
#'             }
#'         }
#'     }
#'     return(object)
#' }


#' .insertEvidence <- function(nodes, states, potlist, rip, details=0){

#'     .createEvidenceTable <- function(node, state, levels){
#'         ##print(node); print(state)
#'         if (is.numeric(state) && length(state)==length(levels)){
#'             .fast.parray(node, list(levels), pmin(state,1) )
#'         } else {
#'             if (is.character(state)){
#'                 ii <- match(state, levels)
#'                 if (any(is.na(ii)))
#'                     stop("Invalid state given...")
#'                 pot   <- rep.int(0,length(levels))
#'                 pot[ ii ] <- 1
#'                 .fast.parray(node, list(levels), pot)
#'             } else {
#'                 cat("Node:", node, "State:", toString(state), "\n")
#'                 stop("Can not create evidence table...\n")
#'             }
#'         }
#'     }

#'     ##cat(".insertEvidence\n"); print(nodes); print(states)
#'     .infoPrint(details, 1, cat(".insertFinding\n"))

#'     for (i in 1:length(nodes)){
#'         currn <- nodes[ i ]
#'         currs <- states[[ i ]]
#'         ##cat("Node:", currn, "State:", toString(currs), "\n")
#'         h.idx  <- rip$host[ match(currn, rip$nodes) ]
#'         cpot   <- potlist[[ h.idx ]]
#'         lev    <- dimnames(cpot)[[currn]]
#'         evtab  <- .createEvidenceTable(currn, currs, lev)
#'         potlist[[ h.idx ]]  <- tableOp(cpot, evtab, "*")
#'     }
#'     potlist
#' }


#' print.grainFinding <- function(x, ...){
#'     cat("Finding: \n")

#'     n.nodes <- length(x$nodes)
#'     len <- min(10, max(unlist(lapply(x$nodes, nchar))))
#'     for(i in 1:n.nodes){
#'         cat(sprintf("%*s: %s\n", len, x$nodes[ i ], toString(x$states[[ i ]])))
#'     }
#'     if (!is.null(attr(x,"pFinding")))
#'         cat("Pr(Finding)=",attr(x,"pFinding"),"\n")
#'     return(x)
#' }
