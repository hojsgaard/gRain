## #################################################################
##
## setEvidence etc.
## 
## The old setFinding type functions (in the Finding.R file) are just
## remapped.
##
## setEvidence with the syntax below is used in Scutari's bnlearn
## book.
##
## As of August 2016, these functions call a new implementation.
##
#################################################################

## FIXME: setEvidence: Could check that invalid nodes and/or invalid states are not given.
## FIXME: setEvidence: Functions .evidenceHasValidStates and .evidenceHasValidNodes
## FIXME: setEvidence: Does the checking. For now it has just been noted in the .Rd files
## FIXME: setEvidence: invalid states / nodes are ignored

#' @title Set evidence.
#' 
#' @description Set, update and remove evidence..
#' 
#' @name grain-evidence
#' 
#' @aliases setEvidence retractEvidence absorbEvidence
#' @param object A "grain" object
#' @param nodes A vector of nodes; those nodes for which the
#'     (conditional) distribution is requested.
#' @param states A vector of states (of the nodes given by 'nodes')
#' @param evidence An alternative way of specifying findings
#'     (evidence), see examples below.
#' @param nslist deprecated
#' @param propagate Should the network be propagated?
#' @param details Debugging information
#' @return A list of tables with potentials.
#' @note \code{setEvidence()} is an improvement of \code{setFinding()}
#'     (and as such \code{setFinding} is obsolete). Users are
#'     recommended to use \code{setEvidence()} in the future.
#' 
#' \code{setEvidence()} allows to specification of "hard evidence" (specific
#' values for variables) and likelihood evidence (also known as virtual
#' evidence) for variables.
#' 
#' The syntax of \code{setEvidence()} may change in the future.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{setFinding}} \code{\link{getFinding}}
#'     \code{\link{retractFinding}} \code{\link{pFinding}}
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords models utilities
#' @examples
#' 
#' testfile <- system.file("huginex", "chest_clinic.net", package = "gRain")
#' chest <- loadHuginNet(testfile, details=0)
#' qb <- querygrain(chest)
#' qb
#' 
#' lapply(qb, as.numeric) # Safe
#' sapply(qb, as.numeric) # Risky
#' 
#' ## setFinding / setEvidence
#' 
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
#' 
#' ## 1) These two forms are identical
#' setEvidence(chest, c("asia","xray"), c("yes", "yes"))
#' setFinding(chest, c("asia","xray"), c("yes", "yes"))
#' 
#' ## 2) Suppose we do not know with certainty whether a patient has
#' ## recently been to Asia. We can then introduce a new variable
#' ## "guess.asia" with "asia" as its only parent. Suppose
#' ## p(guess.asia=yes|asia=yes)=.8 and p(guess.asia=yes|asia=no)=.1
#' ## If the patient is e.g. unusually tanned we may set
#' ## guess.asia=yes and propagate. This corresponds to modifying the
#' ## model by the likelihood (0.8, 0.1) as
#' setEvidence(chest, c("asia","xray"), list(c(0.8,0.1), "yes"))
#' 
#' ## 3) Hence, the same result as in 1) can be obtained with
#' setEvidence(chest, c("asia","xray"), list(c(1, 0), "yes"))
#' 
#' ## 4) An alternative specification using evidence is
#' setEvidence(chest, evidence=list("asia"=c(1, 0), "xray"="yes"))
#' 

#' @name grain-evidence
setEvidence <- function(object, nodes=NULL, states=NULL, evidence=NULL, nslist=NULL,
                        propagate=TRUE, details=0){
    
    if (!is.null(nslist))
        stop("Argument 'nslist' has been deprecated; please use 'evidence' instead\n")

    if ( is.null( evidence ) && is.null( nodes ) )
        stop( "Evidence is not given; nothing to do...")
    
    if ( is.null( evidence ) ) ## Then 'nodes' must be given
        evidence <- .nodes.states2evidence( nodes, states )      
    
    ##setEvidence_(object, evidence, propagate=propagate, details=details)
    setEvi_(object, evidence, propagate=propagate, details=details)
    
}

#' @name grain-evidence
retractEvidence <- function(object, nodes=NULL, propagate=TRUE){
    ##.retractEvidence_internal(object, nodes=nodes, propagate=propagate)
    retractEvi_(object, items=nodes, propagate=propagate)
}

#' @name grain-evidence
absorbEvidence <- function(object, propagate=TRUE ){
    absorbEvi( object, propagate=propagate )
}




### #############################################################
##
##  PROBABLY OBSOLETE BELOW HERE
##
### #############################################################



## ## @rdname grain-evidence
## ## @param evi.list A list of evidence.
## ## @param levels A list of levels of variables.
## newEvidence <- function(evi.list, levels){

##     ## First remove all evidence specified as NA
##     not.na <- !unlist(lapply(lapply(evi.list,is.na), any), use.names=FALSE)
##     if (length( not.na ) > 0)
##         evi.list <- evi.list[ not.na ]

##     evidence           <- vector("list", length(evi.list))
##     is.hard.evidence   <- rep.int(TRUE,  length(evi.list))
##     hard.state         <- rep.int(NA,    length(evi.list))

##     ## cat("newEvidence\n"); print(evi.list)

##     for (i in seq_along(evi.list)){
##         ev <- evi.list[i]
##         v <- ev[[1]]

##         if( is.array(v)){
##             n <- names(dimnames(v))
##             is.hard.evidence[i]    <- FALSE
##             evidence[[i]] <- v
##             next
##         }
        
##         if (is.character(v)){
##             n <- names(evi.list)[i]
##             hard.state[i]  <- v
##             evidence[[i]]  <- .hard.state2parray(n, v, levels[[n]])
##             next
##         }
        
##         if (is.numeric(v)){
##             n <- names(evi.list)[i]
##             is.hard.evidence[i] <- FALSE
##             evidence[[i]] <- .soft.state2parray(n, v, levels[[n]])
##         }
##     }

##     ## print(evidence)
##     ## If evidence is zero on all states or negative on some (or all) states then it is invalid
##     keep <- unlist(lapply(evidence, function(e){ sum(e) !=0 && all(e>=0) }), use.names=FALSE)
##     ## print(keep)

##     nodes <- unique.default( unlist(lapply(evidence, .namesDimnames)),
##                             use.names=FALSE )

##     out <- list(summary=list(
##                     nodes=nodes[keep],
##                     is.hard.evidence=is.hard.evidence[keep],
##                     hard.state=hard.state[keep]),
##                 evidence=evidence[keep])

##     class( out ) <- "grainEvidence_"
##     ## print.default( out )
##     out
## }



## ## setEvidence_: The workhorse
## setEvidence_ <- function(object, evidence=NULL, propagate=TRUE, details=0){
##                                         #
##     #details=1
##     # cat("++++ setEvidence_\n"); print(evidence)
    
##     old.ev <- getEvidence(object)
    
##     if (.is.JointEvidence( old.ev ))
##         stop("JointEvidence (multivariate evidence) has been given;\n   to update evidence use setJointEvidence\n")
    
##     if ( !.is.Evidence(evidence) ){
##         ## Relevant if someone sticks in a row of a dataframe; not sure if this is a good solution.
##         if (class(evidence)=="data.frame"){
##             if (nrow(evidence)>1)
##                 stop("evidence is a data.frame with more than one row; only one row is allowed\n")
##             evidence <- lapply(evidence, as.character)
##         }
        
##         ## Strip any NA's in the evidence; hmmm - maybe I do this twice ?? FIXME ??
##         idx <- !unlist(lapply(evidence, function(e) any(is.na(e))), use.names = FALSE)
##         if (length(idx)>0)
##             evidence <- evidence[ idx ]
        
##         new.ev <- newEvidence(evidence, object$universe$levels)
##     } else {
##         new.ev <- evidence
##     }

##     ## cat("Evidence - later\n"); print(new.ev)

##     tot.ev <- new.ev # can be changed below

##     if (details>0){
##         cat("new.ev (nodes):"); print(new.ev$summary$nodes)
##         cat("old.ev (nodes):"); print(old.ev$summary$nodes)
##     }

##     if (!object$isCompiled){
##         object <- compile(object)
##     }
    
##     object$isInitialized  <- FALSE
##     object$isPropagated   <- FALSE

##     if ( length(old.ev) > 0 ){
##         ## nodes on which there is already evidence will not be given
##         ## new evidence
##         idx <- match( intersect(new.ev$summary$nodes, old.ev$summary$nodes), new.ev$summary$nodes)
##         if ( length( idx ) > 0 ){
##             new.ev <- .delete.evidence( new.ev, idx)
##             if (details>0){
##                 cat("new.ev (nodes - updated):")
##                 print(new.ev$summary$nodes)
##             }
##         }
##         tot.ev <- .append.evidence( old.ev, new.ev )
##     }

##     if ( length( new.ev ) > 0 ){
##         if (details>0){
##             cat("tot.ev (nodes):"); print(tot.ev$summary$nodes)
##             cat("new.ev (to be inserted):"); print(new.ev$summmary$nodes)
##         }
##         host  <- .get.host.clique( new.ev$evidence, object$rip )
##         object$temppot <- .insert.evidence.in.potential( object$temppot, new.ev, host )
##         object$evidence <- tot.ev
##     }

##     if (details>0){
##         cat("after insertion:\n"); print( object$evidence )
##     }

##     if (propagate) propagate(object) else object
## }




## .insert.evidence.in.potential <- function( pot, evi.list, hostclique ){
##     ## if (any(is.na(hostclique))) stop("NAs in hostclique...")
##     for (i in seq_along( evi.list$evidence ) ){
##         j <- hostclique[ i ]
##         p <- evi.list$evidence[[ i ]]
##         pot[[j]] <- tabMult__( pot[[ j ]], p )
##     }
##     pot
## }





## .retractEvidence_internal <- function(object, nodes=NULL, propagate=TRUE){
##     #cat("++++ retractEvidence_\n")
##     .resetgrain <- function(x){
##         x$temppot       <- x$origpot
##         x$evidence       <- NULL
##         x$isPropagated  <- FALSE
##         x
##     }

##     if ( is.null( nodes ) ){
##         object <- .resetgrain( object )
##     } else {
##         old.ev <- getEvidence(object)
##         if ( .is.JointEvidence(old.ev) )
##             stop("JointEvidence (Multivariate evidence) has been set;\n   to retract retractJointEvidence instead...")

##         #cat("old.ev:\n"); print(old.ev)
##         idx <- match(intersect(old.ev$summary$nodes, nodes), old.ev$summary$nodes)
##         #print(idx)
##         if (length(idx)>0){
##             new.ev <- .delete.evidence( old.ev, idx )
##             #cat("new.ev:\n"); print(new.ev)
##             object <- .resetgrain(object)
##             if ( length(new.ev$summary$nodes) > 0 ){
##                 object <- setEvidence_(object, evidence=new.ev, propagate=FALSE)
##             }
##         }
##     }

##     if (propagate){
##         propagate(object)
##     } else {
##         object
##     }
## }







## ## Check for invalid states
## .evidenceHasValidStates <- function( evidence, universe ){
##     out <- TRUE
##     for (i in 1:length(evidence)){
##         n <- names(evidence)[ i ]
##         e <- evidence[[ i ]]
        
##         (j <- match( e, universe$levels[[ n ]] ) )
##         if ( any( is.na( j ) ) ){
##             out <- FALSE
##             st <- sprintf("Invalid states; node: %s states: %s\n", n,
##                           toString( e[ which( is.na(j) ) ] ))
##             cat(st)
##         }
##     }
##     out
## }

## ## Check for invalid nodes
## .evidenceHasValidNodes <- function(evidence, universe){
##     out <- TRUE
##     (j <- match( names( evidence ), universe$nodes ))
##     if ( any( is.na( j ) ) ){
##         st <- sprintf("Node(s) are not in the network: %s\n",
##                       toString( names( evidence )[ which( is.na(j) ) ] ))
##         cat(st)
##         out <- FALSE
##     }
##     out
## }




## ## ####################################################################
## ##
## ## Evidence related utility functions
## ##
## ## ####################################################################


## .is.Evidence <- function(x){
##     !is.null(x) && class(x)=="grainEvidence_"
## }

## .is.JointEvidence <- function(x){
##     !is.null(x) && class(x)=="grainJointEvidence_"
## }


## .append.evidence <- function(ev1, ev2){
##     if (is.null(ev1))
##         return (ev2)
##     if (is.null(ev2))
##         return (ev1)

##     #cat(".append.evidence:\n"); cat("ev1:\n"); print(ev1); cat("ev2:\n"); print(ev2)
##     summary <- lapply(seq_along(ev1$summary),
##                       function(i) c(ev1$summary[[i]], ev2$summary[[i]]))
##     names(summary) <- names(ev1$summary)
##     evidence <- c(ev1$evidence, ev2$evidence)
##     out <- list(summary=summary, evidence=evidence)
##     class(out) <- "grainEvidence_"
##     out
## }


## .delete.evidence <- function(ev, idx){
##     if (length(idx)>0){
##         summary <- lapply(ev$summary, function(x) x[-idx])
##         names(summary) <- names(ev$summary)
##         ev <- list(summary=summary, evidence=ev$evidence[-idx])
##         #ev <- lapply(ev, function(x) x[-idx])
##         class(ev) <- "grainEvidence_"
##     }
##     ev
## }




## ## ####################################################################
## ##
## ## Ting der skal omdoebes / et andet sted hen ved lejlighed
## ##
##  ## ####################################################################











