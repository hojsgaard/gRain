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
#' @description Set, update and remove evidence.
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
#'
#' @return A list of tables with potentials.
#'
#' @note \code{setEvidence()} is an improvement of \code{setFinding()}
#'     (and as such \code{setFinding} is obsolete). Users are
#'     recommended to use \code{setEvidence()} in the future.
#' 
#' \code{setEvidence()} allows to specification of "hard evidence" (specific
#' values for variables) and likelihood evidence (also known as virtual
#' evidence) for variables.
#' 
#' The syntax of \code{setEvidence()} may change in the future.
#'
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{setFinding}}, \code{\link{getFinding}},
#'     \code{\link{retractFinding}}, \code{\link{pFinding}}
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

#' @rdname grain-evidence
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

#' @rdname grain-evidence
retractEvidence <- function(object, nodes=NULL, propagate=TRUE){
    ##.retractEvidence_internal(object, nodes=nodes, propagate=propagate)
    retractEvi_(object, items=nodes, propagate=propagate)
}

#' @rdname grain-evidence
absorbEvidence <- function(object, propagate=TRUE ){
    absorbEvi( object, propagate=propagate )
}



