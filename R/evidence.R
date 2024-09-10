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

## FIXME: setEvidence: Could check that invalid nodes and/or invalid
## states are not given.

## FIXME: setEvidence: Functions .evidenceHasValidStates and
## .evidenceHasValidNodes

## FIXME: setEvidence: Does the checking. For now it has just been
## noted in the .Rd files

## FIXME: setEvidence: invalid states / nodes are ignored

#' @title Set, update and remove evidence.
#' 
#' @description Set, update and remove evidence.
#' 
#' @name grain_evidence
#' 
#' @aliases setEvidence retractEvidence absorbEvidence evidence_get
#' 
#' @param object A "grain" object
#' @param evidence An alternative way of specifying findings
#'     (evidence), see examples below.
#' @param nodes A vector of nodes; those nodes for which the
#'     (conditional) distribution is requested. Now deprecated; use argument 'evidence' instead. 
#' @param states A vector of states (of the nodes given by 'nodes'). Now deprecated; use argument 'evidence' instead. 
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
#'     46(10), 1-26.  \url{https://www.jstatsoft.org/v46/i10/}.
#' @keywords models utilities
#'
#' @examples
#' 
#' example("example_chest_cpt")
#' chest.bn <- grain(compileCPT(chest_cpt))
#' chest.bn <- compile(chest.bn)
#' 
#' ## 1) These two forms are identical
#' setEvidence(chest.bn, c("asia", "xray"), c("yes", "yes"))
#' setFinding(chest.bn, c("asia", "xray"), c("yes", "yes"))
#' 
#' ## 2) Suppose we do not know with certainty whether a patient has
#' ## recently been to Asia. We can then introduce a new variable
#' ## "guess.asia" with "asia" as its only parent. Suppose
#' ## p(guess.asia=yes|asia=yes)=.8 and p(guess.asia=yes|asia=no)=.1
#' ## If the patient is e.g. unusually tanned we may set
#' ## guess.asia=yes and propagate.
#' ##
#' ## This corresponds to modifying the model by the likelihood (0.8,
#' ## 0.1) as
#' 
#' setEvidence(chest.bn, c("asia", "xray"), list(c(0.8, 0.1), "yes"))
#' 
#' ## 3) Hence, the same result as in 1) can be obtained with
#' setEvidence(chest.bn, c("asia", "xray"), list(c(1, 0), "yes"))
#' 
#' ## 4) An alternative specification using evidence is
#' setEvidence(chest.bn, evidence=list(asia=c(1, 0), xray="yes"))
#' 

#' @rdname grain_evidence
#' @export 
setEvidence <- function(object, nodes=NULL, states=NULL, evidence=NULL, 
                        propagate=TRUE, details=0){
    stopifnot_grain(object)    

    ## if (!is.null(nodes))
    ##     warning("Argument 'nodes' and 'states' deprecated, please use argument 'evidence' instead.")

    ## if (!is.null(states))
    ##     warning("Argument 'nodes' and 'states' deprecated, please use argument 'evidence' instead.")

    
    if (is.null(evidence) && is.null(nodes))
        stop("Evidence is not given; nothing to do...")
    
    if (is.null(evidence)) ## Then 'nodes' must be given
        evidence <- nodes_states_to_evidence(nodes, states)   
    
    set_evidence_worker(object, evidence, propagate=propagate, details=details)    
}


set_evidence_worker <- function(object, evidence=NULL, propagate=TRUE, details=0){
    ## details=10
    ## cat("++++ set_evidence_worker input evidence: \n"); str(evidence)

    insertEvi <- function(evi_object, pot, hostclique){
        if ( !inherits(evi_object, "grain_evidence") )
            stop("'object' is not a 'grain_evidence' object")
        
        for (i in seq_along( evi_object$evi_weight) ){
            p <- evi_object$evi_weight[[ i ]]
            j <- hostclique[ i ]
            pot[[j]] <- tabMult(pot[[ j ]], p)
        }
        pot
    }
    
    if (is.null_evi( evidence )){
        cat("Nothing to do\n")
    } else {
        if (!isCompiled(object)){
            object <- compile(object)
        } 
        
        oe <- getEvidence( object ) # Der er noget med typen (old/new)
        ne <- new_evi( evidence, universe(object)$levels )
        
        if (details > 0){
            cat("old.evidence:\n"); print.default(oe)
            cat("new evidence:\n"); print.default(ne)
        }
        
        ## Hvis der er eksisterende evidens, så skal det tages ud af det nye
        if (!is.null_evi( oe ) ){
            ne <- setdiff_evi( ne, oe )
            if (details>0) {
                cat("new evidence - after modification:\n"); print( ne )
            }
        }

        if (length(varNames(ne)) > 0){
            rp  <- getgrain(object, "rip")    
            host  <- get_superset_list(varNames(ne), rp$cliques)
            object$potential$pot_temp <- insertEvi(ne, getgrain(object, "pot_temp"), host)
            
            te <- if (is.null_evi(oe))
                      ne
                  else
                      union_evi(oe, ne)
            object$evidence <- te

        }         
    } 
    if (propagate) propagate(object) else object
}


#' @rdname grain_evidence
#' @export 
retractEvidence <- function(object, nodes=NULL, propagate=TRUE){
    stopifnot_grain(object)
    retract_evidence_worker(object, items=nodes, propagate=propagate)
}


## #' @rdname grain-evi
retract_evidence_worker <- function(object, items=NULL, propagate=TRUE){
    ##cat("++++ retractEvidence_\n")
    .resetgrain <- function(x){
        x$potential$pot_temp <- getgrain(x, "pot_orig")
        x$evidence       <- NULL
        isPropagated(x)  <- FALSE
        x
    }

    if ( is.null( items ) ){
        object <- .resetgrain( object )
    } else {
        if (!(is.character( items ) || is.numeric( items )) )
            stop("'items' must be a character or numeric vector")
        oe <- getEvidence( object )
        if (!is.null_evi( oe )){
            vn <- varNames( oe )
            if (is.numeric(items)){
                items <- vn[items]
            }
            keep <- setdiff( vn, items)
            ## NB: keep <= varNames            
            ## hvis keep==varNames(oe) så gør intet
            ## hvis keep=Ø så bare reset
            ## hvis Ø < keep < varNames så gør som nedenfor
            if ( length( keep ) < length( vn ) ){
                object <- .resetgrain( object )                
                if (length( keep ) > 0){
                    ne <- subset( oe, select=keep )
                    object <- set_evidence_worker( object, evidence = ne )
                }
            }
        }
    }
    if (propagate) propagate(object) else object
}


#' @export 
#' @rdname grain_evidence
absorbEvidence <- function(object, propagate=TRUE ){
    stopifnot_grain(object)
    absorbEvi_( object, propagate=propagate )
}

## #' @name grain-evi
absorbEvi_<- function(object, propagate=TRUE ){
    stopifnot_grain(object)

    ## Update 'object' as
    ## 1) set pot_orig <- pot_temp
    ## 2) ignore any finding
    object$potential$pot_orig <-  ## FRAGILE assignment
        getgrain(object, "pot_temp")
    object$evidence <- NULL
    isPropagated(object) <- FALSE

    if (propagate) propagate(object) else object
}


#' @export 
#' @name grain_evidence
pEvidence <- function(object, evidence=NULL){
    stopifnot_grain(object)
    
    has_evidence <- function(x) {
        !is.null(x$evidence)
    }

    if (is.null(evidence)){
        if (!has_evidence(object))
            return(NULL)
        else {
            attr(getgrain(object, "pot_equi"), "pEvidence")
        }
    } else {
        if (has_evidence(object)){
            stop("argument 'evidence' only allowed on networks without existing evidence\n")
        } else {
            compute_p_evidence(setEvidence(object, evidence=evidence,
                                           propagate=FALSE))
        }
    }
    
}

#' @export
#' @param short If TRUE a dataframe with a summary is returned;
#'     otherwise a list with all details.
#' @name grain_evidence
getEvidence <- function(object, short=TRUE){
    stopifnot_grain(object)    
    ev <- object$evidence
    if (is.null(ev)){
        return(NULL)
    }
    ev
}





## ## #' @name grain-evi
## addEvi  <- set_evidence_worker

#' @rdname grain_evidence
#' @export
evidence_add  <- set_evidence_worker

#' @rdname grain_evidence
evidence_get  <- getEvidence

#' @rdname grain_evidence
evidence_drop  <- retractEvidence

#' @rdname grain_evidence
evidence_p  <- pEvidence


## ## #' @name grain-evi
## dropEvi <- retractEvi



## #############################################################
##
## UTILITIES
##
## #############################################################

nodes_states_to_evidence <- function(nodes, states){
    if (!is.null( states ) && length( nodes ) == length( states )){
        evidence <- as.vector(states, "list")
        names(evidence) <- nodes
        return(evidence)
    } else {
        stop( "Either states are not given or nodes and states do not have same length" )
    }
}



## ## #' @name grain-evi
## "dropEvi<-" <- function(object, value=NULL){
##     retractEvi(object, value)
## }

## ## #' @name grain-evi
## "addEvi<-" <- function(object, value=NULL) {
##     set_evidence_worker(object, value)
## }

## #' @export 
## #' @name grain_evidence
## evidence <- function(object, short=TRUE) {
##     UseMethod("evidence")
## }

## #' @export 
## #' @name grain_evidence
## evidence.grain <- function(object, short=TRUE) {
##     getEvidence(object, short)
## }

## #' @export 
## #' @name grain_evidence
## #' @param value The evidence in the form of a named list or an evidence-object. 
## "evidence<-" <- function(object, value=NULL){
##     UseMethod("evidence<-")
## }

## #' @export 
## #' @name grain_evidence
## "evidence<-.grain" <- function(object, value=NULL) {
##     if (is.null( value )) 
##         retractEvi( object )
##     else
##         setEvidence(object, evidence=value)
## }

## Alternative names

## retractEvi <- function(object, nodes=NULL, propagate=TRUE){
##     stopifnot_grain(object)
##     retract_evidence_worker(object, items=nodes, propagate=propagate)
## }
