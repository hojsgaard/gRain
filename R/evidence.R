#' @title Set, update and remove evidence.
#' @description Set, update and remove evidence.
#' @name grain_evidence
#' @aliases evidence_get
#' 
#' @param object A "grain" object
#' @param evidence A list of name=value. See examples below.
#' @param nodes A vector of nodes. 
#' @param propagate Should the network be propagated?
#' @param details Debugging information
#' @param short If TRUE a dataframe with a summary is returned;
#'     otherwise a list with all details.
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
#' @examples
#' 
#' example("grain")
#' chest_bn <- grain(compileCPT(chest_cpt))
#'
#' bn2 <- chest_bn |> evidence_add(list(asia="yes", xray="yes"))
#' bn3 <- chest_bn |> evidence_add(list(asia=c(0.8, 0.1), xray="yes"))
#' 
#' bn2 |> evidence_get()
#' bn3 |> evidence_get()
#'
#' bn2 |> evidence_prob()
#' bn3 |> evidence_prob()
#'
#' bn2 |> evidence_drop("xray")
#' bn3 |> evidence_drop("xray")
#'
#' bn2 |> evidence_drop("xray") |> evidence_get()
#' bn3 |> evidence_drop("xray") |> evidence_get()
#' 
#' 
#' ## For backward compatibility these functions are available now but
#' # may be deprecated later.
#' bb2 <- setEvidence(chest_bn, c("asia", "xray"), c("yes", "yes"))
#' bb3 <- setEvidence(chest_bn, c("asia", "xray"), list(c(0.8, 0.2), "yes"))
#' bb4 <- setFinding(chest_bn, c("asia", "xray"), c("yes", "yes"))
#'
#' bb2 |> getEvidence()
#' bb3 |> getEvidence()
#'
#' bb2 |> retractEvidence("xray")
#' bb3 |> retractEvidence("xray")
#'
#' bb2 |> pEvidence()
#' bb3 |> pEvidence()
#' 
#' bb2 |> retractEvidence("xray") |> getEvidence()
#' bb3 |> retractEvidence("xray") |> getEvidence()
NULL

### WORKER FUNCTIONS ###

p_evidence_worker <- function(object, evidence=NULL){
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

absorb_evidence_worker <- function(object, propagate=TRUE ){
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

set_evidence_worker <- function(object, evidence=NULL, propagate=TRUE, details=0) {
    ## details=10
    ## cat("++++ set_evidence_worker input evidence: \n"); str(evidence)

    stopifnot_grain(object)    
    
    if (is.null_evi( evidence )){
        cat("Nothing to do\n")
    } else {
        if (!isCompiled(object)){
            object <- compile(object)
        } 
        
        old_evi_object <- getEvidence( object ) # Der er noget med typen (old/new)
        ## cat("old_evi_object:\n"); print(old_evi_object)

        ## if (inherits(evidence, "data.frame")) {
        ##     if (nrow(evidence) !=1) 
        ##         stop("'evidence' is dataframe but must have exactly one row\n")
        ##     evidence <- lapply(evidence, as.character)         
        ## }

        new_evi_object <- grain_evidence_new(evidence, universe(object)$levels)
        
        ## str(list(evidence=evidence, new_evi_object=new_evi_object))
        
        if (details > 0){
            cat("old.evidence:\n"); print.default(old_evi_object)
            cat("new evidence:\n"); print.default(new_evi_object)
        }
        
        ## Hvis der er eksisterende evidens, så skal det tages ud af det nye
        if (!is.null_evi( old_evi_object ) ){
            new_evi_object <- grain_evidence_setdiff( new_evi_object, old_evi_object )
            if (details>0) {
                cat("new evidence - after modification:\n"); print( new_evi_object )
            }
        }

        if (length(grain_evidence_names(new_evi_object)) > 0){
            rp  <- getgrain(object, "rip")    
            host  <- get_superset_list(grain_evidence_names(new_evi_object), rp$cliques)
            object$potential$pot_temp <- 
              insertEvi(new_evi_object, getgrain(object, "pot_temp"), host)
            
            te <- if (is.null_evi(old_evi_object))
                      new_evi_object
                  else
                      grain_evidence_union(old_evi_object, new_evi_object)
            object$evidence <- te
        }         
    } 
    if (propagate) propagate(object) else object
}

insertEvi <- function(evi_object, pot, hostclique) {
  if ( !inherits(evi_object, "grain_evidence") )
    stop("'object' is not a 'grain_evidence' object")
  
  for (i in seq_along( evi_object$evi_weight) ){
    p <- evi_object$evi_weight[[ i ]]
    j <- hostclique[ i ]
    pot[[j]] <- tabMult(pot[[ j ]], p)
  }
  pot
}


retract_evidence_worker <- function(object, nodes=NULL, propagate=TRUE) {
    ##cat("++++ retractEvidence_\n")
    .resetgrain <- function(x){
        x$potential$pot_temp <- getgrain(x, "pot_orig")
        x$evidence       <- NULL
        isPropagated(x)  <- FALSE
        x
    }

    if ( is.null( nodes ) ){
        object <- .resetgrain( object )
    } else {
        if (!(is.character( nodes ) || is.numeric( nodes )) )
            stop("'nodes' must be a character or numeric vector")
        old_evi_object <- getEvidence( object )
        if (!is.null_evi( old_evi_object )){
            #vn <- varNames( old_evi_object )
            vn <- grain_evidence_names(old_evi_object)
            if (is.numeric(nodes)){
                nodes <- vn[nodes]
            }
            keep <- setdiff( vn, nodes)
            ## str(list(keep=keep, vn=vn, nodes=nodes))
            ## NB: keep <= varNames            
            ## hvis keep==varNames(old_evi_object) så gør intet
            ## hvis keep=Ø så bare reset
            ## hvis Ø < keep < varNames så gør som nedenfor
            if ( length( keep ) < length( vn ) ){
                object <- .resetgrain( object )                
                if (length( keep ) > 0){
#                    ne <- subset( old_evi_object, select=keep )
                    idx <-vn %in% keep
                    new_evi_object <- old_evi_object[idx,]
                    print(new_evi_object)
                    object <- set_evidence_worker( object, evidence = new_evi_object )
                }
            }
        }
    }
    if (propagate) propagate(object) else object
}

get_evidence_worker <- function(object, short=TRUE) {
    stopifnot_grain(object)    
    ev <- object$evidence
    if (is.null(ev)){
        return(NULL)
    }
    ev
}

### WORKER FUNCTIONS END ###


#' @rdname grain_evidence
#' @export
evidence_add  <- function(object, evidence, propagate=TRUE, details=0){
    set_evidence_worker(object, evidence, propagate=propagate, details=details)    
}

#' @rdname grain_evidence
#' @export
evidence_get  <- get_evidence_worker

#' @rdname grain_evidence
#' @export
evidence_drop  <- retract_evidence_worker

#' @rdname grain_evidence
#' @export
evidence_prob  <- p_evidence_worker


### OLD NAMES ###

#' @name old_grain_evidence
#' @concept old_names
#' @inherit grain_evidence
#' @param states A vector of states (of the nodes given by
#'     'nodes'). Now deprecated; use argument 'evidence' instead.
#' @export 
setEvidence <- function(object, nodes=NULL, states=NULL, evidence=NULL, 
                        propagate=TRUE, details=0) {
    stopifnot_grain(object)    
    
    if (is.null(evidence) && is.null(nodes))
        stop("Evidence is not given; nothing to do...")
    
    if (is.null(evidence)) ## Then 'nodes' must be given
        evidence <- nodes_states_to_evidence(nodes, states)   
    
    set_evidence_worker(object, evidence, propagate=propagate, details=details)    
}


#' @rdname old_grain_evidence
#' @export 
retractEvidence <- function(object, nodes=NULL, propagate=TRUE) {
    stopifnot_grain(object)
    retract_evidence_worker(object, nodes=nodes, propagate=propagate)
}

#' @rdname old_grain_evidence
#' @export 
absorbEvidence <- function(object, propagate=TRUE ) {
    stopifnot_grain(object)
    absorb_evidence_worker( object, propagate=propagate )
}

#' @rdname old_grain_evidence
#' @export
getEvidence <- function(object, short=TRUE) {
    get_evidence_worker(object, short=short)
}

#' @rdname old_grain_evidence
#' @export 
pEvidence <- function(object, evidence=NULL) {
    p_evidence_worker(object, evidence=evidence)
}



nodes_states_to_evidence <- function(nodes, states) {
    if (!is.null( states ) && length( nodes ) == length( states )){
        evidence <- as.vector(states, "list")
        names(evidence) <- nodes
        return(evidence)
    } else {
        stop( "Either states are not given or nodes and states do not have same length" )
    }
}
