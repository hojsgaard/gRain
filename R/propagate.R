#' @title Propagate a graphical independence network (a Bayesian network)
#' 
#' @description Propagation refers to calibrating the cliques of the
#'     junction tree so that the clique potentials are consistent on
#'     their intersections; refer to the reference below for details.
#'
#' @name grain_propagate
#' 
#' @details The \code{propagate} method invokes \code{propagateLS}
#'     which is a pure R implementation of the Lauritzen-Spiegelhalter
#'     algorithm. The c++ based version is several times faster than
#'     the purely R based version.
#'
#' 
#' @aliases propagate.grain propagateLS propagateLS__
#' @param object A grain object
#' @param details For debugging info
#' @param engine Either "R" or "cpp"; "cpp" is the default and the
#'     fastest.
#' @param ... Currently not used
#' @return A compiled and propagated grain object.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{grain}}, \code{\link[gRbase]{compile}}
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{https://www.jstatsoft.org/v46/i10/}.
#' @keywords utilities models
#'
#' @examples
#'
#' example("grain")
#'
#' ## Uncompiled and unpropageted network:
#' bn0  <- grain(chest_cpt, compile=FALSE)
#' bn0
#' ## Compiled but unpropageted network:
#' bn1  <- compile(bn0, propagate=FALSE)
#' ## Compiled and propagated network
#' bn2  <- propagate(bn1)
#' bn2
#' ## Default is that networks are compiled but not propagated at creation time:
#' bn3  <- grain(chest_cpt) 
#' bn3 

#' @rdname grain_propagate
#' @export 
propagate.grain <- function(object, details=object$details, engine="cpp", ...){

    t0 <- proc.time()
    engine <- match.arg(tolower(engine), c("r", "cpp"))
    
    propfun <- switch(engine,
                      "r"   = {propagateLS},
                      "cpp" = {propagateLS__})
    
    object$potential$pot_equi <- ## FRAGILE assignment
        propfun(getgrain(object, "pot_temp"), rip=rip(object))
    
    isPropagated(object) <- TRUE
    .timing(" Time: propagation:", object$control, t0)
    object
}






## Lauritzen Spiegelhalter propagation
##

## Don't remember the idea behind the 'initialize' argument; should always be true

#' @rdname grain_propagate
#' @param cq_pot_list List of clique potentials
#' @param rip A rip ordering
#' @param initialize Always true.
#' 
#' @export
propagateLS <- function(cq_pot_list, rip, initialize=TRUE, details=0){
    cat(".Propagating BN: [propagateLS]\n")
    ## .infoPrint(details, 1, cat(".Propagating BN: [propagateLS]\n"))
    
    cliq       <- rip$cliques
    seps       <- rip$separators
    pa         <- rip$parent
    child_list <- rip$childList
    ncliq      <- length(cliq)

    ## Needed because RIP returns 0 as parent index for clique 1:
    pa[pa == 0] <- NA
    
    ## Backward propagation (collectEvidence) towards root of junction tree
    ##
    t0 <- proc.time()
    if (ncliq > 1){
        for (i in ncliq:2){
            sp   <- seps[[ i ]]
            cq.pot   <- cq_pot_list[[ i ]]
            pa.pot   <- cq_pot_list[[pa[ i ]]]

            if (length(sp) > 0){
                sp.pot               <- tableMargin(cq.pot, sp)
                cq_pot_list[[ i ]]     <- tableOp2(cq.pot, sp.pot, `/`)
                cq_pot_list[[pa[ i ]]] <- tableOp2(pa.pot, sp.pot, `*`)
            } else{
                zzz               <- sum(cq.pot)
                cq_pot_list[[1]]    <- cq_pot_list[[1]] * zzz
                cq_pot_list[[ i ]]  <- cq.pot / zzz
            }
        }
    }

    tmpd           <- cq_pot_list[[1]]
    norm_const     <- sum(tmpd)
    tmpd           <- tmpd / norm_const
    cq_pot_list[[1]] <- tmpd
        
    ## Forward propagation (distributeEvidence) away from root of junction tree
    ##

    t0 <- proc.time()
    for (i in 1:ncliq){
        ch <- child_list[[ i ]]
        if (length(ch) > 0)
        {
            for (j in 1:length(ch))
            {
                if (length(seps[[ch[ j ]]]) > 0)
                {
                    sp.pot            <- tableMargin(cq_pot_list[[ i ]], seps[[ch[ j ]]])
                    cq_pot_list[[ch[ j ]]]  <- tableOp2(cq_pot_list[[ch[ j ]]], sp.pot, `*`)
                }
            }
        }
    }
    
    attr(cq_pot_list, "pEvidence") <- norm_const
    cq_pot_list
}

#' @rdname grain_propagate
#' @export 
compute_p_evidence <- function(object, details=object$details, engine="cpp", ...){

    engine="r" ## FIXME
    t0 <- proc.time()
    engine <- match.arg(tolower(engine), c("r", "cpp"))
    
    compute_fun <- switch(engine,
                      "r"   = {compute_p_evidence_worker},
                      "cpp" = {compute_p_evidence_worker})
    
    out <- compute_fun(getgrain(object, "pot_temp"), rip=rip(object))
    .timing(" Time: propagation:", object$control, t0)
    out
}


## #' @export
compute_p_evidence_worker <- function(cq_pot_list, rip, initialize=TRUE, details=0){
    ## cat("[compute_p_evidence_worker]\n")
    
    cliq       <- rip$cliques
    seps       <- rip$separators
    pa         <- rip$parent
    child_list <- rip$childList
    ncliq      <- length(cliq)

    ## Needed because RIP returns 0 as parent index for clique 1:
    pa[pa == 0] <- NA
    
    ## Backward propagation (collectEvidence) towards root of junction tree
    ##
    t0 <- proc.time()
    if (ncliq > 1){
        for (i in ncliq:2){
            sp   <- seps[[ i ]]
            cq.pot   <- cq_pot_list[[ i ]]
            pa.pot   <- cq_pot_list[[pa[ i ]]]

            if (length(sp) > 0){
                sp.pot               <- tableMargin(cq.pot, sp)
                cq_pot_list[[ i ]]     <- tableOp2(cq.pot, sp.pot, `/`)
                cq_pot_list[[pa[ i ]]] <- tableOp2(pa.pot, sp.pot, `*`)
            } else{
                zzz               <- sum(cq.pot)
                cq_pot_list[[1]]    <- cq_pot_list[[1]] * zzz
                cq_pot_list[[ i ]]  <- cq.pot / zzz
            }
        }
    }

    tmpd           <- cq_pot_list[[1]]
    norm_const     <- sum(tmpd)
    return(norm_const)
}






    ## .infoPrint(details, 2, cat("..BACKWARD:\n"))

## cq   <- cliq[[ i ]]

    ## .infoPrint(details, 2, cat("..FORWARD:\n"))
            ## .infoPrint2(details, 2, "Clique %d: {%s}\n",  i , .colstr( cq ))
            ## .infoPrint2(details, 2, "Marg onto sep {%s}\n", .colstr(sp))
        ## .infoPrint2(details, 2, "Clique %d: {%s}\n",  i , .colstr(cliq[[ i ]]))
            ## .infoPrint2(details,2, "..Children: %s\n", .colstr(ch))
                    ## .infoPrint2(details, 2, "Marg onto sep %i: {%s}\n", ch[ j ], .colstr(seps[[ch[ j ]]]))
                    ## .infoPrint(details, 4, { cat("Marginal:\n"); print (sp.pot) })


                    ##cat(sprintf("......is.na: sp.pot=%i\n", any(is.na(sp.pot))))

            ## str(list(ncliq=ncliq, sp=sp, cq=cq))
            ## if ((length(sp) >= 1) && !is.na(sp)){ ## Changed on May 9, 2022
    
    ## FIXME: propagate.grain : Looks strange
    ## if (!is.null((ev <- getEvidence(object)))){
        ## attr(ev, "pEvidence") <- pEvidence(object)
        ## object$evidence <- ev
    ## }
    
