## ##############################################################
##
#' @title Replace CPTs in Bayesian network
#'
#' @description Replace CPTs of Bayesian network.
#'
#' @name replace-cpt
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' 
## ##############################################################
#'
#' @param object A `grain` object.
#' @param value A named list, see examples below.
#'
#' @details When a Bayesian network (BN) is constructed from a list of
#'     conditional probability tables (CPTs) (e.g. using the function
#'     `grain()`), various actions are taken:
#'
#' 1. It is checked that the list of CPTs define a directed acyclic graph (DAG).
#'
#' 1. The DAG is moralized and triangulated.
#'
#' 1. A list of clique potentials (one for each clique in the
#'    triangulated graph) is created from the list of CPTs.
#'
#' 1. The clique potentials are, by default, calibrated to each other
#'     so that the potentials contain marginal distributions. 
#'
#' The function described here bypass the first two steps which can
#' provide an important gain in speed compared to constructing a new
#' BN with a new set of CPTs with the same DAG.
#' 
#' @seealso \code{\link{grain}}, \code{\link[gRbase]{propagate}},
#'     \code{\link[gRbase]{triangulate}}, \code{\link[gRbase]{rip}},
#'     \code{\link[gRbase]{junctionTree}}
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{https://www.jstatsoft.org/v46/i10/}.
#' 
#' @keywords utilities models
#' @examples
#' ## See the wet grass example at
#' ## https://en.wikipedia.org/wiki/Bayesian_network
#' 
#' yn <- c("yes", "no")
#' p.R    <- cptable(~R, values=c(.2, .8), levels=yn)
#' p.S_R  <- cptable(~S:R, values=c(.01, .99, .4, .6), levels=yn)
#' p.G_SR <- cptable(~G:S:R, values=c(.99, .01, .8, .2, .9, .1, 0, 1), levels=yn)
#' 
#' wet.bn <- compileCPT(p.R, p.S_R, p.G_SR)  |> grain()
#' getgrain(wet.bn, "cpt")[c("R","S")]
#' 
#' # Update some CPTs
#' wet.bn <- replaceCPT(wet.bn, list(R=c(.3, .7), S=c(.1, .9, .7, .3)))
#' getgrain(wet.bn, "cpt")[c("R","S")]
#' 
#' @export 
#' @rdname replace-cpt
replaceCPT <- function(object, value){
    UseMethod("replaceCPT")
}

## Modifies cptlist in object. 

#' @export 
#' @rdname replace-cpt
replaceCPT.cpt_grain <- function(object, value){

    if (!isCompiled(object))
        stop("grain object must be compiled")

    if (!is_named_list(value))
        stop("'value' must be a named list")

    ## cat(".. inserting values in cptlist\n")
    vn <- names(getgrain(object, "cpt"))
    nms <- names(value)
    if (any((id <- is.na(match(nms, vn)))))
        stop("variable(s) ", toString(nms[id]), " not in network")   
    for (i in seq_along(nms)){
        v <- nms[i]
        z <- value[[i]]
        if (length(z) != length(getgrain(object, "cpt")[[v]]))
            stop("replacement value not correct length")                    
        ccc <- object$cptlist[[v]]        
        ccc[] <- z
        ccc <- tabNormalize(ccc, "first")
        object$cptlist[[v]] <- ccc
    }
    ## isCompiled(object) <- FALSE

    ## FIXME: This can be optimized further for speed because only some potentials are typically modified
    pot.1 <- .initialize_array_list(getgrain(object, "pot_orig"), values=1)    
    object$potential$pot_orig <- object$potential$pot_temp <-
        .insert_CPT(getgrain(object, "cpt"), pot.1, details=0)
        
    isPropagated(object) <- FALSE
    object
}

is_named_list <- function(x){
    inherits(x, "list") && !is.null(names(x))
}
















