## ##############################################################
##
#' @title Update CPTs of Bayesian network
#'
#' @description Update CPTs of Bayesian network.
#'
#' @name update-cpt
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' 
## ##############################################################
#'
#' @param object A `grain` object.
#' @param value A named list, see examples below.
#'
#' @details Updates some CPTs in a network but the overhead in redoing
#'     the triangulation and other steps are avoided.
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
#' p.R <- cptable(~R, values=c(.2, .8), levels=yn)
#' p.S_R <- cptable(~S:R, values=c(.01, .99, .4, .6), levels=yn)
#' p.G_SR <- cptable(~G:S:R, values=c(.99, .01, .8, .2, .9, .1, 0, 1), levels=yn)
#' 
#' x <- compileCPT(p.R, p.S_R, p.G_SR)
#' x
#' wet.bn <- grain(x)
#' 
#' getgrain(wet.bn, "cpt")
#' getgrain(wet.bn, "cpt")$R
#' getgrain(wet.bn, "cpt")$S
#'
#' # Now update some cpt's
#' wet.bn2 <- updateCPT(wet.bn, list(R=c(.3, .7), S=c(.1, .9, .7, .3)))
#' 
#' getgrain(wet.bn2, "cpt")$R
#' getgrain(wet.bn2, "cpt")$S
#' 
#' @export 
#' @rdname update-cpt
updateCPT <- function(object, value){
    UseMethod("updateCPT")
}

## Modifies cptlist in object. 

#' @export 
#' @rdname update-cpt
updateCPT.cpt_grain <- function(object, value){

    if (!isCompiled(object))
        stop("grain object must be compiled")

    if (!is_named_list(value))
        stop("'value' must be a named list")

    ## cat(".. inserting values in cptlist\n")
    vn <- names(getgrain(object, "cpt"))
    nn <- names(value)
    if (any((id <- is.na(match(nn, vn)))))
        stop("variable(s) ", toString(nn[id]), " not in network")   
    for (i in seq_along(nn)){
        v <- nn[i]
        z <- value[[i]]
        if (length(z) != length(getgrain(object, "cpt")[[v]]))
            stop("replacement value not correct length")                    
        ccc <- object$cptlist[[v]]        
        ## cat("tab - before\n"); print(ccc)
        ccc[] <- z
        ccc <- tabNormalize(ccc, "first")
        ## cat("tab - after\n"); print(ccc)
        object$cptlist[[v]] <- ccc
        #object$cptlist[[v]][] <- z               
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


























## update-cpt
## "setcpt<-" <- function(object, value){
##     UseMethod("setcpt<-")
## }

## #' @rdname update-cpt
## "setcpt<-.grain" <- function(object, value){
##     updateCPT(object, value)
## }



