#' @title Update components of Bayesian network
#'
#' @description Update components of Bayesian network.
#'
#' @name set-slot
#' 
#' @param object A grain object
#' @param list Entries to be replaced.
#' @param value A named list
#' 
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{grain}}, \code{\link[gRbase]{propagate}},
#'     \code{\link[gRbase]{triangulate}}, \code{\link[gRbase]{rip}},
#'     \code{\link[gRbase]{junctionTree}}
#' @references Søren Højsgaard (2012). Graphical Independence
#'     Networks with the gRain Package for R. Journal of Statistical
#'     Software, 46(10), 1-26.
#'     \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords utilities models
#'
#' @examples
#'
#' yn = c("yes", "no")
#' universe = list(flu=yn, temp=yn, headache=yn)
#' p.f   = tab(~flu, levels=universe, values=c(.1, .9))
#' p.t_f = tab(~temp|flu, levels=universe, values=c(.9, .1, .01, .99))
#' p.h_t = tab(~headache|temp, levels=universe, values=c(.8, .1, .1, .9))
#'
#' pl = compileCPT(list(p.f, p.t_f, p.h_t))
#' pl
#' bn = compile(grain(pl))
#' bn = propagate(bn)
#' bn
#' cp = cpt(bn)
#' cp1 = cp[[1]]
#' cp1

#' @rdname set-slot
set_rip <- function(object, value){
    if (!inherits(object, "cpt_grain")) stop("not a cpt_grain object\n")
    object$rip <- value
    object$ug  <- ugList(rip$cliques)
    object <- add_potential(object)
    object
}

#' @rdname set-slot
set_ug <- function(object, value){
    if (!inherits(object, "cpt_grain")) stop("not a cpt_grain object\n")
    ## FIXME check på value
    object$rip <- rip(value)
    object$ug  <- value
    object <- add_potential(object)
    object
}

#' @rdname set-slot
set_cpt <- function(object, list, value){
    if (!inherits(object, "cpt_grain")) stop("not a cpt_grain object\n")
    if (length(list) == 1 && !is.list(value))
        value <- list(value)
    object$cptlist[list] <- value
    object$isCompiled <- object$isPropagated <- FALSE
    object
}



