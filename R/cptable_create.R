#' @title Create conditional probability tables (CPTs)
#' 
#' @description Creates conditional probability tables of the form
#'     p(v|pa(v)).
#'
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#'
#' @name cpt
#' 
#' @param names Specifications of the names in P(v|pa1,...pak). See
#'     section 'details' for information about the form of the
#'     argument.
#'
#' @param levels 1) a list with specification of the levels of the
#'     factors in \code{names} or 2) a vector with number of levels of
#'     the factors in \code{names}. See 'examples' below.
#' 
#' @param values Probabilities; recycled if necessary. Regarding the
#'     order, please see section 'details' and the examples.
#' 
#' @param normalize See 'details' below.
#'
#' @param smooth Should values be smoothed, see 'Details' below.
#'
#' @details
#'
#' `cptable` is simply a wrapper for `cpt` and the functions can hence
#' be used synonymously.
#'
#' If `smooth` is non--zero, then this value is added to all cells __before__
#' normalization takes place.
#' 
#' Regarding the form of the argument \code{names}: To specify
#' \eqn{P(a|b,c)} one may write `~a|b:c`, `~a:b:c`,
#' `~a|b+c`, `~a+b+c` or `c("a","b","c")`. Internally,
#' the last form is used. Notice that the \code{+} and \code{:}
#' operator are used as a separators only. The order of the variables IS
#' important so the operators DO NOT commute.
#'
#' The first variable in `levels` varies fastest. 
#'
#' @return An array.
#' @keywords utilities
#' 
#' @seealso \code{\link{andtable}}, \code{\link{ortable}},
#'     \code{\link{extract_cpt}}, \code{\link{compileCPT}},
#'     \code{\link{extract_cpt}}, \code{\link{compilePOT}},
#'     \code{\link{grain}}
#' @references Søren Højsgaard (2012). Graphical Independence
#'     Networks with the gRain Package for R. Journal of Statistical
#'     Software, 46(10), 1-26.
#'     \url{https://www.jstatsoft.org/v46/i10/}.
#' @keywords models
#' @examples
#'
#' ## See the wet grass example at
#' ## https://en.wikipedia.org/wiki/Bayesian_network
#' 
#' yn <- c("yes", "no")
#' ssp <- list(R=yn, S=yn, G=yn) # state space
#' 
#' ## Different forms
#' t1 <- cpt(c("S", "R"), levels=ssp,     values=c(.01, .99, .4, .6))
#' t2 <- cpt(~S:R,        levels=ssp,     values=c(.01, .99, .4, .6))
#' t3 <- cpt(~S:R,        levels=c(2, 2), values=c(.01, .99, .4, .6))
#' t4 <- cpt(~S:R,        levels=yn,      values=c(.01, .99, .4, .6))
#' t1; t2; t3; t4
#'
#' varNames(t1)
#' valueLabels(t1)
#' 
#' ## Wet grass example
#' ssp <- list(R=yn, S=yn, G=yn) # state space
#' p.R    <- cptable(~R,     levels=ssp, values=c(.2, .8))
#' p.S_R  <- cptable(~S:R,   levels=ssp, values=c(.01, .99, .4, .6))
#' p.G_SR <- cptable(~G:S:R, levels=ssp, values=c(.99, .01, .8, .2, .9, .1, 0, 1))
#'
#' wet.cpt <- compileCPT(p.R, p.S_R, p.G_SR)
#' wet.cpt
#' wet.cpt$S # etc
#'
#' # A Bayesian network is created with:
#' wet.bn <- grain(wet.cpt)
#' 
#' @rdname cpt
#' @export
cpt <- function(names, levels, values, normalize = "first", smooth=0){
    names <- c(.formula2char(names))
    ## cat("cpt................")
    ## str(list(names=names, levels=levels, values=values, normalize = normalize, smooth=smooth))
    tabNew(names=names, levels=levels, values=values, normalize = normalize, smooth=smooth)
}

## #' @rdname cpt
## #' @export
## cptable <- cpt




#' @rdname cpt
#' @param vpar node an its parents
#' @export
cptable <- function(vpar, levels=NULL, values=NULL, normalize=TRUE,  smooth=0 ){
    vpa  <- c(.formula2char(vpar))        
    if (is.list(levels)){
        v <- vpa[1]
        if (!(v %in% names(levels)))
            stop(paste0("Name ", v, " is not in the 'levels' list\n"))
        levels <- levels[[v]]
    }
    ##str(list(vpa=vpa, xlevels=levels))
    ## if (is.null(values))
    ##     values <- rep(1.0, length(levels))
    out  <- values
    attributes(out) <-
        list(vpa=vpa, normalize=normalize,
             smooth=smooth, levels=levels)
    class(out) <- "cptable_class"
    out
}



