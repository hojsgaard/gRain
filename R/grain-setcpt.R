#' @title Update components of Bayesian network
#'
#' @description Update components of Bayesian network.
#'
#' @name cpt-update
#' 
#' @param object A grain object
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
#' yn   <- c("yes","no")
#' a    <- cptable(~asia,        values=c(1,99), levels=yn)
#' t.a  <- cptable(~tub + asia,  values=c(5,95,1,99), levels=yn)
#'
#' plist <- compileCPT(list(a, t.a )) 
#' bn    <- grain(plist)
#' bnc   <- compile(bn, propagate=FALSE)
#' bncp  <- compile(bn, propagate=TRUE)
#' 
#' ## New p(tub | asia)
#' z <- c(20, 80, 1, 99) 
#'
#' bn2   <- setCPT(bn, list(tub=z))
#' bnc2   <- setCPT(bnc, list(tub=z))
#' bncp2   <- setCPT(bncp, list(tub=z))
#' 

#' @rdname cpt-update
"setcpt<-" <- function(object, value){
    UseMethod("setcpt<-")
}

#' @rdname cpt-update
"setcpt<-.grain" <- function(object, value){
    setCPT(object, value)
}

#' @rdname cpt-update
setCPT <- function(object, value){
    UseMethod("setCPT")
}

.is.named.list <- function(x){
    inherits(x, "list") && !is.null(names(x))
}

#' @rdname cpt-update
setCPT.grain <- function(object, value){
    if (!.is.named.list(value))
        stop("'value' must be a named list")

    vn <- names(cpt(object))
    nn <- names(value)
    if (any((id <- is.na(match(nn, vn)))))
        stop("variable(s) ", toString(nn[id]), " not in network")   
    for (i in seq_along(nn)){
        v <- nn[i]
        z <- value[[i]]
        if (length(z) != length(cpt(object)[[v]]))
            stop("replacement value not correct length")                    
        object$cptlist[[v]][] <- z               
    }

    object$isCompiled <- object$isPropagated <- FALSE
    object
    
    ## if (!getgin(object, "isCompiled")){
    ##     ##cat("object is not compiled\n")
    ##     object
    ## } else {
    ##     isp <- getgin(object, "isPropagated")
    ##     ##cat("object IS compiled; propated?", isp, "\n")
    ##     ## RIP exists; just need to update potentials
    ##     object <- update_pot(object)
        
    ##     if (isp)
    ##         object <- propagate(object)
    ##     object
    ## } 
}
