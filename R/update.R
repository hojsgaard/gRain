#' @title Update a Bayesian network
#' 
#' @description Update a Bayesian network
#' 
#' @param object A Bayesian network of class \code{cpt_grain}
#' @param \dots If \code{CPTlist} is a name in the dotted list, then
#'     the object will be update with this value (which is assumed to
#'     be a list of conditional probabilities).
#' 
#' @return A Bayesian network.
#' 
#' @note There is NO checking that the input matches the settings in
#'     the Bayesian network.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords utilities
#' @examples
#'
#' ## FIXME : Improve example; avoid compileCPT
#' 
#' ## Network for Bernulli experiment; two nodes: X and thetaX
#' yn  <- c("yes", "no")    # Values for X
#' thX.val <- c(.3, .5, .7) # Values for thetaX
#' prX.val <- rep(1, length(thX.val)) # Probabilities for thetaX values
#' 
#' thX <- cptable(~thetaX, values=prX.val, levels=thX.val)
#' X   <- cptable(~X|thetaX, values=rbind(thX.val,1-thX.val), levels=yn)
#' 
#' cptlist <- compileCPT( list(thX, X) )
#' bn  <- compile( grain( cptlist ) )
#' querygrain( setEvidence(bn, nodes="X", states="yes") )
#' 
#' ## To insert a new prior distribution we may do as follows
#' ## (where we can omit the process of recompiling the network)
#' prX.val2 <- c(.2,.3,.5)
#' thX2 <- cptable(~thetaX, values=prX.val2, levels=thX.val)
#' bn2 <- update(bn, CPTlist=compileCPT( list(thX2, X)))
#' querygrain( setEvidence(bn2, nodes="X", states="yes") )
#' 
#' @export update.cpt_grain
"update.cpt_grain" <- function(object,  ...){
    
    if (!.isComp(object))
        object <- compile( object )
    
    ##cl <- match.call(expand.dots=TRUE)
    args <- list(...)
    arg.names <- names(args)
    
    if ("CPTlist" %in% arg.names){
        object$cptlist[names(args$CPTlist)] <- args$CPTlist
        pot.with.1        <- .mkArrayList( rip(object), uni(object) )
        newpot            <- .insertCPT(object$cptlist, pot.with.1, details=0)
        object$potential <-
            list(pot_orig=newpot,
                 pot_temp=newpot,
                 pot_equi=.initArrayList(pot.with.1, NA))
        object$isPropagated <- FALSE
    }
    object
}





















