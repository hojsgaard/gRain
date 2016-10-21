### #####################################################
###
### Creating grain objects
###
### #####################################################

#' @title Graphical Independence Network
#' 
#' @description The 'grain' builds a graphical independence network.
#'
#' @name grain-main
#' 
#' @details If 'smooth' is non-zero then entries of 'values' which a
#'     zero are replaced by the value of 'smooth' - BEFORE any
#'     normalization takes place.
#' 
#' @aliases grain grain.CPTspec grain.POTspec grain.graphNEL
#'     grain.dModel plot.grain iplot.grain
#' @param x An argument to build an independence network
#'     from. Typically a list of conditional probability tables, a DAG
#'     or an undirected graph. In the two latter cases, data must also
#'     be provided.
#' @param data An optional data set (currently must be an array/table)
#' @param control A list defining controls, see 'details' below.
#' @param smooth A (usually small) number to add to the counts of a
#'     table if the grain is built from a graph plus a dataset.
#' @param details Debugging information.
#' @param ... Additional arguments, currently not used.
#' @return An object of class "grain"
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{cptable}}, \code{\link{compile.grain}},
#'     \code{\link{propagate.grain}}, \code{\link{setFinding}},
#'     \code{\link{setEvidence}}, \code{\link{getFinding}},
#'     \code{\link{pFinding}}, \code{\link{retractFinding}}
#'     %\code{\link[gRbase]{gmData}}
#' @references S<f8>ren H<f8>jsgaard (2012). Graphical Independence
#'     Networks with the gRain Package for R. Journal of Statistical
#'     Software, 46(10), 1-26.
#'     \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords models
#' @examples
#' 
#' ## Asia (chest clinic) example:
#' yn   <- c("yes","no")
#' a    <- cptable(~asia,              values=c(1,99), levels=yn)
#' t.a  <- cptable(~tub+asia,          values=c(5,95,1,99), levels=yn)
#' s    <- cptable(~smoke,             values=c(5,5), levels=yn)
#' l.s  <- cptable(~lung+smoke,        values=c(1,9,1,99), levels=yn)
#' b.s  <- cptable(~bronc+smoke,       values=c(6,4,3,7), levels=yn)
#' e.lt <- cptable(~either+lung+tub,   values=c(1,0,1,0,1,0,0,1), levels=yn)
#' x.e  <- cptable(~xray+either,       values=c(98,2,5,95), levels=yn)
#' d.be <- cptable(~dysp+bronc+either, values=c(9,1,7,3,8,2,1,9), levels=yn)
#' plist <- compileCPT(list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be))
#' bn    <- grain(plist)
#' bn
#' summary(bn)
#' plot(bn)
#' bnc <- compile(bn, propagate=TRUE)
#' 
#' ## If we want to query the joint distribution of the disease nodes,
#' ## computations can be speeded up by forcing these nodes to be in
#' ## the same clique of the junction tree:
#' 
#' bnc2 <- compile(bn, root=c("lung", "bronc", "tub"), propagate=TRUE)
#' 
#' system.time({
#'   for (i in 1:200)
#'     querygrain(bnc, nodes=c("lung","bronc", "tub"), type="joint")})
#' system.time({
#'   for (i in 1:200)
#'     querygrain(bnc2, nodes=c("lung","bronc", "tub"), type="joint")})
#'
#'
#' ## Simple example - one clique only in triangulated graph:
#' plist.s <- compileCPT( list(a, t.a) )
#' bn.s <- grain( plist.s )
#' querygrain( bn.s )
#' 
#' ## Simple example - disconnected network:
#' plist.d <- compileCPT( list(a, t.a, s) )
#' bn.d <- grain( plist.d )
#' querygrain( bn.d )
#' 
#' 
#' ## Create network from data and graph specification.
#' ## There are different ways:
#' data(HairEyeColor)
#' hec <- HairEyeColor
#' daG <- dag( ~Hair + Eye:Hair + Sex:Hair )
#' class( daG )
#' uG <- ug( ~Eye:Hair + Sex:Hair )
#' class( uG )
#' 
#' ## Create directly from dag:
#' b1  <- grain( daG, hec )
#' class( b1 )
#' 
#' ## Build model from undirected (decomposable) graph
#' b3  <- grain( uG, hec )
#' class( b3 )
#' 
#' @export grain
grain <- function(x, data=NULL, control=list(), smooth=0, details=0,...){
  UseMethod("grain")
}

## A list of cpt's
##
#' @rdname grain-main
grain.CPTspec <- function(x, data=NULL, control=list(), smooth=0, details=0,...){
  ##cat("grain.CPTspec\n")
  control  <- .setControl(control)
  ans  <- c(list(universe    = attr(x,"universe"),
                 data        = data,
                 dag         = attr(x,"dag"),    ## Needed to save network in Hugin format
                 cptlist     = c(x)              ## Needed to save network in Hugin format
                 ),
            .setExtraComponents(control, details))

  class(ans) <- c("CPTgrain","grain")
  return(ans)
}

#' @rdname grain-main
grain.POTspec <- function(x, data=NULL, control=list(), smooth=0, details=0,...){
  ## cat("grain.POTspec\n")
  control  <- .setControl(control)
  ans  <- c(list(universe    = attr(x, "universe"),
                 data        = data,
                 equipot     = c(x),
                 ug          = attr(x, "ug"),
                 rip         = attr(x, "rip"),
                 dag         = attr(x, "dag"),    ## Needed to save network in Hugin format
                 cptlist     = attr(x, "cptlist") ## Needed to save network in Hugin format
                 ),
            .setExtraComponents(control, details))
  class(ans) <- c("POTgrain","grain")
  ans
}

## A graph + data (wrappers for calling grain.POTspec and grain.CPTspec)
#' @rdname grain-main
grain.graphNEL <- function(x, data=NULL, control=list(), smooth=0, details=0,...){
    if (missing(data))
        stop("Data must be given to create grain from graph\n")
    if (!(is.array(data) || is.data.frame(data)))
        stop("Data must be an array or a dataframe\n")

  if (is.DAG(x)){
    ans <- grain(compileCPT(extractCPT(data, x, smooth=smooth)),
                 data=data, control=control, details=details)
  } else {
    if (is.TUG(x)){
      ans <- grain(compilePOT(extractPOT(data, x, smooth=smooth)),
                   data=data, control=control, details=details)
    } else {
      stop("graph 'x' is neither a directed acyclic graph or a triangulated undirected graph")
    }
  }
  return(ans)
}

#' @rdname grain-main
grain.dModel <- function(x, data=NULL, control=list(), smooth=0, details=0,...){

    if (!x$isDecomposable)
        stop("Model must be decompsable graphical model\n")
    g <- ugList( terms( x ) )
    if (is.null(data))
        data <- x$datainfo$data
    grain(g, data=data, smooth=smooth, details=details, ...)
}


## Printing grain
##
print.grain <- function(x,...){
    cat("Independence network: Compiled:", x$isCompiled,
        "Propagated:", x$isPropagated, "\n")
    cat("  Nodes:"); str(unname(nodeNames(x)))
    if ( !is.null(x$evidence) ){
        cat("  Evidence:\n");
        print(as.data.frame( getEvidence(x)$summary) )
        if (!is.null((p <- pEvidence(x))))
            cat(sprintf("  pEvidence: %f\n", p))
    }
    return(invisible(x))
}

## FIXME: this print.grain is for new type of evidence...
print.grain <- function(x,...){
    cat("Independence network: Compiled:", x$isCompiled,
        "Propagated:", x$isPropagated, "\n")
    cat("  Nodes:"); str(unname(nodeNames(x)))
    if ( !is.null((ev <- evidence(x))) ){
        cat("  Evidence:\n");
        ##½        print(as.data.frame( ev ) )
        print( ev )
        if (!is.null((p <- pEvidence(x))))
            cat(sprintf("  pEvidence: %f\n", p))
    }
    return(invisible(x))
}

#' @rdname grain-main
#' @param object Any R object.
is.grain <- function(object){
    "grain" %in% class(object)
}


.setControl <- function(control){
  con <- list(timing=0)
  con[(namc <- names(control))] <- control
  con
}

.setExtraComponents <- function(control, details){
  list(
      ## FIXME Do we isInitialized? I doubt! Status: Now removed!
      ## isInitialized = FALSE,
      isCompiled    = FALSE,
      isPropagated  = FALSE,
      evidence      = NULL,
      control       = .setControl(control),
      details       = details
      )
}





## NOTICE:
##
## extractPOT() generates
## {p(C1), p(R2|S1), ..., p(Rn|Sn)}
## so these are clique potentials but they are not equilibrated.
## extractPOT() also generates  {p(v|pa(v)} which is not needed except to be
## able to save a network as a hugin net file.
##
## extractCPT() generates
## {p(v|pa(v))}
##
## compilePOT() and compileCPT() only makes 'internal' computations/setups
## and do not fundamentally change the above.
##

