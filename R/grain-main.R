## cpt_spec
## - universe
## - cptlist
## - dag
##
## pot_spec
## - universe
## - cqpot (er det egentlig ikke klike marginaler? NEJ, det er betingede fordelinger)
## - ug
## - rip
##

## compile
## -------
##   cpt_spec
##   - rip
##   - ug
##   - potlist
##
##   pot_spec
##   - potlist

## Efter compilering:
## ------------------
## cptspec findes slet ikke, men selve cpt'erne er i et slot for sig selv (cpt)
## - cptlist (kopi)
## - dag (kopi)
##
## - universe (kopi)
## - ug (compile result)
## - rip (compile result)
## - potlist (compile result)
##
## cqpot findes slet ikke, men selve pot'erne er i et slot for sig selv (cqpot)
## - cqpot  (kopi)
##
## - universe  (kopi)
## - ug  (kopi)
## - rip  (kopi)
## - potlist (compile result)
##
## Der skal være metoder
##   get_dag, der tager en hvis den findes og ellers laver den
##   get_cqpot (get_cqmarg) der tager en hvis den findes og ellers laver den
##
## Hvis man ændrer cpt/pot:
## -- indenfor univers er det ok
## -- man må ændre værdi af cpt/pot men ikke ændre retning på pile / ikke ændre domæne
## -- man må ikke slette cpt/pot, men man må gerne give det uniform fordeling
##
## Hvis man ændrer ug:
## -- for cpt_spec: efterfølgende laves der: rip/potlist
## -- for pot_spec: ikke lovligt
##




## ###################################################################
##
#' @title Graphical Independence Network
#' @description Creating grain objects (graphical independence network).
#' @name grain-main
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
##
## ###################################################################

#' @details If 'smooth' is non-zero then entries of 'values' which a
#'     zero are replaced by the value of 'smooth' - BEFORE any
#'     normalization takes place.
#' 
#' @aliases grain grain.cpt_spec grain.pot_spec grain.graphNEL
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
#' @seealso \code{\link{cptable}}, \code{\link{compile.grain}},
#'     \code{\link{propagate.grain}}, \code{\link{setFinding}},
#'     \code{\link{setEvidence}}, \code{\link{getFinding}},
#'     \code{\link{pFinding}}, \code{\link{retractFinding}}
#' @references Søren Højsgaard (2012). Graphical Independence
#'     Networks with the gRain Package for R. Journal of Statistical
#'     Software, 46(10), 1-26.
#'     \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords models
#' 
#' @examples
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
#' \dontrun{
#' if (require(microbenchmark)){
#' microbenchmark(
#'   querygrain(bnc, nodes=c("lung","bronc", "tub"), type="joint"),
#'   querygrain(bnc2, nodes=c("lung","bronc", "tub"), type="joint")
#' )}}
#' 
#'
#' ## Simple example - one clique only in triangulated graph:
#' plist.s <- compileCPT(list(a, t.a))
#' bn.s <- grain(plist.s)
#' querygrain(bn.s)
#' 
#' ## Simple example - disconnected network:
#' plist.d <- compileCPT(list(a, t.a, s))
#' bn.d <- grain(plist.d)
#' querygrain(bn.d)
#' 
#' ## Create network from data and graph specification.
#' ## There are different ways:
#' 
#' data(HairEyeColor)
#' hec <- HairEyeColor
#' daG <- dag(~Hair + Eye:Hair + Sex:Hair)
#' class(daG)
#' uG <- ug( ~Eye:Hair + Sex:Hair)
#' class(uG)
#' 
#' ## Create directly from dag:
#' bn.dag  <- grain(daG, data=hec)
#' class(bn.dag)
#' compile(bn.dag)
#' 
#' ## Build model from undirected (decomposable) graph
#' bn.ug  <- grain(uG, data=hec)
#' class(bn.ug)
#' compile(bn.ug)
#' 
#' @export grain
grain <- function(x, control=list(), smooth=0, details=0, data=NULL, ...){
  UseMethod("grain")
}


#' @rdname grain-main
grain.cpt_spec <- function(x, control=list(), smooth=0, details=0, ...){

        
    ##cat("grain.cpt_spec\n")
    control  <- .setControl(control)
    out  <- c(list(universe    = attr(x, "universe"),
                   cptlist     = as_cpt_spec_simple(x), ## Strips unnecessary stuff                   
                   dag         = attr(x, "dag") ## FIXME Needed to save network in Hugin format                   
                   ),
              .setExtraComponents(control, details))
    ## FIXME: Generate dag if does not exist??
    class(out) <- c("cpt_grain", "grain")
    out
}

#' @rdname grain-main
grain.pot_spec <- function(x, control=list(), smooth=0, details=0,...){
    
    control  <- .setControl(control)
    out  <- c(list(universe    = attr(x, "universe"),              
                   cqpot       = x, ## FIXME: was c(x)...                  
                   ug          = attr(x, "ug"),
                   rip         = attr(x, "rip")
                   ),
              .setExtraComponents(control, details))
    ## FIXME: Generate dag if does not exist??
    class(out) <- c("pot_grain", "grain")
    out
}

#' @rdname grain-main
grain.pot_rep <- function(x, ...){grain(compile(x))}

#' @rdname grain-main
grain.cpt_rep <- function(x, ...){grain(compile(x))}


## A graph + data (wrappers for calling grain.pot_spec and grain.cpt_spec)
#' @rdname grain-main
grain.graphNEL <- function(x, control=list(), smooth=0, details=0, data=NULL, ...){
    if (is.null(data))
        stop("Data must be given to create grain from graph\n")
    if (!(is.named.array(data) || is.data.frame(data)))
        stop("Data must be an array or a dataframe\n")

    if (is_dag(x))
        zz <- extract_cpt(data, x, smooth=smooth)
    else if (is_tug(x))
        zz <- extract_pot(data, x, smooth=smooth)
    else
        stop("graph 'x' is neither a directed acyclic graph or a triangulated undirected graph")

    ## zz is either cpt_spec or pot_spec
    zz <- compile(zz)
    grain(zz, data=data, control=control, details=details)
}

#' @rdname grain-main
grain.dModel <- function(x, control=list(), smooth=0, details=0, data=NULL, ...){
    if (!x$isDecomposable)
        stop("Model must be decompsable\n")
    if (is.null(data)) ## FIXME grain.dModel: Need to check data
        data <- x$datainfo$data

    gg <- ugList(terms(x))
    grain(gg, data=data, smooth=smooth, details=details, ...)
}


## Printing grain
##
print.grain <- function(x, ...){
    cat("Independence network: Compiled:", .isComp(x),
        "Propagated:", .isProp(x), "\n")
    cat("  Nodes:"); str(unname(nodeNames(x)))
    if ( !is.null(x$evidence) ){
        cat("  Evidence:\n");
        print(as.data.frame( getEvidence(x)$summary) )
        if (!is.null((p <- pEvidence(x))))
            cat(sprintf("  pEvidence: %f\n", p))
    }
    invisible(x)
}

## FIXME: this print.grain is for new type of evidence...
print.grain <- function(x,...){
    cat("Independence network: Compiled:", .isComp(x),
        "Propagated:", .isProp(x), "\n")
    cat("  Nodes:"); str(unname(nodeNames(x)))
    if ( !is.null((ev <- evidence(x))) ){
        cat("  Evidence:\n");
        ##½        print(as.data.frame( ev ) )
        print( ev )
        if (!is.null((p <- pEvidence(x))))
            cat(sprintf("  pEvidence: %f\n", p))
    }
    invisible(x)
}


.setExtraComponents <- function(control, details){
  list(
      isCompiled    = FALSE,
      isPropagated  = FALSE,
      evidence      = NULL,
      control       = .setControl(control),
      details       = details
      )
}

.setControl <- function(control){
  con <- list(timing = 0)
  con[(namc <- names(control))] <- control
  con
}



