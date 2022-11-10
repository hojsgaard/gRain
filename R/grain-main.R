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
#' @param x An argument to build an independence network
#'     from. Typically a list of conditional probability tables, a DAG
#'     or an undirected graph. In the two latter cases, data must also
#'     be provided.
#' @param data An optional data set (currently must be an array/table)
#' @param control A list defining controls, see 'details' below.
#' @param smooth A (usually small) number to add to the counts of a
#'     table if the grain is built from a graph plus a dataset.
#' @param compile Should network be compiled.
#' @param details Debugging information.
#' @param ... Additional arguments, currently not used.
#' @return An object of class "grain"
#'
#' @note A change from earlier versions of this package is that grain
#'     objects are now compiled upon creation.
#' @seealso \code{\link{cptable}}, \code{\link{compile.grain}},
#'     \code{\link{propagate.grain}}, \code{\link{setFinding}},
#'     \code{\link{setEvidence}}, \code{\link{getFinding}},
#'     \code{\link{pFinding}}, \code{\link{retractFinding}},
#'     \code{\link{extractCPT}}, \code{\link{extractPOT}},
#'     \code{\link{compileCPT}}, \code{\link{compilePOT}}
#' 
#' @references Søren Højsgaard (2012). Graphical Independence
#'     Networks with the gRain Package for R. Journal of Statistical
#'     Software, 46(10), 1-26.
#'     \url{https://www.jstatsoft.org/v46/i10/}.
#' @keywords models
#' 
#' @examples
#' 
#' ## Create network from conditional probability tables CPTs:
#' 
#' yn   <- c("yes", "no")
#' a    <- cpt(~asia,              values=c(1,99), levels=yn)
#' t.a  <- cpt(~tub+asia,          values=c(5,95,1,99), levels=yn)
#' s    <- cpt(~smoke,             values=c(5,5), levels=yn)
#' l.s  <- cpt(~lung+smoke,        values=c(1,9,1,99), levels=yn)
#' b.s  <- cpt(~bronc+smoke,       values=c(6,4,3,7), levels=yn)
#' e.lt <- cpt(~either+lung+tub,   values=c(1,0,1,0,1,0,0,1), levels=yn)
#' x.e  <- cpt(~xray+either,       values=c(98,2,5,95), levels=yn)
#' d.be <- cpt(~dysp+bronc+either, values=c(9,1,7,3,8,2,1,9), levels=yn)
#' cpt_list  <- list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be)
#' chest_cpt <- compileCPT(cpt_list)
#' ## Alternative: chest_cpt <- compileCPT(a, t.a, s, l.s, b.s, e.lt, x.e, d.be)
#' 
#' chest_bn  <- grain(chest_cpt)
#'
#' ## Create network from data and graph specification.
#'
#' data(lizard, package="gRbase")
#'
#' ## From a DAG: height <- species -> diam
#' daG <- dag(~species + height:species + diam:species)
#'
#' ## From an undirected graph UG : [height:species][diam:species]
#' uG  <- ug(~height:species + diam:species)
#' 
#' liz.ug   <- grain(uG, data=lizard)
#' liz.dag  <- grain(daG, data=lizard)

#' @rdname grain-main
#' @export 
grain <- function(x, ...){
  UseMethod("grain")
}

#' @rdname grain-main
#' @export
grain.cpt_spec <- function(x, control=list(), smooth=0, compile=TRUE, details=0, ...){
    ##cat("grain.cpt_spec\n")
    control  <- .setControl(control)
    out  <- c(list(universe    = attr(x, "universe"),
                   cptlist     = as_cpt_spec_simple(x), ## Strips unnecessary stuff                   
                   dag         = attr(x, "dag")), ## FIXME Needed to save network in Hugin format                   
              .setExtraComponents(control, details))
    class(out) <- c("cpt_grain", "grain")
    if (compile) compile(out) else out
}


## For backward compatibility with bnlearn who calls methods and not generic functions...
## #' @method grain CPTspec

#' @export grain.CPTspec
grain.CPTspec <- grain.cpt_spec

#' @rdname grain-main
#' @export
grain.CPTspec <- grain.cpt_spec

#' @export
#' @rdname grain-main
grain.pot_spec <- function(x, control=list(), smooth=0, compile=TRUE, details=0,...){
    
    control  <- .setControl(control)
    out  <- c(list(universe    = attr(x, "universe"),              
                   cqpot       = x, ## FIXME: was c(x)...                  
                   ug          = attr(x, "ug"),
                   rip         = attr(x, "rip")),
              .setExtraComponents(control, details))
    ## FIXME: Generate dag if does not exist??
    class(out) <- c("pot_grain", "grain")
    if (compile) compile(out) else out
}

## A graph + data (wrappers for calling grain.pot_spec and grain.cpt_spec)
#' @export
#' @rdname grain-main
grain.graphNEL <- function(x, control=list(), smooth=0, compile=TRUE, details=0, data=NULL, ...){
    ## cat("grain.graphNEL\n")
    if (is.null(data))
        stop("Data must be given to create grain from graph\n")
    if (!(is.named.array(data) || is.data.frame(data)))
        stop("Data must be an array or a dataframe\n")

    if (is_dag(x)){
        zz <- extractCPT(data, x, smooth=smooth)
        zz <- compileCPT(zz)
    } else if (is_tug(x)){
        zz <- extractPOT(data, x, smooth=smooth)
        zz <- compilePOT(zz)
    }
    else
        stop("graph 'x' is neither a directed acyclic graph or a triangulated undirected graph")

    grain(zz, data=data, control=control, compile=compile, details=details)
}

#' @export
#' @rdname grain-main
grain.dModel <- function(x, control=list(), smooth=0, compile=TRUE, details=0, data=NULL, ...){
    if (!x$isDecomposable)
        stop("Model must be decompsable\n")
    if (is.null(data)) ## FIXME grain.dModel: Need to check data
        data <- x$datainfo$data

    graph_ <- ugList(terms(x))
    grain(graph_, data=data, smooth=smooth, compile=compile, details=details, ...)
}

#' @export
## #' @rdname grain-main
grain.pot_representation <- function(x, ...){
    grain(compilePOT(x))
}

#' @export
## #' @rdname grain-main
grain.cpt_representation <- function(x, ...){
    grain(compileCPT(x))
}

#' @export
print.grain <- function(x, ...){

    hasEvidence <- !is.null(evidence(x))
    cat("Independence network: Compiled:", isCompiled(x),
        "Propagated:", isPropagated(x), "Evidence:", hasEvidence, "\n")
    ## cat("  Nodes:"); str(unname(nodeNames(x)))
    invisible(x)
}

#' @export
summary.grain <- function(object, type='std', ...){

    type <- match.arg(type, c("std", "cliques", "rip", "configurations"))
    
    cat("Independence network: Compiled:", isCompiled(object),
        "Propagated:", isPropagated(object), "\n")
    
    if (length(object$evidence)) getEvidence(object)
    
    cat(" Nodes :")
    utils::str(nodeNames(object)) ## $universe$nodes)
    
    if (isCompiled(object)){
        
        cq.length <- sapply(rip(object)$cliques, length)
        
        cat(sprintf(" Number of cliques:              %4d \n",  length(cq.length)))
        cat(sprintf(" Maximal clique size:            %4d \n",  max(cq.length)))
        cat(sprintf(" Maximal state space in cliques: %4d \n",
                    max(unlist(lapply(getgrain(object, "pot_equi"), length)) )))
        ##                    max(unlist(lapply(pot(object)$pot_equi, length)) )))
      
        if(length(e <- getEvidence(object))){
            print(e)
        }
        
        switch(type,
               "rip"={
                   cat("\nRIP ordering:\n")
                   print(rip(object))
               },
               "cliques"={
                   cat("\nCliques:\n")
                   .printList(rip(object)$cliques)
               },
               "configurations"={
                   cat("\nConfigurations:\n")
                   nc <- seq_along(rip(object)$cliques)
                   for (i in nc){
                       ##cat(" clique ", i, ":", length(getgrain(object, "pot_equi")[[i]]), "\n")
                       ##cat(" clique ", i, ":", length(pot(object)$pot_equi[[i]]), "\n")
                   }
               })
    }
    invisible(object)
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





## #' 
## #' plot(bn.uG)
## #' plot(bn.daG)
## #' 
## #' querygrain(bn.uG)
## #' querygrain(bn.daG)
## #' 
## #' # Sanity: Are the distributions identical?
## #' t1 <- querygrain(bn.uG, type="joint")
## #' t2 <- querygrain(bn.daG, type="joint")
## #' t1 %a/% t2
## #' 
## #' # At a lower level
## #' bn2.uG <- extractPOT(lizard, ~height:species + diam:species)  %>% grain
## #' bn2.daG <- extractCPT(lizard, ~species + height:species + diam:species)  %>% grain
## #' 
## #' plot(bn2.uG)
## #' plot(bn2.daG)
## #' 
## #' t1 <- querygrain(bn2.uG, type="joint")
## #' t2 <- querygrain(bn2.daG, type="joint")
## #' t1 %a/% t2
## #' 




## #' @rdname grain-main    
## grain.marg_rep <- function(x, ...){grain(compileCPT(x))} FIXME to implement


## FIXME Rethink print.grain
