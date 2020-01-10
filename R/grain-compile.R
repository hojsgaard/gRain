## NOTICE: the compiled object will contain a dag and a cptlist.
## These are not used for any calculations; only used for saving
## the network in Hugin format...

#' @title Compile a graphical independence network (a Bayesian network)
#' 
#' @description Compiles a Bayesian network. This means creating a
#'     junction tree and establishing clique potentials.
#'
#' @name grain_compile
#' 
#' @aliases compile.grain compile.cpt_grain compile.pot_grain
#' @param object A grain object.
#' @param propagate If TRUE the network is also propagated meaning
#'     that the cliques of the junction tree are calibrated to each
#'     other.
#' @param root A set of variables which must be in the root of the
#'     junction tree
#' @param control Controlling the compilation process.
#' @param details For debugging info. Do not use.
#' @param \dots Currently not used.
#' @return A compiled Bayesian network; an object of class
#'     \code{grain}.
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
#' ## To be written
#' 



#' @rdname grain_compile
add_jtree <- function(object, root=NULL){
    UseMethod("add_jtree")
}

#' @rdname grain_compile
add_jtree.cpt_grain <- function(object, root=NULL){
    if (!inherits(object, "cpt_grain")) stop("Not a cpt_grain object\n") 
    object[c("rip", "ug")] <- .create_jtree(object, root) 
    object
}
#' @rdname grain_compile
add_jtree.pot_grain <- function(object, root=NULL){
    if (is.null(rip(object)))
        stop("No rip component in object \n")
    object
}


.create_jtree <- function(object, root=NULL, update=TRUE){

    mdag <- moralize(getgin(object, "dag"), result="dgCMatrix")
    if (length(root) > 1) mdag <- .setRoot(mdag, root)
        
    ug_   <- triangulateMAT(mdag)  ## FIXME : Coercions are a MESS
    rp_   <- ripMAT( ug_ )
    ug_   <- as(ug_, "graphNEL")

    list(rip=rp_, ug=ug_)
}


#' @rdname grain_compile
add_potential <- function(object){
    UseMethod("add_potential")
}

#' @rdname grain_compile
add_potential.cpt_grain <- function(object){
    object$potential <- .create_potential(object)
    object
}

#' @rdname grain_compile
add_potential.pot_grain <- function(object){
    if (is.null(object$cqpot))
        stop("No cqpot component in object \n")
    object$potential <-
        list(pot_orig=object$cqpot,
             pot_temp=object$cqpot,
             pot_equi=.initArrayList(object$cqpot, NA))
    ##class(pot_orig) <- class(pot_temp) <- class(pot_equi) <- "pot_spec"
    object
}


.create_potential <- function(object){

    if (is.null(rip(object)))
        stop("No rip slot (junction tree) in object\n")
    if (is.null(cpt(object)))
        stop("No cpt slot in object\n")
    
    pot.1    <- .mkArrayList(rip(object), universe(object))
    pot_orig <- pot_temp <- .insertCPT(cpt(object), pot.1, details=0)
    pot_equi <- .initArrayList(pot.1, NA)

    ##class(pot_orig) <- class(pot_temp) <- class(pot_equi) <- "pot_spec"
    list(pot_orig=pot_orig, pot_temp=pot_temp, pot_equi=pot_equi)
}

#' @rdname grain_compile
compile.grain <- function(object, propagate=FALSE, root=NULL,
           control=object$control, details=0, ...) {    
    NextMethod("compile")
}


#' @rdname grain_compile
compile.cpt_grain <- function(object, propagate=FALSE, root=NULL, control=object$control,
                              details=0, ...){
        
        object <- add_jtree(object, root)
        object <- add_potential(object)
        
        object$isCompiled   <- TRUE
        object$isPropagated <- FALSE
        object$control      <- control
        
        if (propagate) propagate(object) else object
    }


#' @rdname grain_compile
compile.pot_grain <-
    function(object, propagate=FALSE, root=NULL, control=object$control,
             details=0,...) {
      
        object <- add_jtree(object, root)
        object <- add_potential(object)
        
        object$isCompiled   <- TRUE
        object$isPropagated <- FALSE
        object$control      <- control
        
        if (propagate) propagate(object) else object
    }


.createJTreeGraph <- function(rip){
    if (length(rip$cliques) > 1){
        ft <- cbind(rip$parents, 1:length(rip$parents))
        ft <- ft[ft[, 1] != 0, ,drop=FALSE]
        V <- seq_along(rip$parents)
        if (nrow(ft) == 0){
            jt <- new("graphNEL", nodes = as.character(V), edgemode = "undirected")
        } else {
            jt <- graph::ftM2graphNEL(ft, V=as.character(V), edgemode="undirected")
        }
    } else {
        jt <- new("graphNEL", nodes = "1", edgemode = "undirected")
    }
    jt
}


.timing <- function(text, control, t0){
  if (!is.null(control$timing) && control$timing)
    cat(sprintf("%40s", text), proc.time()-t0,"\n")

}



## setRoot: Completes the variables in <root> in the graph,
## FIXME: setRoot: Assumes sparse matrix. A MESS
.setRoot <- function(mdagM, root){
    vn  <- colnames(mdagM)
    dn  <- dimnames(mdagM)
    ft  <- names2pairs(match(root, vn), sort=FALSE, result="matrix")
    ft  <- rbind(ft, ft[, 2:1, drop=FALSE])
    mdagM <- .sparse_setXtf1(mdagM, ft)
    dimnames(mdagM) <- dn
    mdagM
}






    ##gg <- as(getgin(object, "dag"), "dgCMatrix")
    ##dagList(vpar(getgin(object, "cptlist")), result="dgCMatrix")
    ##gg <- moralizeMAT(gg)
    ##mdagM <- moralizeMAT( as(object$dag, "dgCMatrix") )
    ## Force variables in <root> to be complete in mdagM
    ## if (length(root) > 1) mdagM <- .setRoot( mdagM, root )




