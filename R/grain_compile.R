## NOTICE: the compiled object will contain a dag and a cptlist.
## These are not used for any calculations; only used for saving
## the network in Hugin format...

## FIXME compile examples to be written

#' @title Compile Bayesian network.
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
#' @param tug A triangulated undirected graph. 
#' @param root A set of variables which must be in the root of the
#'     junction tree
#' @param control Controlling the compilation process.
#' @param details For debugging info. Do not use.
#' @param \dots Currently not used.
#' @return A compiled Bayesian network; an object of class
#'     \code{grain}.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' 
#' @seealso \code{\link{grain}}, \code{\link[gRbase]{propagate}},
#'     \code{\link{propagate.grain}},
#'     \code{\link[gRbase]{triangulate}}, \code{\link[gRbase]{rip}},
#'     \code{\link[gRbase]{junctionTree}}
#'
#' @references Søren Højsgaard (2012). Graphical Independence
#'     Networks with the gRain Package for R. Journal of Statistical
#'     Software, 46(10), 1-26.
#'     \url{https://www.jstatsoft.org/v46/i10/}.
#' @keywords utilities models

#' @rdname grain_compile
#' @export
compile.grain <- function(object, propagate=FALSE, tug=NULL, root=NULL,
           control=object$control, details=0, ...) {    

    object <- .add_jtree(object, tug=tug, root=root)
    object <- .add_potential(object)
    
    isCompiled(object) <- TRUE
    isPropagated(object) <- FALSE
    
    object$control      <- control
    if (propagate) propagate(object) else object
}



## #############################################
##
## dot-functions only below here
##
## #############################################


.timing <- function(text, control, t0) {
  if (!is.null(control$timing) && control$timing)
    cat(sprintf("%40s", text), proc.time()-t0,"\n")

}

## Completes the variables in <set> in the graph,
.make_set_complete <- function(tugM, set) {
    vn   <- colnames(tugM)
    dn   <- dimnames(tugM)
    ft   <- names2pairs(match(set, vn), sort=FALSE, result="matrix")
    ft   <- rbind(ft, ft[, 2:1, drop=FALSE])
    tugM <- .sparse_setXtf1(tugM, ft)
    dimnames(tugM) <- dn
    tugM
}

.add_jtree <- function(object, tug=NULL, root=NULL) {
    UseMethod(".add_jtree")
}

## #' @rdname grain_compile
.add_jtree.cpt_grain <- function(object, tug=NULL,root=NULL) {

    if (is.null(tug)) {
        tug <- moralize(getgin(object, "dag"), result="dgCMatrix")        
    } 
    object[c("rip", "ug")] <- .create_jtree(tug, root) 
    object
}

## #' @rdname grain_compile
.add_jtree.pot_grain <- function(object, tug=NULL, root=NULL) {
    if (is.null(rip(object)))
        stop("No rip component in object \n")
    object
}

.create_jtree <- function(tugM, root=NULL, update=TRUE) {

    tugM <- as(tugM, "dgCMatrix")
    if (length(root) > 1)
        tugM <- .make_set_complete(tugM, root)
    ug_   <- triangulateMAT(tugM)  ## FIXME : Coercions are a MESS
    rp_   <- ripMAT(ug_)
    ug_   <- as(ug_, "igraph")
    
    list(rip=rp_, ug=ug_)
}



## #' @rdname grain_compile
.add_potential <- function(object){
    UseMethod(".add_potential")
}


.create_potential <- function(object) {

    pot.1    <- .make_array_list(getgrain(object, "rip"), universe(object))
    pot_orig <- pot_temp <- .insert_CPT(getgrain(object, "cpt"), pot.1, details=0)
    pot_equi <- .initialize_array_list(pot.1, values=NA)

    list(pot_orig=pot_orig, pot_temp=pot_temp, pot_equi=pot_equi)
}

## #' @rdname grain_compile
.add_potential.cpt_grain <- function(object) {
    object$potential <- .create_potential(object)
    object
}

## #' @rdname grain_compile
.add_potential.pot_grain <- function(object){
    if (is.null(object$cqpot))
        stop("No cqpot component in object \n")

    object$potential <-
        list(pot_orig=object$cqpot,
             pot_temp=object$cqpot,
             pot_equi=.initialize_array_list(object$cqpot, values=NA))
    object
}

## Create potential list (rip, universe)
##
.make_array_list <- function(rip.order, universe, values=1){
    cliques <- rip.order$cliques    
    potlist  <- as.list(rep(NA, length(cliques)))
    
    for ( i in seq_along(cliques)){
        cq    <- cliques[[ i ]]
        vlab  <- universe$levels[cq]
        potlist[[ i ]] <- tabNew(cq, vlab, values)
    }
    potlist
}

.initialize_array_list <- function(x, values=NA){
    lapply(x, function(z) {
        z[] <- values             
        z
    } )
}

## Insert cpt's into potential list (cptlist, APlist)
##
.insert_CPT <- function(cptlist, potlist, details=0) {
    if (details>=1)
        cat(".Inserting cpt's in potential list [.insert_CPT]\n")

    pot_names <- lapply(potlist, function(x) names(dimnames(x)))
    cpt_names <- unname(lapply(cptlist, function(x) varNames(x)))
    hosts    <-  get_superset_list(cpt_names, pot_names)

    ## str(list(cpt_names=cpt_names, pot_names=pot_names, hosts=hosts))
    
    for (i in 1:length(cptlist)) {
            cptc <- cptlist[[ i ]]
            h    <- hosts[ i ]
            ## str(list(i=i, h=h))
            ## print(h); print(potlist[[h]])
            potlist[[ h ]] <- tableOp( potlist[[ h ]], cptc, "*" )
        }
    .infoPrint(details, 4, {cat("....potlist (after insertion):\n"); print(potlist) })
    potlist
}

