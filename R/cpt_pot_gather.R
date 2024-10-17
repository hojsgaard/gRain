
#' @title Compile conditional probability tables / cliques potentials.
#' @description Compile conditional probability tables / cliques
#'     potentials as a preprocessing step for creating a graphical
#'     independence network
#'
#' @name components_gather
#' 
#' @param x To \code{compileCPT} x is a list of conditional
#'     probability tables; to \code{compilePOT}, x is a list of clique
#'     potentials.
#'
## #' @param object A list of potentials or of CPTs.
#'
#' @param forceCheck Controls if consistency checks of the probability
#'     tables should be made.
#' 
#' @param ... Additional arguments; currently not used.
#' 
#' @aliases parse_cpt, parse_cpt.xtabs, parse_cpt.default
#' 
#' @details
#'     * `compileCPT` is relevant for turning a collection of
#'     cptable's into an object from which a network can be built. For
#'     example, when specification of a cpt is made with cptable then
#'     the levels of the node is given but not the levels of the
#'     parents. `compileCPT` checks that the levels of variables in
#'     the cpt's are consistent and also that the specifications
#'     define a dag.
#' 
#'     * `compilePOT` is not of direct relevance for the
#'     user for the moment. However, the elements of the input should
#'     be arrays which define a chordal undirected graph and the
#'     arrays should, if multiplied, form a valid probability density.
#'  
#' @return A list with a class attribute.
#' 
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#'
#' @seealso \code{\link{extract_cpt}}, \code{\link{extract_pot}}, \code{\link{extract_marg}}
#' 
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{https://www.jstatsoft.org/v46/i10/}.
#'
#' @keywords utilities
#'
#' @examples
#'
#' example("example_chest_cpt")
#' x <- compile_cpt(chest_cpt)
#' class(x)
#' grain(x)
#' 


## FIXME: This listify stuff caused trouble as the class attribute was lost. For potentials got rid of it all. 

#' @rdname components_gather
#' @export
compile_cpt <- function(x, ..., forceCheck=TRUE) {
    args <- c(list(x), list(...))
    args <- listify_dots(args)
    compile_cpt_worker(args, forceCheck=forceCheck)
    ## compile_cpt_worker(x, forceCheck=forceCheck)
}


#' @rdname components_gather
#' @export
compile_pot <- function(x, ..., forceCheck=TRUE){
  # doit(lapply(x, .namesDimnames))
  #   dots <- list(...)
  #   args <- c(list(x), dots)
  #   args <- listify_dots(args)
  #   args100 <<- args
    compile_pot_worker(x, forceCheck=forceCheck)   
}



## -------------------------------------------------------------
## For backward compatibility; deprecate in future release
## -------------------------------------------------------------

#' @rdname old_components_gather
#' @inherit components_gather
#' @concept old_names
#' @export
compileCPT <- compile_cpt

#' @rdname old_components_gather
#' @export
compilePOT <- compile_pot



## #############################################################

#' @export
print.cpt_spec <- function(x, ...){
    ## cat("cpt_spec with probabilities:\n")
    lapply(x,
           function(xx){
               vn <- varNames(xx)
               .print_probability(vn)
           })
    invisible(x)
}

#' @export
print.pot_spec <- function(x, ...){
    ## cat("pot_spec with potentials:\n")
    lapply(x,
           function(xx){
               vn <- names(dimnames(xx))
               cat("(", paste(vn, collapse=' '),") \n")
           })    
    invisible(x)
}

summary.cpt_spec <- function(object, ...){
    ## cat("cpt_spec with probabilities:\n")
    lapply(object,
           function(xx){
               vn <- varNames(xx)
               .print_probability(vn)               
           })
    invisible(object)

}

## ###########################################################
## Helper functions  -- used only in grain-main.R
## ###########################################################

as_cpt_spec_simple <- function(x){
    z <- c(x)
    attr(z, "universe") <- attr(x, "universe")
    class(z) <- "cpt_spec_simple"
    z
}

print.cpt_spec_simple <- function(x,...){
    cat("cpt_spec_simple with probabilities:\n")
    lapply(x,
           function(xx){
               vn <- varNames(xx)
               .print_probability(vn)               
           })
  invisible(x)
}



## ###################################################
##
## dot functions below here
##
## ###################################################

compile_cpt_worker <- function(x, forceCheck=TRUE) {
    ## x: A list of cpts (arrays)
    ## type <- is.list(x) ##&& all(sapply(x, is.named.array))
    ## if (!type) stop("A list of named arrays is expected")
        
    ## zz: Internal representation of cpts
    zz  <- lapply(x, parse_cpt)
    
    universe <- .create_universe(zz)

    ## Given node names; need to check that they are not replicated
    vn_given <- sapply(zz, "[[", "vnam")
    if (length(vn_given) != length(unique(vn_given)))
        stop("Some nodes specified more than once: ", toString(vn_given))
    
    ## Are all cpts defined?
    ss <- setdiff(unique(unlist(vn_given)),  universe$nodes)
    if (length(ss) > 0)
        stop(paste("Distribution not specified for nodes(s):", toString(ss)))


    ## Does specification define a DAG? If x is cpt_representation the answer is yes
    if (inherits(x, "cpt_representation")) {
        graph <- attr(x, "graph")
    } else {
        vp <- lapply(zz, "[[", "vpar")
        graph <- dagList(vp, forceCheck=forceCheck, result="igraph")
    }

    ## Need list of cpts (each represented as an array)
    out <- lapply(seq_along(zz), .create_array, zz, universe)    
    names(out) <- universe$nodes
    
    attr(out, "universe") <- universe
    attr(out, "dag")      <- graph
    class(out)            <- "cpt_spec"
    out
}

compile_pot_worker <- function(x, ...){
    ## x: a list of arrays, and a rip attribute

    type <- is.list(x) && all(sapply(x, is.named.array))
    if (!type) stop("A list of named arrays is expected")    

    universe  <- .make.universe(x)
    
    
    attr(x, "universe") <- universe
    attr(x, "ug")    <- attr(x, "graph") 
    attr(x, "graph") <- NULL
#    attr(x, "rip")   <- rp
    class(x) <- "pot_spec"
    x
}

# if (inherits(x, "pot_representation")){ 
#   cat("Result of extract_pot\n")
#     graph <- attr(x, "graph")
#     rp    <- attr(x, "rip")
# } else {
#   cat("Not Result of extract_pot\n")
#   graph <- lapply(x, .namesDimnames)
#   graph <- ug(graph)
#   rp    <- rip(graph)
# }


.create_universe <- function(zz){
    vn <- unlist(lapply(zz, "[[", "vnam"))
    vl <- lapply(zz, "[[", "vlev")
    di <- unlist(lapply(vl, length))
    names(vl) <- vn        
    universe  <- list(nodes = vn, levels = vl, nlev = di)
    universe
}

.create_array <- function(i, zz, universe){
    cp <- zz[[i]]
    dn <- universe$levels[cp$vpar]
    di <- sapply(dn, length)
    val <- array(rep(1.0, prod(di)), dim=di, dimnames=dn)
    if (length(cp$values) > 0)
        val[] <- cp$values + cp$smooth
    val
}

.make.universe <- function(x){
    lll       <- unlist(lapply(x, dimnames), recursive=FALSE)
    nnn       <- names(lll)
    iii       <- match(unique(nnn), nnn)
    levels    <- lll[iii]
    vn        <- nnn[iii]
    di        <- c(lapply(levels, length), recursive=TRUE)
    names(di) <- vn
    universe  <- list(nodes = vn, levels = levels, nlev   = di)
    universe
}



## ##################################################################
##
## INTERNAL UTILITIES
##
## Used only in compileCPT
##
## ##################################################################

#' @rdname components_gather
#' @param xi cpt in some representation
#' @export
parse_cpt <- function(xi){
    UseMethod("parse_cpt")
}

#' @export
parse_cpt.xtabs <- function(xi){
    ## cat("parse_cpt.xtabs\n")
    NextMethod("parse_cpt")
}


#' @export
parse_cpt.default <- function(xi){
    ## cat("parse_cpt.default\n")
    if (!is.named.array(xi)) stop("'xi' must be a named array")
    .parse_cpt_finalize(varNames(xi), valueLabels(xi)[[1]],
                        as.numeric(xi), 0)
}


.parse_cpt_finalize <- function(vpar, vlev, values, smooth){

    ## str(list(vpar=vpar, vlev=vlev, values=values, smooth=smooth))
    ## Normalization of CPTs happen here
    ## str(list(vpar=vpar, vlev=vlev, values=values, smooth=smooth))
    values <- matrix(values, nrow=length(vlev))
    s  <- colSums(values)
    for (j in 1:ncol(values)) values[, j] <- values[, j] / s[j]
    values <- as.numeric(values)

    out <- list(vnam=vpar[1], vlev=vlev, vpar=vpar, values=values,
                normalize="first", smooth=smooth)
    class(out) <- "cpt_generic"
    out    
}


## NOTE : All this cptable stuff must be kept in the package in order
## to make bnlearn work (and moreover, cptable is used in the original
## JSS paper)

#' @export
parse_cpt.cptable_class <- function(xi) {
    ## cat("parse_cpt.cptable_class\n")
    .parse_cptable_finalize(attr(xi, "vpa"), attr(xi, "levels"), 
                        as.numeric(xi), attr(xi, "smooth"))
}

.parse_cptable_finalize <- function(vpar, vlev, values, smooth) {

    ## str(list(vpar=vpar, vlev=vlev, values=values, smooth=smooth))
    ## Normalization of CPTs happen here
    ## str(list(vpar=vpar, vlev=vlev, values=values, smooth=smooth))
    values <- matrix(values, nrow=length(vlev))
    s  <- colSums(values)
    for (j in 1:ncol(values)) values[, j] <- values[, j] / s[j]
    values <- as.numeric(values)

    out <- list(vnam=vpar[1], vlev=vlev, vpar=vpar, values=values,
                normalize="first", smooth=smooth)
    class(out) <- "cpt_generic"
    out    
}


## ##################################################################
##
## Extend compilation function 
##
## ##################################################################

## compile.cpt_representation <- function(object, ...)
##     compileCPT(object, ...)

## compile.pot_representation <- function(object, ...)
##     compilePOT(object, ...)

## #################################################################
