#' @title Compile conditional probability tables / cliques potentials.
#' 
#' @description Compile conditional probability tables / cliques
#'     potentials as a preprocessing step for creating a graphical
#'     independence network
#'
#' @name compile_components
#' 
#' @aliases compileCPT compilePOT
#' 
#' @param x To \code{compile_cpt} x is a list of conditional
#'     probability tables; to \code{compile_pot}, x is a list of clique
#'     potentials.
#'
## #' @param object A list of potentials or of CPTs.
#'
#' @param forceCheck Controls if consistency checks of the probability
#'     tables should be made.
#' @param ... Additional arguments; currently not used.
#' 
#' @details \code{compileCPT}, \code{compilePOT} are wrappers for
#'     \code{compile_cpt} and \code{compile_pot} and are kept for
#'     backward compatibility.
#'
#' @return A list with a class attribute.
#' 
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#'
#' @seealso \code{\link{extract_cpt}}, \code{\link{extract_pot}}
#' 
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{http://www.jstatsoft.org/v46/i10/}.
#'
#' @keywords utilities
#' 

#' @rdname compile_components
compile_cpt <- function(x, forceCheck=TRUE){
        
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

    if (!is.list(x)) stop("A list is expected")    

    ## zz: Internal representation of cpts
    zz  <- lapply(x, parse_cpt)

    uni <- .create_universe(zz)

    ## Given node names; need to check that they are not replicated
    vn_given <- sapply(zz, "[[", "vnam")
    if (length(vn_given) != length(unique(vn_given)))
        stop("Some nodes specified more than once: ", toString(vn_given))
    
    ## Are all cpts defined?
    ss <- setdiff(unique(unlist(vn_given)),  uni$nodes)
    if (length(ss) > 0)
        stop(paste("Distribution not specified for nodes(s):", toString(ss)))
    
    ## Does specification define a DAG?
    vp <- lapply(zz, "[[", "vpar")
    dg <- dagList(vp, forceCheck=forceCheck)

    ## Need list of cpts (each represented as an array)
    out <- lapply(seq_along(zz), .create_array, zz, uni)    
    names(out) <- uni$nodes

    ## Wrap it up
    attr(out, "universe") <- uni
    attr(out, "dag") <- dg
    class(out) <- "cpt_spec"
    out
}


#' @rdname compile_components
compile_pot <- function(x){
    
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

    if (!inherits(x, "pot_rep")) stop("can not compile 'x'\n")
    if (is.null(attr(x, "rip"))) stop("no rip attribute; not a proper POT_spec object")

    attr(x, "universe") <- .make.universe(x)
    attr(x, "ug")       <- ug(attr(x, "rip")$cliques)
    class(x) <- "pot_spec"
    x
}



## ##################################################################
## These parse_cpt functions are used only in compile_cpt
## ##################################################################

parse_cpt <- function(xi){
    UseMethod("parse_cpt")
}

parse_cpt.xtabs <- function(xi){
    NextMethod("parse_cpt")
}

parse_cpt.cptable_class <- function(xi){
    .parse_cpt_finalize(varNames(xi), valueLabels(xi)[[1]],
                        as.numeric(xi), attr(xi, "smooth"))
}

parse_cpt.default <- function(xi){
    if (!is.named.array(xi)) stop("'xi' must be a named array")
    .parse_cpt_finalize(varNames(xi), valueLabels(xi)[[1]],
                        as.numeric(xi), 0)
}

.parse_cpt_finalize <- function(vpar, vlev, values, smooth){
    out <- list(vnam=vpar[1], vlev=vlev, vpar=vpar, values=values,
                normalize="first", smooth=smooth)
    class(out) <- "cpt_generic"
    out    
}

## #############################################################

print.cpt_spec <- function(x, ...){
    cat("cpt_spec with probabilities:\n")
    lapply(x,
           function(xx){
               vn <- varNames(xx)
               .print_probability(vn)
           })
  invisible(x)
}

print.pot_spec <- function(x, ...){
    cat("pot_spec with potentials:\n")
    lapply(x,
           function(xx){
               vn <- names(dimnames(xx))
               cat("(", paste(vn, collapse=' '),") \n")
           })    
  invisible(x)
}

## FIXME: print.marg_spec missing
## FIXME: summary.pot_spec missing
## FIXME: summary.marg_spec missing


summary.cpt_spec <- function(object, ...){
    cat("cpt_spec with probabilities:\n")
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


## ##################################################################

compile.cpt_rep <- function(object, ...)
    compile_cpt(object)

compile.pot_rep <- function(object, ...)
    compile_pot(object)


## #################################################################

## For bacward compatibility (book, paper etc)
compileCPT <- compile_cpt
compilePOT <- compile_pot
















