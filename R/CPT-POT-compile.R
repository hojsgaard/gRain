#' @title Compile conditional probability tables / cliques potentials.
#' 
#' @description Compile conditional probability tables / cliques
#'     potentials as a preprocessing step for creating a graphical
#'     independence network
#'
#' @name compile_potentials
#' 
#' @aliases compileCPT summary.cpt_spec compilePOT print.cpt_spec
#'     summary.cpt_spec
#' @param x To \code{compileCPT} x is a list of conditional
#'     probability tables; to \code{compilePOT}, x is a list of clique
#'     potentials.
#' @param object A list of potentials or of CPTs.
#' @param forceCheck Controls if consistency checks of the probability
#'     tables should be made.
#' @param ... Additional arguments; currently not used.
#' 
#' @return \code{compileCPT} returns a list of class 'cptspec'
#'     \code{compilePOT} returns a list of class 'potspec'
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{extractCPT}}, \code{\link{extractPOT}}
#' 
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{http://www.jstatsoft.org/v46/i10/}.
#'
#' @keywords utilities
#' 

#' @rdname compile_potentials
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



#' @rdname compile_potentials
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

## OLD NAMES - KEEP THESE 
#' @rdname compile_potentials
compileCPT <- compile_cpt

#' @rdname compile_potentials
compilePOT <- compile_pot


## ##################################################################

#' @rdname compile_potentials
parse_cpt <- function(xi){
    UseMethod("parse_cpt")
}

#' @rdname compile_potentials
parse_cpt.cptable_class <- function(xi){
    .parse_cpt_finalize(varNames(xi), valueLabels(xi)[[1]],
                        as.numeric(xi), attr(xi, "smooth"))
}

#' @rdname compile_potentials
#' @param xi Conditional probability table
parse_cpt.xtabs <- function(xi){
    NextMethod("parse_cpt")
}

#' @rdname compile_potentials
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

## ##################################################################

#' @rdname compile_potentials
compile.pot_rep <- function(object, ...)
    compile_pot(object)

#' @rdname compile_potentials
compile.cpt_rep <- function(object, ...)
    compile_cpt(object)


## ##################################################################




print.cpt_spec <- function(x,...){
    cat("cpt_spec with probabilities:\n")
    lapply(x,
           function(xx){
               vn <- varNames(xx)
               .print_probability(vn)
           })
  invisible(x)
}

print.pot_spec <- function(x,...){
    cat("pot_spec with potentials:\n")
    lapply(x,
           function(xx){
               vn <- names(dimnames(xx))
               cat("(", paste(vn, collapse=' '),") \n")
           })    
  invisible(x)
}

summary.cpt_spec <- function(object, ...){
    cat("cpt_spec with probabilities:\n")
    lapply(object,
           function(xx){
               vn <- varNames(xx)
               .print_probability(vn)               
           })
  invisible(object)
}

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


















