## FIXME compileCPT: Det skal rettes...
## compileCPT( list(p1, p2) )
## p1 <- cptable(~a|b, levels=uni2, values=c(7, 3, 9, 1))
## p2 <- cptable(~b, levels=uni2, values=c(7, 3))
## compileCPT( list(p1, p1, p2) )

#' @title Compile conditional probability tables / cliques potentials.
#' 
#' @description Compile conditional probability tables / cliques
#'     potentials as a preprocessing step for creating a graphical
#'     independence network
#'
#' @name compile-cpt
#' 
#' @aliases compileCPT summary.CPTspec compilePOT print.CPTspec
#'     summary.CPTspec
#' @param x To \code{compileCPT} x is a list of conditional
#'     probability tables; to \code{compilePOT}, x is a list of clique
#'     potentials.
#' @param object A list of potentials or of CPTs.
#' @param forceCheck Controls if consistency checks of the probability
#'     tables should be made.
#' @param details Controls amount of print out. Mainly for debugging
#'     purposes.
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

#' @rdname compile-cpt
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
    zz <- lapply(x, parse_cpt)
    uni <- .create_universe(zz)

    ## Does specification define a DAG?
    ## FIXME what if forceCheck = FALSE???
    vp <- lapply(zz, "[[", "vpar")
    dg <- dagList(vp, forceCheck=forceCheck)

    if (forceCheck){
        ## Is distribution specified for all variables?
        ss <- setdiff(unique.default(unlist(vp)),  uni$nodes)
        if (length(ss) > 0)
            stop(paste("distribution not specified for variable(s):", toString(ss)))
    }

    ## Need list of cpts (each represented as an array)
    out <- lapply(seq_along(zz), .create_array, zz, uni)    
    names(out) <- uni$nodes
    
    attr(out, "universe") <- uni
    attr(out, "dag") <- dg
    class(out) <- "CPTspec"
    out
}


parse_cpt <- function(xi){
    UseMethod("parse_cpt")
}

parse_cpt.cptable <- function(xi){
    .parse_cpt_finalize(varNames(xi), valueLabels(xi)[[1]], as.numeric(xi), attr(xi, "smooth"))
}

parse_cpt.xtabs <- function(xi){
    NextMethod("parse_cpt")
}

parse_cpt.default <- function(xi){
    if (!is.named.array(xi)) stop("'xi' must be a named array")
    .parse_cpt_finalize(varNames(xi), valueLabels(xi)[[1]], as.numeric(xi), 0)
}

.parse_cpt_finalize <- function(vpar, vlev, values, smooth){
    out <- list(vnam=vpar[1], vlev=vlev, vpar=vpar, values=values,
                normalize="first", smooth=smooth)
    class(out) <- "cpt_generic"
    out    
}


as_cptlist <- function(x, forceCheck=FALSE){
    if (!(inherits(x, "list") &&
          all(unlist(lapply(x, function(a) is.named.array(a))))))
        stop("input must be named list of named arrays")
    
    if (forceCheck){
        compile_cpt(x, forceCheck=TRUE)
    } else {
        out <- x
        vl  <- lapply(out, function(g) dimnames(g)[[1]])
        vn  <- names(vl)
        di  <- unlist(lapply(vl, length))
        universe   <- list(nodes = vn, levels = vl, nlev = di)
        attr(out, "universe") <- universe
        ## FIXME do we need dag here????
        class(out) <- "CPTspec_simple"
        out
    }    
}





#' @rdname compile-cpt
compilePOT <- function(x){
    
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

    if (!inherits(x, "POT_rep")) stop("can not compile 'x'\n")
    if (is.null(attr(x, "rip"))) stop("no rip attribute; not a proper POT_spec object")

    attr(x, "universe") <- .make.universe(x)
    attr(x, "ug")       <- ug(attr(x, "rip")$cliques)
    class(x) <- "POTspec"
    x
}

#' @rdname compile-cpt
compile.POT_rep <- function(object, ...)
    compilePOT(object)

#' @rdname compile-cpt
compile.CPT_rep <- function(object, ...)
    compileCPT(object)

print.CPTspec <- function(x,...){
    cat("CPTspec with probabilities:\n")
    lapply(x,
           function(xx){
               vn <- varNames(xx)
               if (length(vn) > 1){
                   cat(paste(" P(", vn[1], "|", paste(vn[-1], collapse=' '), ")\n"))
               } else {
                   cat(paste(" P(", vn, ")\n"))
               }
           })
  invisible(x)
}

print.CPTspec <- function(x, ...){
    print.default(c(x))
}


summary.CPTspec <- function(object, ...){
    cat("CPTspec with probabilities:\n")
    lapply(object,
           function(xx){
               vn <- varNames(xx)
               if (length(vn) > 1){
                   cat(paste(" P(", vn[1], "|", paste(vn[-1], collapse=' '), ")\n"))
               } else {
                   cat(paste(" P(", vn, ")\n"))
               }
           })
  invisible(object)
}



as_CPTspec_simple <- function(x){
    z <- c(x)
    attr(z, "universe") <- attr(x, "universe")
    class(z) <- "CPTspec_simple"
    z
}

print.CPTspec_simple <- function(x,...){
    cat("CPTspec_simple with probabilities:\n")
    lapply(x,
           function(xx){
               vn <- varNames(xx)
               if (length(vn) > 1){
                   cat(paste(" P(", vn[1], "|", paste(vn[-1], collapse=' '), ")\n"))
               } else {
                   cat(paste(" P(", vn, ")\n"))
               }
           })
  invisible(x)
}

print.POTspec <- function(x,...){
    cat("POTspec with potentials:\n")
    lapply(x,
           function(xx){
               vn <- names(dimnames(xx))
               cat(   "(", paste(vn, collapse=' '),") \n")
           })
    
  invisible(x)
}

#' @rdname compile-cpt
compileCPT <- compile_cpt

.compileCPT <- function(x, forceCheck=TRUE, details=0){

    zz <- lapply(x, .parseCPT_item)

    vnamList <- lapply(zz, "[[", "vnam") ## variable name
    vparList <- lapply(zz, "[[", "vpar") ## variable and parent names
    vlevList <- lapply(zz, "[[", "vlev") ## variable level
    vn       <- unlist(vnamList, use.names=FALSE)
    names( vlevList ) <- vn
    
    ## Check for acyclicity
    if (details>=1) cat(". creating dag and checking for acyclicity...\n")
    dg  <- dagList(vparList)
    oo  <- topoSort(dg)
    if (length(oo) == 0)
        stop("Graph defined by the cpt's is not acyclical...\n");

    if (details>=1) cat(". creating probability tables ...\n")
    
    out        <- vector("list", length(vn))
    names(out) <- vn
    
    for ( i  in seq_along( vnamList ) ){
        if (!forceCheck && class(zz[[ i ]])[1]=="parray"){
            out[[ i ]] <- zz[[ i ]]
        } else {
            vpar  <- zz[[  i  ]]$vpar
            lev   <- vlevList[ vpar ]
            val   <- zz[[ i ]]$values

            if (forceCheck){
                mm    <- match(vpar, vn)
                if (any(is.na(mm))){
                    sss <- sprintf("compileCPT: Distribution not specified for node(s)\n %s \n",
                                   toString(vpar[which(is.na(mm))]))
                    stop(sss, call.=FALSE)
                }
                if ( prod( unlist(lapply(lev, length))) != length(val) ){
                    cat(sprintf("Error for v,pa(v): %s\n", toString(vpar)))
                    str( lev )
                    str( val )
                    stop("Table dimensions do not match!")
                }
            }

            ## FIXME: DIrty hack!
            uu <- parray(varNames  = vpar,
                         levels    = lev,
                         values    = val,
                         normalize = zz[[ i ]]$normalize,
                         smooth    = zz[[ i ]]$smooth
                         )
            class(uu) <- "array"
            out[[ i ]] <- uu                
        }
    }

    di         <- unlist(lapply(vlevList, length))

    universe        <- list(nodes = vn, levels = vlevList, nlev = di)
    attr(out, "universe") <- universe
    attr(out, "dag") <- dg ## FIXME: Not really needed to store dag
    class(out) <- "CPTspec"
    out
}




.parseCPT_item <- function(xi){ ## Create intermediate form of CPTs
    
    cls <- c("cptable", "parray", "array", "matrix", "xtabs", "table") ## FIXME: Fragile
    j   <- inherits(xi, cls, which=T)

    if (!any(u <- j > 0)) stop("xi not a valid object")
    
    cls <- cls[u]    
    vpar <- varNames(xi)
    vlev <- valueLabels(xi)[[1]]
    smooth <- if (cls[1] == "cptable") attr(xi, "smooth") else 0
    
    out <- list(vnam=vpar[1], vlev=vlev, vpar=vpar, values=as.numeric(xi),
                normalize="first", smooth=smooth)
    out
}





#' @rdname compile-cpt
.compile_cpt <- function(x, forceCheck=TRUE){

    if (!is.list(x)) stop("A list is expected")    
    ## FIXME: Fragile that there is no control over input
    ##zz <- lapply(x, .parseCPT_item)
    zz <- lapply(x, parse_cpt)

    vn <- unlist(lapply(zz, "[[", "vnam"))
    vl <- lapply(zz, "[[", "vlev")
    vd <- unlist(lapply(vl, length))
    vp <- lapply(zz, "[[", "vpar")
    di <- unlist(lapply(vl, length))
    
    names(vl) <- vn
    names(vd) <- vn

    universe  <- list(nodes = vn, levels = vl, nlev = di)

    ## FIXME: Perhaps not create dag here at all unless check is required
    dg <- dagList(vp, forceCheck = forceCheck)
    
    if (forceCheck){        
        ## Distribution specified for all variables?
        ss <- setdiff(unique.default(unlist(vp)),  vn)
        if (length(ss) > 0)
            stop(paste("distribution not specified for variable(s):", toString(ss)))

        ## Are dimensions consistent?
        lapply(1:length(vn), function(i){
            ##cat(c(i, prod(vd[vp[[i]]]), length(zz[[i]]$values)), "\n")
            if (length(zz[[i]]$values > 0)) { ## =0 can happen if cptable(..., values=NULL)
                if (prod(vd[vp[[i]]]) != length(zz[[i]]$values)){
                    cat(paste("problem here:\n",
                              "cpt:", toString(vp[[i]]), "\n",
                              "number of values given: ", length(zz[[i]]$values),"\n",
                              "number of values expected: ", prod(vd[vp[[i]]]), "\n"))
                    stop("can not proceed!")
                }
            }
            })
    }
    
    out <-
        lapply(1:length(vn), function(i){
            dn <- vl[vp[[i]]]
            di <- c(vd[vp[[i]]], use.names=FALSE)
            if (length(zz[[i]]$values) > 0)
                ar <- array(zz[[i]]$values + zz[[i]]$smooth, dim=di, dimnames=dn)
            else
                ar <- array(rep(1.0, prod(di)) + zz[[i]]$smooth, dim=di, dimnames=dn)
            ##ar <- tabNormalize(ar, type="first") ## FIXME Should happen at compile time???
            ar
        })
    names(out) <- vn

    attr(out, "universe") <- universe
    attr(out, "dag") <- dg 
    class(out) <- "CPTspec"
    out
}







    ## print(object)
  ## cat(sprintf("attributes:\n %s\n", toString(names(attributes(object)))))
  ## cat(sprintf("names:\n %s \n", toString(attributes(object)$names)))
  ## cat(sprintf("nodes:\n %s \n", toString(attributes(object)$nodes)))
  ## cat(sprintf("levels: \n"))
  ## str(attributes(object)$levels)
  ## cat(sprintf("vparList: \n"))
  ## str(attributes(object)$vparList)
  ## cat(sprintf("nlev: \n"))
  ## str(attributes(object)$nlev)
  ## cat(sprintf("dag: \n"))
  ## print(attributes(object)$dag)
  ## cat(sprintf("dagM: \n"))
  ## print(attributes(object)$dagM)
  ## return(invisible(object))
