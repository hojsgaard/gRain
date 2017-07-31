#' @title Extract conditional probabilities and clique potentials from
#'     data.
#' 
#' @description Extract list of conditional probability tables and
#'     list of clique potentials from data.
#'
#' @name extract-cpt
#' 
#' @details If \code{smooth} is non--zero then \code{smooth} is added
#'     to all cell counts before normalization takes place.
#' 
#' @aliases extractCPT extractCPT.table extractCPT.data.frame
#'     extractPOT extractPOT.table extractPOT.data.frame
#' @param data_ A named array or a dataframe.
#' @param graph A graphNEL object or a list or formula which can be
#'     turned into a graphNEL object by calling \code{ug} or
#'     \code{dag}. For extractCPT, graph must be/define a DAG while for
#'     extractPOT, graph must be/define undirected triangulated graph.
#' @param smooth See 'details' below.
#' @return extractCPT: A list of conditional probability tables
#'     extractPOT: A list of clique potentials.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{compileCPT}}, \code{\link{compilePOT}},
#'     \code{\link{grain}}
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords utilities
#' @examples
#' 
#' 
#' ## Asia (chest clinic) example:
#' 
#' ## Version 1) Specify conditional probability tables.
#' yn <- c("yes","no")
#' a    <- cptable(~asia, values=c(1,99),levels=yn)
#' t.a  <- cptable(~tub+asia, values=c(5,95,1,99),levels=yn)
#' s    <- cptable(~smoke, values=c(5,5), levels=yn)
#' l.s  <- cptable(~lung+smoke, values=c(1,9,1,99), levels=yn)
#' b.s  <- cptable(~bronc+smoke, values=c(6,4,3,7), levels=yn)
#' e.lt <- cptable(~either+lung+tub,values=c(1,0,1,0,1,0,0,1),levels=yn)
#' x.e  <- cptable(~xray+either, values=c(98,2,5,95), levels=yn)
#' d.be <- cptable(~dysp+bronc+either, values=c(9,1,7,3,8,2,1,9), levels=yn)
#' plist <- compileCPT(list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be))
#' pn1 <- grain(plist)
#' q1 <- querygrain(pn1)
#' 
#' ## Version 2) Specify DAG and data
#' data(chestSim100000, package="gRbase")
#' dgf   <- ~asia + tub * asia + smoke + lung * smoke +
#'          bronc * smoke + either * tub * lung +
#'          xray * either + dysp * bronc * either
#' dg    <- dag(dgf)
#' pp    <- extractCPT(chestSim100000, dg)
#'
#' pn2   <- grain(pp)
#' ## Same as:
#' cpp2  <- compileCPT(pp)
#' pn2   <- grain(cpp2)
#' 
#' q2    <- querygrain(pn2)
#' 
#' ## Version 2) Specify triangulated undirected graph and data
#' ugf <- list(c("either", "lung", "tub"), c("either", "lung", "bronc"), 
#'     c("either", "xray"), c("either", "dysp", "bronc"), c("smoke", 
#'     "lung", "bronc"), c("asia", "tub"))
#' gg    <- ug(ugf)
#' pp    <- extractPOT(chestSim100000, gg)
#'
#' pn3   <- grain(pp)
#' ## Same as:
#' cpp3  <- compilePOT(pp)
#' pn3   <- grain(cpp3)
#'
#' q3    <- querygrain(pn3)
#' 
#' ## Compare results:
#' str(q1)
#' str(q2[names(q1)])
#' str(q3[names(q1)])
#' 
#' @rdname extract-cpt

extractCPT <- function(data_, graph, smooth=0){
    .is.valid.data(data_)

    if (inherits(graph, c("formula", "list")))
        graph <- dag(graph)

    if (!is.DAG(graph)) stop("'graph' not a DAG")

    vpa <- vpar(graph)
    out <- .extractCPT_(data_, vpa=vpa, smooth=smooth)
    ##FIXME: Should any info be stored in the output? vpa for example?
    class(out) <- "CPT_rep"
    out
}




.is.valid.data <- function(data_){
    if ( !(is.data.frame(data_) || is.named.array(data_)) )
        stop("'data_' must be dataframe or array.")
}

.extractCPT_ <- function(data_, vpa, smooth=0){

    if (is.data.frame(data_)){
        out <- lapply(vpa,
                      function(ss){ xtabs(~., data=data_[, ss, drop=FALSE]) })
    } else {
        out <- lapply(vpa,
                      function(ss){ tabMarg(data_, ss) })
    }
    
    ## FIXME : Get rid of this parray stuff (at least as a class)
    ## FIXME: Normalization takes place here
    out <- lapply(out, as.parray, normalize="first", smooth=smooth)
    
    chk <- unlist(lapply(out, function(zz) any(is.na(zz))))
    nnn <- names(chk)[which(chk)]
    if (length(nnn) > 0){
        cat(sprintf("NAs found in conditional probability table(s) for nodes: %s\n", toString(nnn)))
        cat(sprintf("  ... consider using the smooth argument\n"))
    }

    out
}

#' @rdname extract-cpt
extractPOT <- function(data_, graph, smooth=0){
    .is.valid.data(data_)

    if (inherits(graph, c("formula", "list")))
        graph <- ug(graph)
    
    if (!is.TUG(graph)) stop("'graph' not undirected and triangulated")
    rp  <- rip( graph )
    
    out <- .extractPOT_(data_, rip=rp, smooth=smooth)
    attr(out, "rip")     <- rp
    
    class(out) <- "POT_rep"
    out
}

.extractPOT_ <- function(data_, rip, smooth=0){

    .normalize <- function(tt, sp){
        if (length(sp) > 0) tabDiv0(tt, tabMarg(tt, sp))
        else tt / sum(tt)        
    }
    
    .extractPOT_table <- function(data_, cliq, seps=NULL, smooth=0){
        out <- vector("list", length(cliq))
        for ( i  in seq_along(cliq)){
            cq    <- cliq[[ i ]]
            sp    <- seps[[ i ]]
            ##str(list(cq=cq, sp=sp))
            t.cq  <- tabMarg(data_, cq) + smooth
            out[[i]] <- .normalize(t.cq, sp)
        }
        out
    }
    
    .extractPOT_dataframe <- function(data_, cliq, seps=NULL, smooth=0){        
        out <- vector("list", length(cliq))
        for ( i  in seq_along(cliq)){
            cq   <- cliq[[ i ]]
            sp   <- seps[[ i ]]
            ##str(list(cq=cq, sp=sp))
            
            t.cq  <- xtabs(~., data=data_[ , cq, drop=FALSE]) + smooth                       
            out[[i]] <- .normalize(t.cq, sp)
        }
        out
    }
    

    if (is.data.frame(data_)){
        .extractPOT_dataframe(data_, rip$cliques, rip$sep, smooth=smooth)
    } else {
        .extractPOT_table(data_, rip$cliques, rip$sep, smooth=smooth)
    }
}

#' @rdname extract-cpt
extractMARG <- function(data_, graph, smooth=0){
    .is.valid.data(data_)
    if (!is.TUG(graph))
        stop("'graph' not undirected and triangulated")
    
    rp  <- rip( graph )

    out <- .extractMARG_(data_, rip=rp, smooth=smooth)
    attr(out, "rip")     <- rp
        
    class(out) <- "MARG_rep"
    out
}


.extractMARG_ <- function(data_, rip, smooth=0){

    .extractMARG_table <- function(data_, cliq, seps=NULL, smooth=0){
        out <- vector("list", length(cliq))
        for ( i  in seq_along(cliq)){
            cq    <- cliq[[ i ]]
            t.cq  <- tabMarg(data_, cq) + smooth
            out[[i]] <- t.cq / sum(t.cq)
        }
        out
    }

    
    .extractMARG_dataframe <- function(data_, cliq, seps=NULL, smooth=0){        
        out <- vector("list", length(cliq))
        for ( i  in seq_along(cliq)){
            cq   <- cliq[[ i ]]
            t.cq  <- xtabs(~., data=data_[ , cq, drop=FALSE]) + smooth
            out[[i]] <- t.cq / sum(t.cq)
        }
        out
    }
            
    if (is.data.frame(data_)){
        .extractMARG_dataframe(data_, rip$cliques, rip$sep, smooth=smooth)
    } else {
        .extractMARG_table(data_, rip$cliques, rip$sep, smooth=smooth)
    }
}


## FIXME: or cpt_from_data
#' @rdname extract-cpt
data2cpt <- extractCPT

#' @rdname extract-cpt
data2pot <- extractPOT

#' @rdname extract-cpt
data2marg <- extractMARG


#' @rdname extract-cpt
#' @param mg An object of class \code{MARG_rep}
marg2pot <- function(mg){
    if (!inherits(mg, "MARG_rep")) stop("'mg' not a MARG_rep object\n")
    rp <- attr(mg, "rip")
    seps <- rp$separators
    pt <- lapply(seq_along(rp$cliques),
                 function(i){
                     if (length(seps[[i]]) == 0)
                         mg[[i]]
                     else
                         tabDiv0(mg[[i]], tabMarg(mg[[i]], seps[[i]]))               
                 })
    attr(pt, "rip") <- rp
    class(pt) <- "POT_rep"
    pt
}

#' @rdname extract-cpt
#' @param pt An object of class \code{POT_rep}
pot2marg <- function(pt){
    if (!inherits(pt, "POT_rep")) stop("'pt' not a POT_rep object\n")    
    mg <- pt
    rp <- attr(pt, "rip")
    seps <- rp$separators
    par  <- rp$parents
    
    for (i in 2:length(rp$cliques)){
        if (par[i] > 0){
            mg[[i]] <- tabMult(mg[[i]], tabMarg(mg[[par[i]]], seps[[i]]))
        }
    }
    class(mg) <- "MARG_rep"
    mg
}




