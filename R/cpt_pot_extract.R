## ######################################################################
#' @title Extract conditional probabilities and clique potentials from
#'     data.
#' @description Extract list of conditional probability tables and
#'     list of clique potentials from data.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @name components_extract
## ######################################################################
#'
#' @details If \code{smooth} is non-zero then \code{smooth} is added
#'     to all cell counts before normalization takes place.
#' 
#' @param data_ A named array or a dataframe.
#'
#' @param graph An \code{igraph} object or a list or formula which can be
#'     turned into a \code{igraph} object by calling \code{ug} or
#'     \code{dag}. For \code{extract_cpt}, graph must be/define a DAG while for
#'     \code{extract_pot}, graph must be/define undirected triangulated graph.
#' 
#' @param smooth See 'details' below.
#' 
#' @return
#'   * \code{extract_cpt}: A list of conditional probability tables.
#'   * \code{extract_pot}: A list of clique potentials.
#'   * \code{extract_marg}: A list of clique marginals. 
#'
#' @seealso \code{\link{compileCPT}}, \code{\link{compilePOT}},
#'     \code{\link{grain}}
#'
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{https://www.jstatsoft.org/v46/i10/}.
#' @keywords utilities
#'
#' @examples
#' 
#' ## Extract cpts / clique potentials from data and graph
#' # specification and create network. There are different ways:
#'
#' data(lizard, package="gRbase")
#'
#' # DAG: height <- species -> diam
#' daG <- dag(~species + height:species + diam:species, result="igraph")
#'
#' # UG : [height:species][diam:species]
#' uG  <- ug(~height:species + diam:species, result="igraph")
#' 
#' pt <- extract_pot(lizard, ~height:species + diam:species) 
#' cp <- extract_cpt(lizard, ~species + height:species + diam:species)
#'
#' pt
#' cp
#'
#' # Both specify the same probability distribution
#' tabListMult(pt) |> as.data.frame.table()
#' tabListMult(cp) |> as.data.frame.table()
#'
#' \dontrun{
#' # Bayesian networks can be created as
#' bn.uG   <- grain(pt)
#' bn.daG  <- grain(cp)
#'
#' # The steps above are wrapped into a convenience method which
#' # builds a network from at graph and data.
#' bn.uG   <- grain(uG, data=lizard)
#' bn.daG  <- grain(daG, data=lizard)
#' }
#' 
#' @rdname components_extract
#' @export 
extract_cpt <- function(data_, graph, smooth=0) {

    is_valid_data(data_)
    if (inherits(graph, c("formula", "list")))
        graph <- dag(graph)
    if (!is_dag(graph))
        stop("'graph' not a DAG")
    
    vpa <- vpar(graph)
    out <- extract_cpt_worker(data_, vpa_=vpa, smooth=smooth)
    attr(out, "graph") <- graph
    class(out)         <- "cpt_representation"
    out
}

#' @export
#' @rdname components_extract
extract_pot <- function(data_, graph, smooth=0) {
    
    is_valid_data(data_)
    if (inherits(graph, c("formula", "list")))
        graph <- ug(graph)    
    if (!is_tug(graph))
        stop("'graph' not undirected and triangulated")

    rip_  <- rip(graph)
    out   <- extract_pot_worker(data_, rip_$cliques, rip_$sep, smooth=smooth)

    ## FIXME These attributes should be moved elsewhere
    attr(out, "rip")   <- rip_
    attr(out, "graph") <- graph    
    class(out)         <- "pot_representation"
    out
}

#' @export 
#' @rdname components_extract
extract_marg <- function(data_, graph, smooth=0) {

    is_valid_data(data_)
    
    if (inherits(graph, c("formula", "list")))
        graph <- ug(graph)    
    if (!is_tug(graph))
        stop("'graph' not undirected and triangulated")
    
    rip_  <- rip(graph)
    out   <- extract_marg_worker(data_, rip_$cliques, rip_$sep, smooth=smooth)
    attr(out, "rip")   <- rip_
    attr(out, "graph") <- graph    
    class(out)         <- "marg_representation"
    out
}


#' @export 
#' @rdname components_extract
#' @param marg_rep An object of class \code{marg_rep}
marg2pot <- function(marg_rep) {
    
    if (!inherits(marg_rep, "marg_representation"))
        stop("'marg_rep' not a marg_representation object\n")

    rip_ <- attr(marg_rep, "rip")
    seps <- rip_$separators
    pt <- lapply(seq_along(rip_$cliques),
                 function(i) {
                     if (length(seps[[i]]) == 0)
                         marg_rep[[i]]
                     else
                         tabDiv0(marg_rep[[i]], tabMarg(marg_rep[[i]], seps[[i]]))               
                 })
    attr(pt, "rip") <- rip_
    class(pt) <- "pot_representation"
    pt
}


#' @export 
#' @rdname components_extract 
#' @param pot_rep An object of class \code{pot_representation}
pot2marg <- function(pot_rep) {
    if (!inherits(pot_rep, "pot_representation"))
        stop("'pot_rep' not a pot_representation object\n")    
    rip_ <- attr(pot_rep, "rip")
    seps <- rip_$separators
    par  <- rip_$parents

    mg <- pot_rep    
    for (i in 2:length(rip_$cliques)) {
        if (par[i] > 0){
            mg[[i]] <- tabMult(mg[[i]], tabMarg(mg[[par[i]]], seps[[i]]))
        }
    }
    class(mg) <- "marg_rep"
    mg
}






## ##################################################################
##
## worker functions below here
##
## ##################################################################


extract_cpt_worker <- function(data_, vpa_, smooth=0) {
    
    is.df <- is.data.frame(data_)

    out <- lapply(vpa_, function(ss){
        marginal_data(data_, ss, is.df)
    })
    
    out <- vector("list", length(vpa_))
    for (i in seq_along(vpa_)){
        cq   <- vpa_[[ i ]]
        t.cq <- marginal_data(data_, cq, is.df) + smooth       
        out[[i]] <- tabNormalize(t.cq + smooth, type="first")
    }
    chk <- unlist(lapply(out,
                         function(zz) {
                             any(is.na(zz))
                         }))

    nas <- names(vpa_)[chk]

    if (length(nas) > 0) {
        cat(sprintf("NAs found in cpt(s) for node(s): %s\n", toString(nas)))
        cat(sprintf("  ... consider using the smooth argument\n"))
    }
    out
}


extract_pot_worker <- function(data_, cliq, seps=NULL, smooth=0) {        
    
    .normalize <- function(tt, sp) {
        if (length(sp) > 0){ 
          mm <- tabMarg(tt, sp)
          tabDiv(tt, mm)
        } else {
          mm <- sum(tt)
          tt / mm        
        }
    }

    is.df <- is.data.frame(data_)
    
    out <- vector("list", length(cliq))
    for (i in seq_along(cliq)){
        cq   <- cliq[[ i ]]
        sp   <- seps[[ i ]]
        t.cq <- marginal_data(data_, cq, is.df) + smooth       
        out[[i]] <- .normalize(t.cq, sp)
    }
    
    chk <- sapply(out, function(z) any(is.na(z)))
    nms <- lapply(out, function(z) names(dimnames(z)))
    
    if (any(chk)) {
        cat(sprintf("NAs found in clique potential(s) for clique(s): \n"))
        cat(toString(nms[chk]), "\n")
        cat(sprintf("  ... consider using the smooth argument\n"))
    }
    
    out
}


extract_marg_worker <- function(data_, marg, seps=NULL, smooth=0) {        
    out <- vector("list", length(marg))
    is.df <- is.data.frame(data_)
    
    for (i in seq_along(marg)){
        cq   <- marg[[ i ]]
        t.cq <- marginal_data(data_, cq, is.df) + smooth       
        out[[i]] <- t.cq / sum(t.cq)
    }
    out
}




## helper function; can possibly be made faster
marginal_data <- function(data_, set, is.df=NULL) {

    ## .dfMarg can possibly be made faster
    .dfMarg <- function(data_, set) {
        xtabs(~., data=data_[ , set, drop=FALSE])
    }

    if (is.null(is.df))
        is.df <- is.data.frame(data_)

    if (is.df)
        .dfMarg(data_, set)
    else
        tabMarg(data_, set)
        
}

is_valid_data <- function(data_) {
    if (!(is.data.frame(data_) || is.named.array(data_)))
        stop("'data_' must be dataframe or array.")
}


## ---------------------------------------------------------------
## FIXME for backward compatibility; deprecate in future release
## ---------------------------------------------------------------



#' @name old_components_extract
#' @inherit components_extract
#' @concept old_names
#' @export 
extractCPT <- extract_cpt

#' @rdname old_components_extract
#' @export 
extractPOT <- extract_pot

#' @rdname old_components_extract
#' @export 
extractMARG <- extract_marg


