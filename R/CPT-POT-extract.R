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
#' @param x An array or a dataframe.
#' @param graph A graph represented as a graphNEL object.  For
#'     extractCPT, graph must be a DAG while for extractPOT, graph
#'     must be undirected triangulated graph.
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
#' cpp2  <- compileCPT(pp)
#' pn2   <- grain(cpp2)
#' q2    <- querygrain(pn2)
#' 
#' ## Version 2) Specify triangulated undirected graph and data
#' ugf <- list(c("either", "lung", "tub"), c("either", "lung", "bronc"), 
#'     c("either", "xray"), c("either", "dysp", "bronc"), c("smoke", 
#'     "lung", "bronc"), c("asia", "tub"))
#' gg    <- ugList(ugf)
#' pp    <- extractPOT(chestSim100000, gg)
#' cpp3  <- compilePOT(pp)
#' pn3   <- grain(cpp3)
#' q3    <- querygrain(pn3)
#' 
#' ## Compare results:
#' str(q1)
#' str(q2[names(q1)])
#' str(q3[names(q1)])
#' 
#' @rdname extract-cpt
extractCPT <- function(x, graph, smooth=0){
    if ( !(is.data.frame(x) || is.named.array(x)) )
        stop("'x' must be dataframe or array.")
    ## if (!inherits(x, c("data.frame","table", "xtabs")))
    ##   stop("'x' must be one of dataframe, table, xtabs")
    if (!inherits(graph, "graphNEL"))
        stop("'graph' must be a graphNEL object")
    if (!is.DAG(graph))
        stop("'graph' must be a DAG")
    
    V   <- graph::nodes(graph)
    vpa <- vpar(graph)[V]
    
  if (class(x)[1]=="data.frame"){
      ### cat("extractCPT - data.frame\n")
      ans <- lapply(vpa, function(ss){
          ##cat(sprintf("---- %s ----\n", toString( ss )))
          zzz <- xtabs(~., data=x[, ss, drop=FALSE])
          ##print(zzz)
          zzz
      })
  } else {
    ans <- lapply(vpa, function(ss){tableMargin(x, ss)})
  }

  #.ans <<- ans
  ans <- lapply(ans, as.parray, normalize="first", smooth=smooth)
  chk <- unlist(lapply(ans, function(zz) any(is.na(zz))))
  nnn <- names(chk)[which(chk)]
  if (length(nnn)>0){
      cat(sprintf("NAs found in conditional probability table(s) for nodes: %s\n", toString(nnn)))
      cat(sprintf("  ... consider using the smooth argument\n"))
  }

  class(ans) <- c("extractCPT","list")
  ans
}

#' @rdname extract-cpt
extractPOT <- function(x, graph, smooth=0){
    if ( !(is.data.frame(x) || is.named.array(x)) )
        stop("'x' must be dataframe or array.")
    ## if (!inherits(x, c("data.frame","table")))
    ##   stop("'x' must be a dataframe or a table")
  if (!inherits(graph, "graphNEL"))
    stop("'graph' must be a graphNEL object")
  if (!is.TUG(graph))
    stop("'graph' must be a triangulated undirected graph")

  .rip  <- rip( graph )

  if (class(x)[1]=="data.frame"){
    ans <- .extractPOT_dataframe(x, .rip$cliques, .rip$sep, smooth=smooth)
  } else {
    ans <- .extractPOT_table(x, .rip$cliques, .rip$sep, smooth=smooth)
  }

  attr(ans, "rip")     <- .rip

  dg 	    <- ug2dag(graph)
  cptlist <- extractCPT(x, dg, smooth=smooth)
  attr(ans, "dag")     <- dg        ## Needed to save network in Hugin format
  attr(ans, "cptlist") <- cptlist   ## Needed to save network in Hugin format

  class(ans) <- c("extractPOT","list")
  ans
}

.extractPOT_table <- function(x, cliq, seps=NULL, smooth=0){
  ans <- vector("list", length(cliq))
  for ( i  in seq_along(cliq)){
    cq    <- cliq[[ i ]]
    sp    <- seps[[ i ]]
    t.cq  <- tableMargin(x, cq) + smooth
    names(dimnames(t.cq)) <- cq
    if (!is.null(seps) && length(sp)>0){
      t.sp      <- tableMargin(t.cq, sp)
      ans[[ i ]] <- tableOp2(t.cq, t.sp, op=`/`)
    } else {
      ans[[ i ]] <- t.cq / sum(t.cq)
    }
  }
  ans
}

.extractPOT_dataframe <- function(x, cliq, seps=NULL, smooth=0){
  ans <- vector("list", length(cliq))
  for ( i  in seq_along(cliq)){
    cq   <- cliq[[ i ]]
    sp   <- seps[[ i ]]
    ## FIXME: Isn't that the same 
    xxx  <- xtabs(~., data=x[ , cq, drop=FALSE])  ## cross classfy data in dataframe
    t.cq <- tableMargin(xxx, cq) + smooth         ## then marginalize
    ## FIXME: As
    t.cq  <- xtabs(~., data=x[ , cq, drop=FALSE]) + smooth

    
    names(dimnames(t.cq)) <- cq
    if (!is.null(seps) && length(sp)>0){
      t.sp       <- tableMargin(t.cq, sp)
      ans[[ i ]] <- tableOp2(t.cq, t.sp, op=`/`)
    } else {
      ans[[ i ]] <- t.cq / sum(t.cq)
    }
  }
  ans
}
