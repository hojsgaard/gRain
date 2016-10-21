#' @title Compile a graphical independence network (a Bayesian network)
#' 
#' @description Compiles a Bayesian network. This means creating a
#'     junction tree and establishing clique potentials.
#'
#' @name grain-compile
#' 
#' @aliases compile.grain compile.CPTgrain compile.POTgrain
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
#' @export compile.grain
compile.grain <-
  function(object, propagate=FALSE, root=NULL,
           control=object$control, details=0,...) {    #method <- match.arg(tolower(method), c("mcwh","r"))
    NextMethod("compile")
}

#' @rdname grain-compile
compile.CPTgrain <-
  function(object, propagate=FALSE, root=NULL, control=object$control, details=0, ...){

    mdagM <- moralizeMAT(as(object$dag,"Matrix"))
    vn    <- colnames(mdagM)
    nlev  <- object$universe$nlev[vn]

    if (!is.null(root) && length(root)>1){
        mdagM <- .setRoot( mdagM, root )
    }

    ugM   <- triangulateMAT(mdagM)
    .rip   <- ripMAT( ugM )
    ## FAST jt    <- .createJTreeGraph(.rip)
    ug    <- as(ugM,   "graphNEL")
    mdag  <- as(mdagM, "graphNEL")

### Insert potentials; ## Input: rip, universe, cptlist
    pot.with.1   <- .createPotList( .rip, object$universe )
    origpot    <- temppot <- .insertCPT(object$cptlist, pot.with.1, details)
    equipot   <- .insertNA(pot.with.1)

### Collect results
    ans  <- list(rip         = .rip,
                 ug          = ug,
                 equipot     = equipot,
                 temppot     = temppot,
                 origpot     = origpot,
                 details     = details )
    ans        <- c(object, ans)
    class(ans) <- class(object)

    ans$isCompiled   <- TRUE
    ans$isPropagated <- FALSE
    ans$control      <- control

    if (propagate){     ## Propagate if asked to
      .infoPrint(details, 1, cat (".Initializing network\n"))
      ans <- propagate(ans)
    }
    return(ans)
}

## NOTICE: the compiled object will contain a dag and a cptlist.
## These are not used for any calculations; only used for saving
## the network in Hugin format...

#' @rdname grain-compile
compile.POTgrain <-
  function(object, propagate=FALSE, root=NULL,
           control=object$control, details=0,...) {

      t00 <-  proc.time()
      ##FAST jt  <- .createJTreeGraph(object$rip)
      ans     <- list(## FAST jt          = jt,
                      temppot   = object$equipot,
                      origpot   = object$equipot,
                      ## FAST mdag        = object$ug,
                      details     = details )

      object$equipot   <- .insertNA(object$equipot)
                                        #object$details <- NULL
      ans            <- c(object, ans)
      class(ans)     <- class(object)

      ans$isCompiled   <- TRUE
      ans$isPropagated <- FALSE
      ans$control      <- control
      .timing(" Time: (total) compile:", control, t00)

      if (propagate){     ## Propagate if asked to
          .infoPrint(details, 1, cat(".Initializing network\n"))
          ans             <- propagate(ans)
      }
      ans$dag <- NULL
      return(ans)
  }



.setRoot <- function(mdagM, root){
    vn    <- colnames(mdagM)
    dn <- dimnames(mdagM)
    ft <- names2pairs(match(root, vn),sort=FALSE, result="matrix")
    ft <- rbind(ft,ft[,2:1,drop=FALSE])
    mdagM <- .sparse_setXtf1(mdagM, ft)
    dimnames(mdagM) <- dn
    mdagM
}


.createJTreeGraph <- function(rip){
  if (length(rip$cliques)>1){
    ft <-cbind(rip$parents, 1:length(rip$parents))
    ft <- ft[ft[,1]!=0,, drop=FALSE]
    V <- seq_along(rip$parents)
    if (nrow(ft)==0){
      jt <- new("graphNEL", nodes = as.character(V), edgemode = "undirected")
    } else {
      jt <- graph::ftM2graphNEL(ft, V=as.character(V), edgemode="undirected")
    }
  } else {
    jt <- new("graphNEL", nodes = "1", edgemode = "undirected")
  }
  return(jt)
}


.timing <- function(text, control, t0){
  if (!is.null(control$timing) && control$timing)
    cat(sprintf("%40s", text), proc.time()-t0,"\n")

}








