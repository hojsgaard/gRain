#' @title Propagate a graphical independence network (a Bayesian network)
#' 
#' @description Propagation refers to calibrating the cliques of the
#'     junction tree so that the clique potentials are consistent on
#'     their intersections
#'
#' @name grain-propagate
#' 
#' @details The \code{propagate} method invokes \code{propagateLS}
#'     which is a pure R implementation of the Lauritzen-Spiegelhalter
#'     algorithm.
#' 
#' The function \code{propagate__} invokes \code{propagateLS__} which
#' is a c++ implementation of the Lauritzen-Spiegelhalter algorithm.
#' 
#' The c++ based version is several times faster than the purely R
#' based version, and after some additional testing the c++ based
#' version will become the default.
#' 
#' @aliases propagate.grain propagateLS propagate__ propagateLS__
#' @param object A grain object
#' @param details For debugging info
#' @param ... Currently not used
#' @return A compiled and propagated grain object.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{grain}}, \code{\link[gRbase]{compile}}
#' @references Søren Højsgaard (2012). Graphical Independence
#'     Networks with the gRain Package for R. Journal of Statistical
#'     Software, 46(10), 1-26.
#'     \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords utilities models
#' @examples
#' 
#' 
#' yn   <- c("yes","no")
#' a    <- cptable(~asia,              values=c(1,99), levels=yn)
#' t.a  <- cptable(~tub+asia,          values=c(5,95,1,99), levels=yn)
#' s    <- cptable(~smoke,             values=c(5,5), levels=yn)
#' l.s  <- cptable(~lung+smoke,        values=c(1,9,1,99), levels=yn)
#' b.s  <- cptable(~bronc+smoke,       values=c(6,4,3,7), levels=yn)
#' e.lt <- cptable(~either+lung+tub,   values=c(1,0,1,0,1,0,0,1), levels=yn)
#' x.e  <- cptable(~xray+either,       values=c(98,2,5,95), levels=yn)
#' d.be <- cptable(~dysp+bronc+either, values=c(9,1,7,3,8,2,1,9), levels=yn)
#' plist <- compileCPT(list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be))
#' pn    <- grain(plist)
#' pnc  <- compile(pn, propagate=FALSE)
#' 
#' if (require(microbenchmark))
#'     microbenchmark(
#'         propagate(pnc),
#'         propagate__(pnc) )
#' 
#' @export propagate.grain
propagate.grain <- function(object, details=object$details, ...){

  t0 <- proc.time()
  ## propagate.grain: equipot is updated after propagation on temppot
  ## such that equipot will contain the updated potentials.
  ## object$equipot <- propagateLS(object$temppot,
  ##                               rip=object$rip, initialize=TRUE, details=details)

  object$equipot <-
      propagateLS(object$temppot, rip=object$rip)

  ## object$isInitialized <- TRUE
  object$isPropagated  <- TRUE

  ## FIXME: propagate.grain : Looks strange
  if ( !is.null(getEvidence(object)) ){
      ev <- getEvidence(object)
      attr(ev, "pEvidence") <- pEvidence(object)
      object$evidence <- ev
  }

  .timing(" Time: propagation:", object$control, t0)
  return(object)
}


#' @rdname grain-propagate
propagate__ <- function(object, details=object$details, ...){

  t0 <- proc.time()
  ## propagate.grain: equipot is updated after propagation on temppot
  ## such that equipot will contain the updated potentials.
  ## object$equipot <- propagateLS(object$temppot,
  ##                               rip=object$rip, initialize=TRUE, details=details)

  object$equipot <-
      propagateLS__(object$temppot, rip=object$rip)

  ## object$isInitialized <- TRUE
  object$isPropagated  <- TRUE

  ## FIXME: propagate.grain : Looks strange
  if ( !is.null(getEvidence(object)) ){
      ev <- getEvidence(object)
      attr(ev, "pEvidence") <- pEvidence(object)
      object$evidence <- ev
  }

  .timing(" Time: propagation:", object$control, t0)
  return(object)
}

## Lauritzen Spiegelhalter propagation
##

#' @rdname grain-propagate
#' @param cqpotList Clique potential list
#' @param rip A rip ordering
#' @param initialize Always true
propagateLS <- function(cqpotList, rip, initialize=TRUE, details=0){
    #details=20
    #cat(".Propagating BN: [propagateLS]\n")
    .infoPrint(details, 1, cat(".Propagating BN: [propagateLS]\n"))
    ## FIXME: Don't remember the idea behind the 'initialize' argument; should always be true

    cliq       <- rip$cliques
    seps       <- rip$separators
    pa         <- rip$parent
    childList  <- rip$childList

    ncliq      <- length(cliq)

    ## This assignment is needed because RIP now returns 0 as the
    ## parent index for the first clique
    pa[pa==0]<-NA

    ## Backward propagation (collect evidence) towards root of junction tree
    ##
    .infoPrint(details,2, cat("..BACKWARD:\n"))
    t0 <- proc.time()
    if (ncliq>1){
        for ( i  in ncliq:2){
            cq   <- cliq[[ i ]]
            sp   <- seps[[ i ]]
            .infoPrint2(details, 2, "Clique %d: {%s}\n",  i , .colstr( cq ))
            cq.pot   <- cqpotList[[ i ]]
            pa.pot   <- cqpotList[[pa[ i ]]]

            if (length(sp)>=1 && !is.na(sp)){
                .infoPrint2(details, 2, "Marg onto sep {%s}\n", .colstr(sp))
                sp.pot               <- tableMargin(cq.pot, sp)
                cqpotList[[ i ]]     <- tableOp2(cq.pot, sp.pot, `/`)
                cqpotList[[pa[ i ]]] <- tableOp2(pa.pot, sp.pot, `*`)
            } else{
                zzz               <- sum(cq.pot)
                cqpotList[[1]]    <- cqpotList[[1]] * zzz
                cqpotList[[ i ]]  <- cq.pot / zzz
            }
        }
    }

    normConst      <- sum(cqpotList[[1]])
    cqpotList[[1]] <- cqpotList[[1]] / normConst

    ## Forward propagation (distribute evidence) away from root of junction tree
    ##
    .infoPrint(details,2,cat("..FORWARD:\n"))
    t0 <- proc.time()
    for ( i  in 1:ncliq){
        .infoPrint2(details, 2, "Clique %d: {%s}\n",  i , .colstr(cliq[[ i ]]))
        ch <- childList[[ i ]]
        if (length(ch)>0)
            {
                .infoPrint2(details,2, "..Children: %s\n", .colstr(ch))
                for ( j  in 1:length(ch))
                    {
                        if (length(seps[[ch[ j ]]])>0)
                            {
                                .infoPrint2(details, 2, "Marg onto sep %i: {%s}\n", ch[ j ], .colstr(seps[[ch[ j ]]]))
                                sp.pot            <- tableMargin(cqpotList[[ i ]], seps[[ch[ j ]]])
                                ##cat(sprintf("......is.na: sp.pot=%i\n", any(is.na(sp.pot))))
                                cqpotList[[ch[ j ]]]  <- tableOp2(cqpotList[[ch[ j ]]], sp.pot, `*`)
                                .infoPrint(details, 4, { cat("Marginal:\n"); print (sp.pot) })
                            }
                    }
            }
    }

    attr(cqpotList, "pEvidence") <- normConst
    cqpotList
}
