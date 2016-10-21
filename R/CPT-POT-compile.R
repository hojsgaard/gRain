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
#'     potentials
#' @param forceCheck Controls if consistency checks of the probability
#'     tables should be made.
#' @param details Controls amount of print out. Mainly for debugging
#'     purposes
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
compileCPT <- function(x, forceCheck=TRUE, details=0){

    ## if (class(x)[1] != "extractCPT")
    ##     stop("can not compile 'x'\n")
    #details=1

    zz <- lapply(x, .parseCPTlist)

    vnamList <- lapply(zz, "[[", "vnam") ## variable name
    vparList <- lapply(zz, "[[", "vpar") ## variable and parent names
    vlevList <- lapply(zz, "[[", "vlev") ## variable level
    vn       <- unlist(vnamList, use.names=FALSE)
    names( vlevList ) <- vn

    out        <- vector("list", length(vn))
    names(out) <- vn
    di         <- unlist( lapply(vlevList, length) )

    if (details>=1) cat(". creating probability tables ...\n")

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

            out[[ i ]] <-
                parray(varNames  = vpar,
                       levels    = lev,
                       values    = val,
                       normalize = zz[[ i ]]$normalize,
                       smooth    = zz[[ i ]]$smooth
                       )
        }
    }

    if (details>=1) cat(". creating dag and checking for acyclicity...\n")

    ## Check for acyclicity
    dg  <- dagList(vparList)
    ##dgM <- as(dg, "Matrix")
    oo  <- topoSort(dg)
    if (oo[1]==-1)
        stop("Graph defined by the cpt's is not acyclical...\n");

    ## FIXME: topoSort returnerer character(0) så der skal checkes på length==0
    
    universe        <- list(nodes = vn, levels = vlevList, nlev = di)
    attributes(out) <- list(universe = universe,
                            dag      = dg,
                            vparList = vparList,
                            names    = vn )

    class(out) <- "CPTspec"
    out
}

#' @rdname compile-cpt
compilePOT <- function(x){
    if (class(x)[1] != "extractPOT")
        stop("can not compile 'x'\n")
    uug <- ugList(lapply(x, function(a) names(dimnames(a))))
    if (!is.TUG(uug))
        stop("Graph defined by potentals is not triangulated...\n")

    lll       <- unlist(lapply(x, dimnames),recursive=FALSE)
    nnn       <- names(lll)
    iii       <- match(unique(nnn), nnn)
    levels    <- lll[iii]
    vn        <- nnn[iii]
    di        <- c(lapply(levels, length), recursive=TRUE)
    names(di) <- vn
    universe  <- list(nodes = vn, levels = levels, nlev   = di)
    ans       <- x
    attributes(ans) <- list(universe = universe,
                            ug       = uug,
                            rip      = attr(x, "rip"),
                            dag      = attr(x, "dag"),      ## Needed to save network in Hugin format
                            cptlist  = attr(x, "cptlist"))  ## Needed to save network in Hugin format

    class(ans) <- "POTspec"
    return(ans)
}

.parseCPTlist <- function(xi){ ## Create intermediate form of CPTs

    if (is.array(xi))
        cls <- "array"
    else
        cls <- class(xi)[1]
    
    vpar <- varNames(xi)
    vlev <- valueLabels(xi)[[1]]
    
  switch(cls,
         "cptable"={
             tmp     <- list(vnam=vpar[1], vlev=vlev, vpar=vpar, values=xi$values,
                             ##normalize="first", smooth=0)
                             normalize=if (xi$normalize) "first" else "none", smooth=xi$smooth)
         },
         "parray"={
             tmp     <- list(vnam=vpar[1], vlev=vlev, vpar=vpar, values=as.numeric(xi),
                             normalize="first", smooth=0)
         },
         "array"={
             tmp     <- list(vnam=vpar[1], vlev=vlev, vpar=vpar, values=as.numeric(xi),
                             normalize="first", smooth=0)
         })
  tmp
}


print.CPTspec <- function(x,...){
  cat("CPTspec with probabilities:\n")
  lapply(x,
         function(xx){
           vn <- varNames(xx)
           if (length(vn)>1){
             cat(paste(" P(",vn[1],"|",paste(vn[-1],collapse=' '),")\n"))
           } else {
             cat(paste(" P(",vn,")\n"))
           }
         })
  return(invisible(x))
}

summary.CPTspec <- function(object,...){
    lapply(object, function(x){x})
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
}





print.POTspec <- function(x,...){
  cat("POTspec with potentials:\n")
  lapply(x,
         function(xx){
           vn <- names(dimnames(xx))
           cat(   "(", paste(vn, collapse=' '),") \n")
         })

  return(invisible(x))
}


