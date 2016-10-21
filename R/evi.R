## ###############################################################
##
#' @title Set evidence in grain objects
#'
#' @description Setting and removing evidence in grain objects. 
#'
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' 
## ###############################################################
#'
#' @name grain-evi
#'
#' @param object A "grain" object
#' @param nodes A vector of nodes; those nodes for which the
#'     (conditional) distribution is requested.
#' @param states A vector of states (of the nodes given by 'nodes')
#' @param evidence An alternative way of specifying findings
#'     (evidence), see examples below.
#' @param propagate Should the network be propagated?
#' @param details Debugging information
#'
#' @seealso \code{\link{setEvidence}} \code{\link{getEvidence}}
#'     \code{\link{retractEvidence}} \code{\link{pEvidence}}
#'     \code{\link{setFinding}} \code{\link{getFinding}}
#'     \code{\link{retractFinding}} \code{\link{pFinding}}
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords models utilities
#' @examples
#'
#' ## setFinding / setEvidence
#' 
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
#' bn <- grain(plist)
#'  
#' ## 1) These forms are identical
#'
#' e1 <- list(dysp="no", xray="no")
#' 
#' setEvi(bn, evidence=e1)
#' setEvi(bn, nodes=c("dysp","xray"), states=c("no", "no"))
#' setEvidence(bn, nodes=c("dysp","xray"), states=c("no", "no"))
#'
#' # Notice: setFinding is old school but it was used in the
#' # "Graphical Models with R" book.
#' setFinding(bn, nodes=c("dysp","xray"), states=c("no", "no"))
#'
#' ## 2) Updating evidence
#' # Notice that only 'asia' is set because 'dysp' was set earlier
#' 
#' e2 <- list(dysp="yes", asia="yes")
#' bn1 <- setEvi(bn, evidence=e1)
#' bn1
#' bn2 <- setEvi(bn1, evidence=e2)
#' bn2
#'
#' ## 3) Shorter forms
#'
#' bn2 <- bn
#' evidence(bn2)
#' evidence(bn2) <- e1
#' evidence(bn2)
#' evidence(bn2) <- e2
#' evidence(bn2)
#' evidence(bn2) <- NULL
#' evidence(bn2)
#'
#'
#' ## 4) Alternative forms:
#' 
#' setEvi(bn, evidence=list("asia"=c(1, 0), "xray"="yes"))
#'
#' ## 5) Suppose we do not know with certainty whether a patient has
#' ## recently been to Asia. We can then introduce a new variable
#' ## "guess.asia" with "asia" as its only parent. Suppose
#' ## p(guess.asia=yes|asia=yes)=.8 and p(guess.asia=yes|asia=no)=.1
#' ## If the patient is e.g. unusually tanned we may set
#' ## guess.asia=yes and propagate. This corresponds to modifying the
#' ## model by the likelihood (0.8, 0.1) as
#'
#' b =setEvi(bn, nodes=c("asia","xray"), states=list(c(0.8,0.1), "yes"))
#' as.data.frame( evidence( b ) )
#'



#' @rdname grain-evi
setEvi <- function(object, nodes=NULL, states=NULL, evidence=NULL, 
                   propagate=TRUE, details=0){
    if ( !inherits(object, "grain") )
        stop("'object' is not a 'grain' object")
    
    if ( is.null( evidence ) && is.null( nodes ) )
        stop( "Evidence is not given; nothing to do...")
    
    if ( is.null( evidence ) ) ## Then 'nodes' must be given
        evidence <- .nodes.states2evidence( nodes, states )      
    
    setEvi_( object, evidence=evidence, propagate=propagate, details=details)
}

## setEvi_: New implementation; august 2016
#' @rdname grain-evi
setEvi_ <- function(object, evidence=NULL, propagate=TRUE, details=0){
#    details=0
#    cat("++++ setEvi_ input evidence: \n"); str(evidence)

    ## If object is not compiled then do so
    
    if ( is.null_ev( evidence ) ){
        cat("Nothing to do\n")
    } else {
        if (!object$isCompiled){
            object <- compile(object)
            object$isInitialized  <- FALSE
            object$isPropagated   <- FALSE            
        } 
        
        oe <- getEvidence( object ) # Der er noget med typen (old/new)
        if (details>0){ print("old.evidence"); print(oe) }
        ne <- new_ev( evidence, universe(object)$levels ) 
        if (details>0){ print("new evidence"); print( ne ) }
        
        ## Hvis der er eksisterende evidens, så skal det tages ud af det nye
        if (!is.null_ev( oe ) ){
            ne <- setdiff_ev( ne, oe )
            if (details>0) { print("new evidence - after modif"); print( ne ) }
        }

        if ( length( varNames( ne ) ) > 0 ){
            # host  <- .get.host.clique( ne$evidence, object$rip )
            # object$temppot <- .insert.evidence.in.potential( object$temppot, ne, host )
            # str(list(host=host))
            ## FIXME: 3 lines: Replaces .get.host.clique
            rp  <- getgrain(object, "rip")    
            host  <- getHostClique(varNames( ne ), rp$cliques)
            #str(list(ne, host))
            #pot.b <<- object$temppot
            object$temppot <- insertEvi( ne, object$temppot, host )
            #pot.a <<- object$temppot
            
            te <- if (is.null_ev( oe )) ne else union_ev( oe, ne )
            ##te <- union_ev( oe, ne )
            if (details>0) {print("total evidence"); print( te )}
            object$evidence <- te
        }         
    } 
    if (propagate) propagate(object) else object
}


#' @rdname grain-evi
#' @param items Items in the evidence list to be removed. Here,
#'     \code{NULL} means remove everything. If \code{items} is a
#'     character vector (of nodes) then evidence on these nodes is
#'     removed. If \code{items} is a numeric vector then those items
#'     in the evidence list is removed. Notice that \code{0} means nothing is
#'     removed.
retractEvi <- function(object, items=NULL, propagate=TRUE){
    if ( !inherits(object, "grain") )
        stop("'object' is not a 'grain' object")
    retractEvi_(object, items=items, propagate=propagate)
}

#' @rdname grain-evi
retractEvi_ <- function(object, items=NULL, propagate=TRUE){
    ##cat("++++ retractEvidence_\n")
    .resetgrain <- function(x){
        x$temppot       <- x$origpot
        x$evidence       <- NULL
        x$isPropagated  <- FALSE
        x
    }

    if ( is.null( items ) ){
        object <- .resetgrain( object )
    } else {
        if (!(is.character( items ) || is.numeric( items )) )
            stop("'items' must be a character or numeric vector")
        oe <- getEvidence( object )
        if (!is.null_ev( oe )){
            vn <- varNames( oe )
            if (is.numeric(items)){
                items <- vn[items]
            }
            keep <- setdiff( vn, items)
            ## NB: keep <= varNames            
            ## hvis keep==varNames(oe) så gør intet
            ## hvis keep=Ø så bare reset
            ## hvis Ø < keep < varNames så gør som nedenfor
            if ( length( keep ) < length( vn ) ){
                object <- .resetgrain( object )                
                if (length( keep ) > 0){
                    ne <- subset( oe, select=keep )
                    object <- setEvi( object, evidence = ne )
                }
            }
        }
    }
    if (propagate) propagate(object) else object
}

#' @name grain-evi
absorbEvi <- function(object, propagate=TRUE ){
    if ( !inherits(object, "grain") )
        stop("'object' is not a 'grain' object")
    absorbEvi_( object, propagate = propagate )
}

#' @name grain-evi
absorbEvi_<- function(object, propagate=TRUE ){
    ## Update 'object' as
    ## 1) set origpot <- temppot
    ## 2) ignore any finding
    object$origpot <- object$temppot
    object$evidence <- NULL
    object$isPropagated <- FALSE

    if (propagate) propagate(object) else object
}



#' @name grain-evi
pEvidence <- function(object){
    if ( !inherits(object, "grain") )
        stop("'object' is not a 'grain' object")
    attr(object$equipot,"pEvidence")
}

#' @name grain-evi
getEvidence <- function(object){
    if ( !inherits(object, "grain") )
        stop("'object' is not a 'grain' object")
    object$evidence
}

#' @name grain-evi
"dropEvi<-" <- function(object, value=NULL){
    retractEvi(object, value)
}

#' @name grain-evi
"addEvi<-" <- function(object, value=NULL){
    setEvi(object, value)
}

#' @name grain-evi
evidence <- function(object){
    UseMethod("evidence")
    
}

#' @name grain-evi
evidence.grain <- function(object){
    getEvidence( object )
}

#' @name grain-evi
#' @param value The evidence in the form of a named list or an evidence-object. 
"evidence<-" <- function(object, value=NULL){
    UseMethod("evidence<-")
}

#' @name grain-evi
"evidence<-.grain" <- function(object, value=NULL){
    if (is.null( value ))
        retractEvi( object )
    else
        setEvidence(object, evidence=value)
}

## Alternative names

#' @name grain-evi
addEvi  <- setEvi

#' @name grain-evi
dropEvi <- retractEvi

#' @name grain-evi
getEvi  <- getEvidence

## #############################################################
##
## UTILITIES
##
## #############################################################


## Bruges af setEvi; FIXME: lav om så vi bruger synlig funktion                     
## .insert.evidence.in.potential <- function( pot, evi.list, hostclique ){
##     ## if (any(is.na(hostclique))) stop("NAs in hostclique...")
##     for (i in seq_along( evi.list$evidence ) ){
##         j <- hostclique[ i ]
##         p <- evi.list$evidence[[ i ]]
##         pot[[j]] <- tabMult__( pot[[ j ]], p )
##     }
##     pot
## }

#' @rdname grain-evi
#' @param evi.list A "grain_ev" object.
#' @param pot A list of clique potentials (a potential is an array).
#' @param hostclique A numerical vector indicating in which element of
#'     'pot' each eviendence item in 'evi.list' should be inserted in.
insertEvi <- function(evi.list, pot, hostclique){
    ##cat("insertEvi\n")
    if ( !inherits(evi.list, "grain_ev") )
        stop("'object' is not a 'grain_ev' object")
    
    for (i in seq_along( evi.list$evidence ) ){
        p <- evi.list$evidence[[ i ]]
        #print(p)
        j <- hostclique[ i ]
        #print( j )
        pot[[j]] <- tabMult__( pot[[ j ]], p )
    }
    pot
}





.nodes.states2evidence <- function(nodes, states){
    if (!is.null( states ) && length( nodes )==length( states )){
        evidence <- as.vector(states, "list")
        names(evidence) <- nodes
        evidence
    } else {
        stop( "Either states are not given or nodes and states do not have same length" )
    }
}

## Bruges af setEvi; FIXME: lav om så vi bruger synlig funktion
## .get.host.clique <- function(evidence, rip){
##     unlist(lapply(evidence,
##                   function(x){
##                       n <- names(dimnames(x))
##                       gRbase::get_superset_(n, rip$cliques, all=FALSE)
##                   }),
##                   use.names = FALSE)
## }

#' @rdname grain-evi
#' @param set.list A list of sets (a set is a character vector).
#' @param cliques A list of sets (a set is a character vector).
getHostClique <- function(set.list, cliques){
    out <- lapply(set.list,
                  function(x){
                      gRbase::get_superset_( x, cliques, all=FALSE)
                  })
    len <- sapply(out, length)
    if (any( len == 0)){
        message("Set(s) not in any clique:")
        str(set.list[ len==0 ])
        stop("exiting...\n")
    }
    unlist( out )
}








