## ###############################################################
##
#' @title Evidence objects
#' @description Functions for defining and manipulating evidence.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @name evidence_object
##
## ###############################################################
#'
#' @aliases subset.grain_evidence print.grain_evidence varNames.grain_evidence
#' 
#' @details Evidence is specified as a list. Internally, evidence is
#'     represented as a grain evidence object which is a list with 4 elements.
#' 
#' @examples
#'
#' ## Define the universe
#' yn <- c("yes", "no")
#' uni <- list(asia = yn, tub = yn, smoke = yn, lung = yn,
#'             bronc = yn, either = yn, xray = yn, dysp = yn)
#'
#' e1 <- list(dysp="no", xray="no")
#' eo1 <- new_evi(e1, levels=uni)
#' eo1  |> as.data.frame()
#' 
#' e2 <- list(dysp="no", xray=c(0, 1))
#' eo2 <- new_evi(e2, levels=uni)
#' eo2 |> as.data.frame()
#'
#' # Above e1 and e2 specifies the same evidence but information about
#' # whether the state has been set definite or as a weight is
#' # maintained.
#' 
#' e3 <- list(dysp="yes", asia="yes")
#' eo3 <- new_evi(e3, uni)
#' eo3 |> as.data.frame()
#' 
#' # If evidence 'e1' is already set in the network and new evidence
#' # 'e3' emerges, then evidence in the network must be updated. But
#' # there is a conflict in that dysp="yes" in 'e1' and
#' # dysp="no" in 'e3'. The (arbitrary) convention is that
#' # existing evidence overrides new evidence so that the only new
#' # evidence in 'e3' is really asia="yes".
#'
#' # To subtract existing evidence from new evidence we can do:
#' setdiff_evi(eo3, eo1) |> as.data.frame()
#'
#' # Likewise the 'union' is
#' union_evi(eo3, eo1) |> as.data.frame()
#'
#' @export 
#' @rdname evidence_object
#' @param evi_list A named list with evidence; see 'examples' below.
#' @param levels A named list with the levels of all variables. 
new_evi <- function(evi_list=NULL, levels){

    if (inherits(evi_list, "grain_evidence")) {
        return(evi_list)
    }

    if (length(evi_list) == 0){
        out <- list(nodes      = character(0),
                    is_hard    = logical(0),
                    hard_state = character(0),
                    evi_weight = list() )
    } else {
        ## First remove all evidence specified as NA
        not.na <- !unlist(lapply(lapply(evi_list,is.na), any), use.names=FALSE)
        if (length( not.na ) > 0)
            evi_list <- evi_list[ not.na ]
        
        evi_weight   <- vector("list", length(evi_list))
        is_hard      <- rep.int(TRUE,  length(evi_list))
        hard_state   <- rep.int(NA,    length(evi_list))

        
        for (i in seq_along(evi_list)){
            ev <- evi_list[i]
            v <- ev[[1]]
            
            if( is.array(v)){
                n <- names(dimnames(v))
                is_hard[i]    <- FALSE
                evi_weight[[i]] <- v
                next
            }
            
            if (is.character(v)){
                n <- names(evi_list)[i]
                hard_state[i]  <- v
                evi_weight[[i]]  <- hard_state_to_array(n, v, levels[[n]])
                next
            }
            
            if (is.numeric(v)){
                n <- names(evi_list)[i]
                is_hard[i] <- FALSE
                evi_weight[[i]] <- soft_state_to_array(n, v, levels[[n]])
            }
        }
        
        ## If evidence is zero on all states or negative on some (or all) states then it is invalid
        keep <- unlist(lapply(evi_weight, function(e){ sum(e) !=0 && all(e>=0) }), use.names=FALSE)
        ## print(keep)
        
        nodes <- unique.default(unlist(lapply(evi_weight, .namesDimnames)),
                                use.names=FALSE )
        out <- list(
            nodes      = nodes[keep],
            is_hard    = is_hard[keep],
            hard_state = hard_state[keep],
            evi_weight = evi_weight[keep])
    }
    class(out) <- c("grain_evidence", "list")
    out
}

#' @export 
#' @rdname evidence_object
#' @param object Some R object.
is.null_evi <- function(object){
    if (missing(object)) TRUE
    else if (length(object) == 0) TRUE
    else if (inherits(object, "grain_evidence") && length(varNames(object)) == 0) TRUE
    else FALSE

}

## #' @rdname evidence_object
## #' @param x Evidence object


## #' @rdname evidence_object
#' @export
varNames.grain_evidence <- function(x) x$nodes

#' @rdname evidence_object
#' @param row.names Not used.
#' @param optional Not used.
#' @param x An evidence object.
#' @param ... Not used.
#' @export 
as.data.frame.grain_evidence <-
    function (x, row.names = NULL, optional = FALSE, ...) {
        is.atom <- sapply(x, is.atomic)
        atom <- x[is.atom]
        n.atom <- length(atom)
        out <- as.data.frame(atom)
        notatom <- x[!is.atom]
        n <- names(notatom)
        for (i in 1:length(notatom)){
            out[i+n.atom] <- notatom[i]
        }
        out
    }

#' @export 
#' @rdname evidence_object
#' @param ev1,ev2 Evidence.
setdiff_evi <- function(ev1, ev2){
    if (length(ev1) == 0) ev1 <- new_evi( ev1 )
    if (length(ev2) == 0) ev2 <- new_evi( ev2 )
    
    nn  <- setdiff( varNames(ev1), varNames(ev2) )
    out <- subset(ev1, select=nn)
    class(out) <- c("grain_evidence", "list")
    out
}

#' @export 
#' @rdname evidence_object
union_evi <- function(ev1, ev2){
    if (length(ev1)==0) ev1 <- new_evi( ev1 )
    if (length(ev2)==0) ev2 <- new_evi( ev2 )
    ev <- setdiff_evi( ev1, ev2 )
    out <- mapply(function(l1, l2){c(l1,l2)},
                  ev, ev2, SIMPLIFY=FALSE, USE.NAMES=TRUE)
    class(out) <- c("grain_evidence", "list")
    out
}

#' @export 
subset.grain_evidence <- function(x, subset, select, drop = FALSE, ...){
    if (missing(select)) x
    else if (length(select)==0) new_evi(list())
    else {
        nl <- as.list(1L:length(varNames(x)))
        names(nl) <- varNames(x)
        vars <- eval(substitute(select), nl, parent.frame())
        if (is.character(vars))
            vars <- match(vars, varNames(x))
        if (any(is.na(vars)))
            stop("'vars' contain NA")
        if (max(vars) > length(varNames(x)))
            stop("'vars' too large")
        out <- lapply(x, "[", vars)
        class(out) <- class(x)
        out
    }
}


## ###############################################
##
## UTILITIES
##
## ###############################################

## Bruges i new_evi
hard_state_to_array <- function(n, v, lev){
    #str(list(n,v,lev))
    tab <- fast_array(n, list(lev), rep.int(0, length(lev)))
    tab[ match( v, lev )] <- 1
    tab
}

## Bruges i new_evi
soft_state_to_array <- function(n, v, lev){
    #str(list(n,v,lev))
    fast_array(n, list(lev), v)
}

## Bruges af ovenstående fns
fast_array <- function(varNames, levels, values=1){
    #str(list(varNames, levels, values))
    dn <- if(is.list(levels)) levels else list(levels)
    names(dn) <- varNames
    array(values, dimnames=dn)
}
