## ###############################################################
##
## #' @title Evidence objects
## #' @description Functions for defining and manipulating evidence.
## #' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
## #' @name evidence_object
##
## ###############################################################
## '
## ' @aliases subset.grain_evidence print.grain_evidence 
## ' 
## ' @details Evidence is specified as a list. Internally, evidence is
## '     represented as a grain evidence object which is a list with 4 elements.
## ' 
## ' @examples
## '
## ' ## Define the universe
## ' yn <- c("yes", "no")
## ' uni <- list(asia = yn, tub = yn, smoke = yn, lung = yn,
## '             bronc = yn, either = yn, xray = yn, dysp = yn)
## '
## ' e1 <- list(dysp="no", xray="no")
## ' eo1 <- grain_evidence_new(e1, levels=uni)
## ' eo1  |> as.data.frame()
## ' 
## ' e2 <- list(dysp="no", xray=c(0, 1))
## ' eo2 <- grain_evidence_new(e2, levels=uni)
## ' eo2 |> as.data.frame()
## '
## ' # Above e1 and e2 specifies the same evidence but information about
## ' # whether the state has been set definite or as a weight is
## ' # maintained.
## ' 
## ' e3 <- list(dysp="yes", asia="yes")
## ' eo3 <- grain_evidence_new(e3, uni)
## ' eo3 |> as.data.frame()
## ' 
## ' # If evidence 'e1' is already set in the network and new evidence
## ' # 'e3' emerges, then evidence in the network must be updated. But
## ' # there is a conflict in that dysp="yes" in 'e1' and
## ' # dysp="no" in 'e3'. The (arbitrary) convention is that
## ' # existing evidence overrides new evidence so that the only new
## ' # evidence in 'e3' is really asia="yes".
## '
## ' # To subtract existing evidence from new evidence we can do:
## ' zz <- grain_evidence_setdiff(eo3, eo1)
## '
## ' # Likewise the 'union' is
## ' zz <- grain_evidence_union(eo3, eo1)
## '
## ' @export 
## ' @rdname evidence_object
## ' @param evi_list A named list with evidence; see 'examples' below.
## ' @param levels A named list with the levels of all variables. 

grain_evidence_new <- function(evi_list=NULL, levels) {

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
    #ooo <<- out
    #print("pppppppppppppppppppp\n")
    out <- grain_evidence2dataframe(out)
    class(out) <- c("grain_evidence", "data.frame")
    out
}

is.null_evi <- function(object) {
    if (missing(object)) TRUE
    else if (length(object) == 0) TRUE
    else if (inherits(object, "grain_evidence") && length(grain_evidence_names(object)) == 0) TRUE
    else FALSE

}

grain_evidence_names <- function(x) x$nodes

## ' @name evidence_object
## ' @param row.names Not used.
## ' @param optional Not used.
## ' @param x An evidence object.
## ' @param ... Not used.
## ' @export 

grain_evidence2dataframe <- function(x) {
  x <<-x
    mm <- lapply(x, function(z) as.data.frame(I(z)))
    mm <- as.data.frame(mm)
    names(mm) <- names(x)
    mm
}

grain_evidence2dataframe <-
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



grain_evidence_setdiff <- function(ev1, ev2) {
    if (length(ev1) == 0) ev1 <- grain_evidence_new( ev1 )
    if (length(ev2) == 0) ev2 <- grain_evidence_new( ev2 )
    
    nn  <- setdiff( grain_evidence_names(ev1), grain_evidence_names(ev2) )
    out <- grain_evidence_subset(ev1, select=nn)
    class(out) <- c("grain_evidence", "list")
    out <- grain_evidence2dataframe(out)
    class(out) <- c("grain_evidence", "data.frame")
    out

}

grain_evidence_union <- function(ev1, ev2) {
    if (length(ev1)==0) ev1 <- grain_evidence_new( ev1 )
    if (length(ev2)==0) ev2 <- grain_evidence_new( ev2 )
    ev <- grain_evidence_setdiff( ev1, ev2 )
    out <- mapply(function(l1, l2){c(l1,l2)},
                  ev, ev2, SIMPLIFY=FALSE, USE.NAMES=TRUE)
    class(out) <- c("grain_evidence", "list")
    out <- grain_evidence2dataframe(out)
    class(out) <- c("grain_evidence", "data.frame")
    out
}


grain_evidence_subset <- function(x, subset, select, drop = FALSE, ...){
    if (missing(select)) x
    else if (length(select)==0) grain_evidence_new(list())
    else {
        nl <- as.list(1L:length(grain_evidence_names(x)))
        names(nl) <- grain_evidence_names(x)
        vars <- eval(substitute(select), nl, parent.frame())
        if (is.character(vars))
            vars <- match(vars, grain_evidence_names(x))
        if (any(is.na(vars)))
            stop("'vars' contain NA")
        if (max(vars) > length(grain_evidence_names(x)))
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

## Bruges i grain_evidence_new
hard_state_to_array <- function(n, v, lev){
    #str(list(n,v,lev))
    tab <- fast_array(n, list(lev), rep.int(0, length(lev)))
    tab[ match( v, lev )] <- 1
    tab
}

## Bruges i grain_evidence_new
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
