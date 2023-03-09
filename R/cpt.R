#' Create conditional probability table CPT.
#'
#' @param names Names of variables defining table; either a character
#'     vector or a right hand sided formula.
#' @param levels 1) a list with specification of the levels of the
#'     factors in \code{names} or 2) a vector with number of levels of
#'     the factors in \code{names}. See 'examples' below.
#' @param values values to go into the array.
#' @param smooth Should values be smoothed, see 'Details' below.
#' @return An array.
#' @keywords utilities
#' @examples
#' 
#' universe <- list(gender=c('male', 'female'),
#'                  answer=c('yes', 'no'),
#'                  rain=c('yes', 'no'))
#' t1 <- cpt(c("gender", "answer"), levels=universe, values=1:4)
#' t1
#' t2 <- cpt(~gender:answer, levels=universe, values=1:4)
#' t2
#' t3 <- cpt(~gender:answer, c(2, 2), values=1:4)
#' t3
#' @export
cpt <- function(names, levels, values, smooth=0){
    names <- c(.formula2char(names))
    tabNew(names=names, levels=levels, values=values, normalize = "first", smooth=smooth)
}

#' Simplify output query to a Bayesian network
#'
#' Simplify output query to a Bayesian network to a dataframe provided
#' that each node has the same levels.
#'
#' @param b Result from running querygrain.
#' 
#' @export
simplify_query <- function(b){
    if (!inherits(b, "list")){
        return(b)
    }
    nms <- lapply(b, names)
    v <- sapply(nms, function(n) {all.equal(n, nms[[1]])})
    
    if (all(v)){
        out <- as.data.frame(do.call(rbind, lapply(b, as.numeric)))
        names(out) <- nms[[1]]
        return(out)
    } else {
        return(b)
    }
}


#' @export
tidy.grain_evidence <- function(x, ...){
    as.data.frame(x)
}
