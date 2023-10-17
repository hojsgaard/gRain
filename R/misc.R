
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
