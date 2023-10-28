#' @title Create repeated patterns in Bayesian networks
#' 
#' @description Repeated patterns is a useful model specification
#'     short cut for Bayesian networks
#'
#' @name repeat_pattern
#' 
#' @param plist A list of conditional probability tables. The variable
#'     names must have the form \code{name[i]} and the \code{i} will
#'     be substituted by the values given in \code{instances} below.
#'     See also the \code{data} argument.
#' @param instances A vector of consecutive integers
#' @param unlist If \code{FALSE} the result is a list in which each
#'     element is a copy of \code{plist} in which \code{name[i]} are
#'     substituted. If \code{TRUE} the result is the result of
#'     applying \code{unlist()}.
#' @param data A two column matrix. The first column is the index /
#'     name of a node; the second column is the index / name of the
#'     node's parent.
#' 
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{grain}}, \code{\link{compile_cpt}}
#' @references Søren Højsgaard (2012). Graphical Independence
#'     Networks with the gRain Package for R. Journal of Statistical
#'     Software, 46(10), 1-26.
#'     \url{https://www.jstatsoft.org/v46/i10/}.
#' @keywords models
#' @examples
#'
#' yn <- c("yes", "no")
#' n <- 3
#'
#' ## Example: Markov chain
#' 
#' x_init  <- cpt(~x0, values=c(1, 9), levels=yn)                  ## p(x0)
#' x_trans <- cpt(~x[i]|x[i-1], values=c(1, 99, 2, 98), levels=yn) ## p(x[i]|x[i-1])
#' pat     <- list(x_trans)                             
#' rep.pat <- repeat_pattern(pat, instances=1:n)
#'
#' mc <- compile_cpt(c(list(x_init), rep.pat))
#' mc
#' mc <- mc |> grain()
#' 
#' ## Example: Hidden markov model:
#' # The x[i]'s are unobserved, the y[i]'s can be observed.
#' 
#' x_init  <- cpt(~x0, values=c(1, 9), levels=yn)                   ##  p(x0)
#' x_trans <- cpt(~x[i]|x[i-1], values=c(1, 99, 2, 98), levels=yn)  ##  p(x[i]|x[i-1])
#' y_emis  <- cpt(~y[i]|x[i], values=c(10, 90, 20, 80), levels=yn)  ##  p(y[i]|x[i]) 
#'
#' pat     <- list(x_trans, y_emis) ## Pattern to be repeated
#' rep.pat <- repeat_pattern(pat, instances=1:n)
#' hmm <- compile_cpt(c(list(x_init), rep.pat)) 
#' hmm
#' hmm <- hmm |> grain()
#' 
#' ## Data-driven variable names
#' 
#' dep <- data.frame(i=c(1, 2, 3, 4, 5, 6, 7, 8),
#'                   p=c(0, 1, 2, 2, 3, 3, 4, 4))
#'
#' x0 <- cpt(~x0, values=c(0.5, 0.5), levels=yn)
#' xa <- cpt(~x[i] | x[data[i, "p"]], values=c(1, 9, 2, 8), levels=yn)
#' xb <- repeat_pattern(list(xa), instances=1:nrow(dep), data=dep)
#' tree <- compile_cpt(c(list(x0), xb))
#' tree
#' tree <- tree |> grain()
#' tree 
#' 


#' @rdname repeat_pattern
#' @export 
repeat_pattern <- function(plist, instances, unlist=TRUE, data=NULL){

    ans <- vector("list", length(instances))
    for (i in seq_along(instances)){
        ans[[ i ]] <- do_pattern_instance(plist, instances[[ i ]], data)
    }
    if (unlist)
        ans <- unlist(ans, recursive=FALSE)
    ans
}


#' @rdname repeat_pattern
#' @export 
repeatPattern <- repeat_pattern

do_pattern_instance <- function(pat_list1, i.val, data=NULL){
    pp <- lapply(pat_list1, function(xx){
        if (inherits(xx, "array")){
            set_dim_names_array(xx, i.val, data)
        } else
            if (inherits(xx, "cptable_class")){
                set_dim_names_cptable(xx, i.val, data)
            } 
        else {
                stop("no behaviour defined\n")
            }
    }
    )
    pp 
}


set_dim_names_array <- function(xx, i, data=NULL){
    nms <- names(dimnames(xx))
    head <- gsub("\\[.*", "", nms)

    if (is.null(data)){
        b1 <- gsub("^.*\\[(.*)\\]", "\\1", nms)
        b2 <- eval(parse(text=paste0("c(", toString(b1), ")")), list(i=i))
    } else {
        b2 <- data[i, ]
    }
    nms2 <- paste0(head, b2)
    names(dimnames(xx)) <- nms2
    xx            
}

set_dim_names_cptable <- function(xx, i, data=NULL){
    nms <- attr(xx, "vpa")
    head <- gsub("\\[.*", "", nms)

    if (is.null(data)){
        b1 <- gsub("^.*\\[(.*)\\]", "\\1", nms)
        b2 <- eval(parse(text=paste0("c(", toString(b1), ")")), list(i=i))
    } else {
        b2 <- data[i, ]
    }
    nms2 <- paste0(head, b2)
    attr(xx, "vpa") <- nms2
    xx            
}



## ' ## All examples from above using cptable
## ' 
## ' z0   <- cptable(~z0, values=c(1, 9), levels=yn)                    ## p(z0)
## ' z_z <- cptable(~z[i]|z[i-1], values=c(1, 99, 2, 98), levels=yn)   ## p(z[i]|z[i-1])
## ' pat  <- list(z_z) ## Pattern to be repeated
## ' rep.pat <- repeat_pattern(pat, instances=1:n)
## ' mc <- compile_cpt(c(list(z0), rep.pat))
## ' mc 
## ' mc |> grain()
## '
## ' #' z0   <- cptable(~z0, values=c(1, 9), levels=yn)                    ##  p(z0)
## ' z_z <- cptable(~z[i]|z[i-1], values=c(1, 99, 2, 98), levels=yn)   ##  p(z[i]|z[i-1])
## ' u_z  <- cptable(~u[i]|z[i], values=c(10, 90, 20, 80), levels=yn)   ##  p(u[i]|z[i])
## '
## ' pat  <- list(z_z, u_z) ## Pattern to be repeated
## ' rep.pat <- repeat_pattern(pat, instances=1:n)
## ' hmm <- compile_cpt(c(list(z0), rep.pat))
## ' hmm
## ' hmm |> grain()
## '
## '
## ' #' z0 <- cptable(~z0, values=c(0.5, 0.5), levels=yn)
## ' za <- cptable(~z[i] | z[data[i, "p"]], values=c(0.5, 0.5), levels=yn)
## ' zb <- repeat_pattern(list(za), instances=1:nrow(dep), data=dep)
## ' tree <- compile_cpt(c(list(z0), zb))  |> grain()
## ' 
