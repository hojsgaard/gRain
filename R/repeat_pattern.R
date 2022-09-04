#' @title Create repeated patterns in Bayesian networks
#' 
#' @description Repeated patterns is a useful model specification
#'     short cut for Bayesian networks
#' 
#' @param plist A list of conditional probability tables. The variable
#'     names must have the form \code{name[i]} and the \code{i} will
#'     be substituted by the values given in \code{instances} below.
#'     See also the \code{data} argument.
#' @param instances A vector of distinct integers
#' @param unlist If \code{FALSE} the result is a list in which each
#'     element is a copy of \code{plist} in which \code{name[i]} are
#'     substituted. If \code{TRUE} the result is the result of
#'     applying \code{unlist()}.
#' @param data Enable variable names of the form \code{name[data[i]]} - to
#'     enable data-driven variable names.
#' 
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{grain}}, \code{\link{compileCPT}}
#' @references Søren Højsgaard (2012). Graphical Independence
#'     Networks with the gRain Package for R. Journal of Statistical
#'     Software, 46(10), 1-26.
#'     \url{https://www.jstatsoft.org/v46/i10/}.
#' @keywords models
#' @examples
#'
#' ## Example: Markov chain
#' yn <- c("yes", "no")
#' 
#' ## Specify p(x0)
#' x.0 <- cptable(~x0, values=c(1, 9), levels=yn)
#' 
#' ## Specify transition density
#' x.x <- cptable(~x[i]|x[i-1], values=c(1, 99, 2, 98), levels=yn)
#'
#' ## Pattern to be repeated
#' pat <- list(x.x)
#' 
#' rep.pat <- repeatPattern(pat, instances=1:5)
#' cpt <- compileCPT(c(list(x.0), rep.pat))
#' mc <- grain(cpt)
#'
#' if (interactive()) iplot(mc)
#' 
#' ## Example: Hidden markov model: The x[i]'s are unobserved, the
#' ## y[i]'s can be observed.
#' 
#' yn <- c("yes", "no")
#' 
#' ## Specify p(x0)
#' x.0 <- cptable(~x0, values=c(1, 9), levels=yn)
#' 
#' ## Specify transition density
#' x.x <- cptable(~x[i]|x[i-1], values=c(1, 99, 2, 98), levels=yn)
#' 
#' ## Specify emission density
#' y.x <- cptable(~y[i]|x[i],   values=c(10, 90, 20, 80), levels=yn)
#' 
#' ## The pattern to be repeated
#' pat <- list(x.x, y.x)
#' 
#' ## Repeat pattern and create network
#' rep.pat <- repeatPattern(pat, instances=1:5)
#' cpt <- compileCPT(c(list(x.0), rep.pat))
#' hmm <- grain(cpt)
#' hmm 
#'
#' if (interactive()) iplot(hmm)
#' 
#' ## Data-driven variable names
#' x0 <- cptable(~x0, values=c(0.5, 0.5), levels=yn)
#' x <- cptable(~x[i] | x[data[i, "p"]], values=c(0.5, 0.5), levels=yn)
#' dep <- data.frame(i=c(1, 2, 3, 4, 5, 6, 7, 8),
#'                   p=c(0, 1, 2, 2, 3, 3, 4, 4))
#' x <- repeatPattern(list(x), instances=dep$i, data=dep)
#' tree <- compileCPT(c(list(x0), x))
#' tree <- grain(tree)
#' tree 
#' 
#' #' if (interactive()) iplot(tree)
#' 
#' @export repeatPattern
#' 
repeatPattern <- function(plist, instances, unlist=TRUE, data=NULL){
    ans <- vector("list", length(instances))
    for (i in seq_along(instances)){
        ans[[ i ]] <- .do.one(plist, instances[[ i ]], data)
    }
    if (unlist)
        ans <- unlist(ans, recursive=FALSE)
    ans
}

.do.one <- function(plist1, i.val, data=NULL){
    pp <- lapply(plist1, function(xx){
        ##xx$vpa <- .subst(xx$vpa, i.val)                       ## FIX 14/8/2017
        attr(xx, "vpa") <- .subst(attr(xx, "vpa"), i.val, data) ## FIX 14/8/2017, 24/7/2022
        xx
    }) 
    pp 
}

.subst <- function(x, i.val, data=NULL){
    ##vv <- c("xyz[i+1]tyu", "xx[i]")
    ##x <- c("xyz[i+1]tyu", "xx[i]","kkkk")
    ##x <- c("xyztyu", "xx","kkkk")
    with.brack <- grep("\\[",x)
    vv <- x[ with.brack ]
    
    if (length(vv)>0){
        idx.vec <- gsub("[^\\[]*\\[([^\\]*)\\].*", "\\1", vv)
        idx.exp <- parse(text=idx.vec)
        idx.val <- unlist(lapply(idx.exp, eval, list(i=i.val, data=data)))
        vv2 <- list()
        for (ii in seq_along(idx.val)){
            vv2[[ii]] <- gsub("\\[([^\\]*)\\]", idx.val[ii], vv[ii])
        }
        vv2 <- unlist(vv2)
        x[ with.brack ] <- vv2
    }
    x
}
