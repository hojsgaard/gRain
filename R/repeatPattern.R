#' @title Create repeated patterns in Bayesian networks
#' 
#' @description Repeated patterns is a useful model specification
#'     short cut for Bayesian networks
#' 
#' 
#' @param plist A list of conditional probability tables. The variable
#'     names must have the form \code{name[i]} and the \code{i} will
#'     be substituted by the values given in \code{instances} below.
#' @param instances A vector of distinct integers
#' @param unlist If \code{FALSE} the result is a list in which each
#'     element is a copy of \code{plist} in which \code{name[i]} are
#'     substituted. If \code{TRUE} the result is the result of
#'     applying \code{unlist()}.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{grain}}, \code{\link{compileCPT}}
#' @references Søren Højsgaard (2012). Graphical Independence
#'     Networks with the gRain Package for R. Journal of Statistical
#'     Software, 46(10), 1-26.
#'     \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords utils
#' @examples
#' 
#' ## Specify hidden markov models. The x[i]'s are unobserved, the
#' ## y[i]'s can be observed.
#' 
#' yn <- c("yes","no")
#' 
#' ## Specify p(x0)
#' x.0 <- cptable(~x0, values=c(1,1), levels=yn)
#' 
#' ## Specify transition density
#' x.x <- cptable(~x[i]|x[i-1], values=c(1,99,2,98),levels=yn)
#' 
#' ## Specify emissiob density
#' y.x <- cptable(~y[i]|x[i],   values=c(1,99,2,98),levels=yn)
#' 
#' ## The pattern to be repeated
#' pp <- list(x.x, y.x)
#' 
#' ## Repeat pattern and create network
#' ppp <- repeatPattern(pp, instances=1:10)
#' qqq <- compileCPT(c(list(x.0),ppp))
#' rrr <- grain(qqq)
#' 
#' 
#' 
#' @export repeatPattern
repeatPattern <- function(plist, instances, unlist=TRUE){
  ans <- list()
  for (ii in seq_along(instances)){
    ans[[ii]] <- .do.one(plist, instances[[ii]])
  }
  if (unlist)
    ans <- unlist(ans, recursive=FALSE)

  ans
}


.subst <- function(x, i.val){
                                        #vv <- c("xyz[i+1]tyu", "xx[i]")
                                        #x <- c("xyz[i+1]tyu", "xx[i]","kkkk")
                                        #x <- c("xyztyu", "xx","kkkk")
  with.brack <- grep("\\[",x)
  vv <- x[with.brack]
  
  if (length(vv)>0){
    idx.vec <- gsub("[^\\[]*\\[([^\\]*)\\].*", "\\1", vv)
    idx.exp <- parse(text=idx.vec)
    idx.val <- unlist(lapply(idx.exp, eval, list(i=i.val)))
    vv2 <- list()
    for (ii in seq_along(idx.val)){
      vv2[[ii]] <- gsub("\\[([^\\]*)\\]", idx.val[ii], vv[ii])
    }
    vv2 <- unlist(vv2)
    x[with.brack] <- vv2
  }
  x
}

.do.one <- function(plist1, i.val){
  pp <- lapply(plist1, function(xx){
    xx$vpa <- .subst(xx$vpa,i.val)
    xx
  }) 
  pp 
}


