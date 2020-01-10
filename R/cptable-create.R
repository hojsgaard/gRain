#' @title Create conditional probability tables (CPTs)
#' 
#' @description Creates conditional probability tables of the form
#'     p(v|pa(v)).
#' 
#' @param vpar Specifications of the names in P(v|pa1,...pak). See
#'     section 'details' for information about the form of the
#'     argument.
#' @param values Probabilities; recycled if necessary. Regarding the
#'     order, please see section 'details' and the examples.
#' @param normalize See 'details' below.
#' @param smooth See 'details' below.
#' @param levels See 'details' below.
#'
#' @details
#' 
#' If \code{normalize=TRUE} then for each configuration of the parents
#' the probabilities are normalized to sum to one.
#' 
#' If \code{smooth} is non--zero then zero entries of \code{values} are
#' replaced with \code{smooth} before normalization takes place.
#' 
#' Regarding the form of the argument \code{vpar}: To specify \eqn{P(a|b,c)}
#' one may write \code{~a|b:c}, \code{~a:b:c}, \code{~a|b+c}, \code{~a+b+c} or
#' \code{c("a","b","c")}. Internally, the last form is used. Notice that the
#' \code{+} and \code{:} operator is used as a separator only. The order of the
#' variables IS important so the operators DO NOT commute.
#' 
#' If \code{a} has levels \code{a1,a2} and likewise for \code{b} and \code{c}
#' then the order of \code{values} corresponds to the configurations
#' \code{(a1,b1,c1)}, \code{(a2,b1,c1)} \code{(a1,b2,c1)}, \code{(a2,b2,c1)}
#' etc. That is, the first variable varies fastest.  Hence the first two
#' elements in \code{values} will be the conditional probabilities of \code{a}
#' given \code{b=b1, c=c1}.
#' 
#' @return A \code{cptable} object (a list).
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{andtable}}, \code{\link{ortable}},
#'     \code{\link{extract_cpt}}, \code{\link{compile_cpt}},
#'     \code{\link{extract_pot}}, \code{\link{compile_pot}},
#'     \code{\link{grain}}
#' @references Søren Højsgaard (2012). Graphical Independence
#'     Networks with the gRain Package for R. Journal of Statistical
#'     Software, 46(10), 1-26.
#'     \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords models
#' @examples
#' 
#' 
#' yn   <- c("yes", "no")
#' ynm  <- c("yes", "no", "maybe")
#' a    <- cptable(~ asia, values=c(1, 99), levels=yn)
#' t.a  <- cptable(~ tub : asia, values=c(5, 95, 1, 99, 1, 999),  levels=ynm)
#' d.a  <- cptable(~ dia : asia, values=c(5, 5, 1, 99, 100, 999), levels=ynm)
#' cptlist <- compileCPT(list(a, t.a, d.a))
#' grain(cptlist)
#' 
#' ## Example: Specifying conditional probabilities as a matrix
#' bayes.levels  <- c('Enzyme', 'Keratine', 'unknown')
#' root.node     <- cptable(~ R, values=c( 1, 1, 1 ), levels=bayes.levels)
#' cond.prob.tbl <- t(matrix(c(1, 0, 0, 0, 1, 0, 0.5, 0.5, 0),
#'    nrow=3, ncol=3, byrow=TRUE, dimnames=list(bayes.levels, bayes.levels)))
#' cond.prob.tbl
#' 
#' ## Notice above: Columns represent parent states; rows represent child states
#' query.node    <- cptable(~ Q | R, values=cond.prob.tbl, levels=bayes.levels)
#' sister.node   <- cptable(~ S | R, values=cond.prob.tbl, levels=bayes.levels)
#'
#' ## Testing 
#' compile(grain(compileCPT(list(root.node, query.node, sister.node))), propagate=TRUE)
#' 
#' @export cptable
#'

cptable <- function(vpar, levels=NULL, values=NULL, normalize=TRUE,  smooth=0 ){
    vpa  <- c(.formula2char(vpar))        
    if (is.list(levels))
        levels <- levels[[vpa[1]]]
    ## if (is.null(values))
    ##     values <- rep(1.0, length(levels))
    out  <- values
    attributes(out) <-
        list(vpa=vpa, normalize=normalize,
             smooth=smooth, levels=levels)
    class(out) <- "cptable_class"
    out
}


#' @rdname cptable
cptab <- cptable


## NORMAL

cdist <- function(vpar, parm=list()){
    vpa  <- c(.formula2char(vpar))        
    out <- parm
    names(out) <- c("intercept", "slope", "sd")
    attr(out, "vpa") <- vpa
    class(out) <- "cdist"
    out
}



print.cptable_class <- function(x, ...){
    ## "print.cptable\n" %>% cat
    v <- c(x)
    dim(v) <- c(length(attr(x,"levels")), length(v) / length(attr(x, "levels")))
    rownames(v) <- attr(x, "levels")
    cat(sprintf("{v, pa(v)} :\n"))
    str(attr(x, "vpa"))
    print(v)    
    ##str(attributes(x))
    invisible(x)
}


summary.cptable_class <- function(object, ...){
    print(object)
    str(attributes(object))
    invisible(object)
}

varNames.cptable_class <- function(x){
    ##x$vpa
    attr(x, "vpa")
}

valueLabels.cptable_class <- function(x){
    out <- list(attr(x, "levels"))
    nam <- attr(x, "vpa")
    names(out) <- attr(x, "vpa")[1] #x$vpa[1]
    out
}



