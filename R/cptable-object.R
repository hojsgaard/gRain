#' @title Object oriented specification of Bayesian network
#' 
#' @description Object oriented specification of Bayesian
#'     network. This is experimental and may well change.
#'
#' @name object-cptable
#' 
#' @param object An object to be given a value.
#' @param bn A Bayesian network
#' @param value A value to be added to an object.
#'
#' @param x  An object
#' @param ... Additional arguments
#'
#' @details This is all experimental and may well change soon.
#'
#' @examples
#'
#' yn <- c("yes", "no")
#' ssp <- list(Rain = yn, Sprinkler = yn, GrassWet = yn)
#'
#' p.R <- parray("Rain", levels=ssp, values=c(.2, .8))
#' p.S_R <- parray(c("Sprinkler", "Rain"), levels = ssp,
#'    values=c(.01, .99, .4, .6))
#' p.G_SR <- parray(~ GrassWet:Sprinkler:Rain, levels = ssp,
#'    values=c(.99, .01, .8, .2, .9, .1, 0, 1))
#'
#' bn <- cpt_domain()
#' add_cpt(bn) <- p.R
#' add_cpt(bn) <- p.S_R
#' add_cpt(bn) <- p.G_SR
#' grain(compileCPT(bn))
#'
#' bn <- cpt_domain()
#' add_cpt(bn) <- cptable(~Rain, levels=ssp, values=c(.2, .8))
#' add_cpt(bn) <- cptable(~Sprinkler|Rain, levels=ssp,
#'    values=c(.01, .99, .4, .6))
#' add_cpt(bn) <- cptable(~GrassWet|Sprinkler:Rain, levels=ssp,
#'    values=c(.99, .01, .8, .2, .9, .1, 0, 1))
#'

#' @rdname object-cptable
cpt_domain <- function(){
    out <- list()
    class(out) <- "cpt_domain"
    out
}

#' @rdname object-cptable
print.cpt_domain <- function(x, ...){
    cat("cpt_domain\n")
    print.default(x)
    invisible(x)
}

#' @rdname object-cptable
add_cpt <- function(object, value){
    UseMethod("add_cpt")
}

#' @rdname object-cptable
"add_cpt.cpt_domain" <- function(object, value){
    out <- c(object, list(value))
    vn <- attr(value, "vpa")[1]
    names(out)[length(out)] <- vn
    class(out) <- "cpt_domain"
    out
}

#' @rdname object-cptable
"add_cpt<-" <- function(object, value)
    UseMethod("add_cpt<-")


#' @rdname object-cptable
"add_cpt<-.cpt_domain" <- function(object, value){
    add_cpt(object, value)
}


#' @rdname object-cptable
mvn_domain <- function(){
    out <- list()
    class(out) <- "mvn_domain"
    out
}

#' @rdname object-cptable
print.mvn_domain <- function(x, ...){
    cat("mvn_domain\n")
    print.default(x)
    invisible(x)
}

#' @rdname object-cptable
add_mvn <- function(object, value){
    UseMethod("add_mvn")
}

#' @rdname object-cptable
"add_mvn.mvn_domain" <- function(object, value){
    out <- c(object, list(value))
    vn <- attr(value, "vpa")[1]
    names(out)[length(out)] <- vn
    class(out) <- "mvn_domain"
    out
}

#' @rdname object-cptable
"add_mvn<-" <- function(object, value)
    UseMethod("add_mvn<-")


#' @rdname object-cptable
"add_mvn<-.mvn_domain" <- function(object, value){
    add_mvn(object, value)
}

#' @rdname object-cptable
compileMVN <- function(bn, ...){

    vpa_list <- lapply(bn, function(z) attr(z, "vpa"))
    vpa_list

    dg <- dagList(vpa_list, forceCheck=TRUE)

    vn <- unique(unlist(vpa_list))
    ss <- setdiff(vn, names(bn))

    if (length(ss) > 0)
        stop("No distribution given for node(s): ", toString(ss))
    
    if (length(names(bn)) != length(unique(names(bn))))
        stop("Nodes specified more than once: ", toString(names(bn)))

    out <- bn
    attr(out, "universe") <- list(nodes=vn)
    attr(out, "dag") <- dg
    class(out) <- "MVNspec"
    out
}







