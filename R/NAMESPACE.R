
#' @useDynLib gRain

## Vanilla R imports (and exports)
## -------------------------------

#' @importFrom stats xtabs runif terms
#'     addmargins as.formula cov.wt fitted formula
#'     ftable getCall logLik loglin na.omit pchisq pf pnorm r2dtable
#'     terms update update.formula 
#' @importFrom utils combn str
#'
#' @importMethodsFrom stats4 plot

#' 
## To make available in vignette 
#' @importFrom magrittr   "%>%"
#' @export "%>%" 

## Miscellaneous
## -------------
#' @importFrom stats simulate predict
#' @importFrom Rcpp evalCpp
#'
#' @import methods
#' @import gRbase
#'
#' @importFrom igraph 
#'     get.adjacency V "V<-" E "E<-" is.directed layout.lgl
#'     layout.graphopt plot.igraph graph.adjacency is.dag
#' 
#' @importFrom broom tidy
NULL

