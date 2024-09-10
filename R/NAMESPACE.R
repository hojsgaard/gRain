
#' @useDynLib gRain

## Vanilla R imports (and exports)
## -------------------------------

#' @importFrom stats xtabs runif terms
#'     addmargins as.formula cov.wt fitted formula
#'     ftable getCall logLik loglin na.omit pchisq pf pnorm r2dtable
#'     terms update update.formula setNames
#' @importFrom utils combn str
#'
#' @importMethodsFrom stats4 plot


## Miscellaneous
## -------------
#' @importFrom stats simulate predict
#' @importFrom Rcpp evalCpp
#'
#' @import methods
#' @import gRbase
#' @importFrom broom tidy
#' 
#' @importFrom igraph 
#'     get.adjacency V "V<-" E "E<-" is.directed layout.lgl
#'     layout.graphopt plot.igraph graph.adjacency is.dag
#' 
### @importFrom bnlearn random.graph hc as.igraph as.grain
#' 
NULL

