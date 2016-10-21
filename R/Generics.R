#' @title gRain generics
#'
#' @description Generic functions etc for the gRain package
#'
#' @name grain-generics
#' 
#' @aliases nodeNames nodeStates nodeNames.grain
#'     nodeStates.grain
#'
#' @param x,object A relevant object.
#' @param nodes Some nodes of the object.
#' @param ... Additional arguments; currently not used.
#' 

#' @rdname grain-generics
nodeNames  <- function(x) UseMethod("nodeNames")

#' @rdname grain-generics
nodeNames.grain  <- function(x)
  x$universe$nodes

#' @rdname grain-generics
nodeStates <- function(x, nodes=nodeNames(x)) UseMethod("nodeStates")

#' @rdname grain-generics
nodeStates.grain <- function(x, nodes=nodeNames(x)){
  x$universe$levels[nodes]
}

#' @rdname grain-generics
universe <- function(object, ...) UseMethod("universe")

#' @rdname grain-generics
universe.grain <- function(object, ...) object$universe

#' @rdname grain-generics
varNames.grainEvidence_ <- function(x) x$summary$nodes

## #' @rdname grain-generics
## rip.grain <- function(object) object$rip
