#' @title gRain generics
#'
#' @description Generic functions etc for the gRain package
#'
#' @name generics
#' 
#' @param object A relevant object.
#' @param nodes Some nodes of the object.
## #' @param value Value to be set for slot in object.
#' @param ... Additional arguments; currently not used.
#' 

#' @export
#' @rdname generics
nodeNames  <- function(object)
{
    UseMethod("nodeNames")
}

#' @export
#' @rdname generics
nodeNames.grain  <- function(object)
{
    getgrain(object, "universe")$nodes
}

#' @export
#' @rdname generics
nodeStates <- function(object, nodes=nodeNames(object))
{
    UseMethod("nodeStates")
}

#' @export
#' @rdname generics
nodeStates.grain <- function(object, nodes=nodeNames(object))
{
  getgrain(object, "universe")$levels[nodes]
}

#' @export
#' @rdname generics
universe <- function(object, ...)
{
    UseMethod("universe")
}

#' @export
#' @rdname generics
universe.grain <- function(object, ...)
{
    getgrain(object, "universe")
}

#' @rdname generics
#' @export
isCompiled <- function(object) {
    getgin(object, "isCompiled")
}

#' @rdname generics
#' @export
isPropagated <- function(object) {
    getgin(object, "isPropagated")
}

## ' @rdname generics
## ' @export
"isCompiled<-" <- function(object, value)
{
    object$isCompiled <- value
    object
}

## ' @rdname generics
## ' @export
"isPropagated<-" <- function(object, value)
{
    object$isPropagated <- value
    object
}

## ---------------------------------------------------------------
##
## Methods where generic function is in gRbase.
##
## ---------------------------------------------------------------


#' @export
#' @rdname generics
vpar.cpt_spec <- function(object, ...)
{
    lapply(object, function(u) names(dimnames(u)))
}

#' @export
#' @rdname generics 
vpar.cpt_grain <- function(object, ...)
{
    lapply(getgin(object, "cptlist"), function(u) names(dimnames(u)))
}

#' @export
#' @rdname generics
rip.grain <- function(object, ...)
{
    getgin(object, "rip")
}

## #' @export
## #' @rdname generics

## varNames.grainEvidence_ <- function(x)
## {
##     getgrain(x, "summary")$nodes
## }




















## #' @rdname generics
## uni <- function(object)
##     UseMethod("uni")

## #' @rdname generics
## uni.grain <- function(object)
##     getgin(object, "universe")

## #' @rdname generics
## pot <- function(object)
##     UseMethod("pot")

## #' @rdname generics
## pot.grain <- function(object)
##     getgin(object, "potential")


## #' @rdname generics
## cpt <- function(object)
##     UseMethod("cpt")

## #' @rdname generics
## cpt.cpt_grain <- function(object)
##     getgin(object, "cptlist")

## #' @rdname generics
## #' @param position Where to insert 'value'
## #' @param value Value to insert at 'position'
## "cpt<-" <- function(object, position, value){
##     UseMethod("cpt<-")
## }

## #' @rdname generics
## "cpt<-.cpt_grain" <- function(object, position, value){
##     object$cptlist[position] <- value
## }


## #' @rdname generics
## potential <- function(object)
##    UseMethod("potential")

## #' @rdname generics
## potential.grain <- function(object)
##    object$potential

