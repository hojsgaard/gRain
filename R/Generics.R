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

#' @rdname grain-generics
rip.grain <- function(object, ...)
    getgin(object, "rip")


#' @rdname grain-generics
uni <- function(object)
    UseMethod("uni")

#' @rdname grain-generics
uni.grain <- function(object)
    getgin(object, "universe")

#' @rdname grain-generics
pot <- function(object)
    UseMethod("pot")

#' @rdname grain-generics
pot.grain <- function(object)
    getgin(object, "potential")

#' @rdname grain-generics
cpt <- function(object)
    UseMethod("cpt")

#' @rdname grain-generics
cpt.cpt_grain <- function(object)
    getgin(object, "cptlist")

#' @rdname grain-generics
#' @param position Where to insert 'value'
#' @param value Value to insert at 'position'
"cpt<-" <- function(object, position, value){
    UseMethod("cpt<-")
}

#' @rdname grain-generics
"cpt<-.cpt_grain" <- function(object, position, value){
    object$cptlist[position] <- value
}



#' @rdname grain-generics
potential <- function(object)
   UseMethod("potential")

#' @rdname grain-generics
potential.grain <- function(object)
   object$potential

#' @rdname grain-generics
vpar.cpt_spec <- function(object, ...){
    lapply(object, function(u) names(dimnames(u)))
}

#' @rdname grain-generics
vpar.cpt_grain <- function(object, ...){
    lapply(getgin(object, "cptlist"), function(u) names(dimnames(u)))
}



.isComp <- function(x) getgin(x, "isCompiled")

.isProp <- function(x) getgin(x, "isPropagated")

