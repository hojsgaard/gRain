#' @method plot grain
#' @export
plot.grain <- function(x, type, ...) {
    ##.primiplot(x$dag)
    .primiplot <- function(ig) {
        V(ig)$label <- V(ig)$name
        V(ig)$size  <- 40
        ig$cex   <-  4
        ig$layout   <- layout.graphopt
        plot(ig)
    }
    
    if (missing(type)) {
        if (isCompiled(x)) {
            .primiplot(x$ug)
        } else {
            if ("pot_grain" %in% class(x)) {
                .primiplot(x$ug)
            } else {
                .primiplot(x$dag)
            }
        }
    } else {
        zz <- x[[type]]
        if (!is.null(zz))
            .primiplot(zz)
        else
            cat("Slot", type, "does not exist \n")
    }
}





