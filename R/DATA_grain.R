#' @title Chest clinic example
#' @description Conditional probability tables for the chest clinic example.
#' @name example_chest
#' @aliases example_chest_cpt
#' @docType data
#' @keywords datasets
#' @examples
#'
#' yn   <- c("yes", "no")
#' a    <- cpt(~asia, values=c(1,99),levels=yn)
#' t.a  <- cpt(~tub|asia, values=c(5,95,1,99),levels=yn)
#' s    <- cpt(~smoke, values=c(5,5), levels=yn)
#' l.s  <- cpt(~lung|smoke, values=c(1,9,1,99), levels=yn)
#' b.s  <- cpt(~bronc|smoke, values=c(6,4,3,7), levels=yn)
#' e.lt <- cpt(~either|lung:tub,values=c(1,0,1,0,1,0,0,1),levels=yn)
#' x.e  <- cpt(~xray|either, values=c(98,2,5,95), levels=yn)
#' d.be <- cpt(~dysp|bronc:either, values=c(9,1,7,3,8,2,1,9), levels=yn)
#'
#' chest_cpt <- list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be)
#' ## bn <- grain(compile_cpt(chest_cpt))
NULL

#' @title Wet grass example
#' @description Conditional probability tables for the wet grass example.
#' @name example_grass
#' @aliases example_grass_cpt
#' @docType data
#' @keywords datasets
#' @examples 
#' yn <- c("yes", "no")
#' p.R    <- cpt(~R, values=c(.2, .8), levels=yn)
#' p.S_R  <- cpt(~S:R, values=c(.01, .99, .4, .6), levels=yn)
#' p.G_SR <- cpt(~G:S:R, values=c(.99, .01, .8, .2, .9, .1, 0, 1), levels=yn)
#'
#' grass_cpt <- list(p.R, p.S_R, p.G_SR)
#' ## bn <- grain(compile_cpt(grass_cpt))
NULL
