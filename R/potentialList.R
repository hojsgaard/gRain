## Create potential list (rip, universe)
##
.mkArrayList <- function(rip.order, universe, values=1){
    cliques <- rip.order$cliques
    
    potlist  <- as.list(rep(NA, length(cliques)))
    
    for ( i in seq_along(cliques)){
        cq    <- cliques[[ i ]]
        vlab  <- universe$levels[cq]
        potlist[[ i ]] <- tab(cq, vlab, values)
    }
    potlist
}

.initArrayList <- function(x, values=NA){
    lapply(x, function(z) {
        z[] <- values             
        z
    } )
}

## Insert cpt's into potential list (cptlist, APlist)
##
.insertCPT <- function(cptlist, potlist, details=0)
{
    if (details>=1) cat(".Inserting cpt's in potential list [.insertCPT]\n")

    APnames <- lapply(potlist, function(x) names(dimnames(x)))
    CPnames <- unname(lapply(cptlist, function(x) varNames(x)))

    ## FIXME .findHost can be replaced by gRbase::get_superset_
    hosts    <- .findHosts( CPnames, APnames )

    for (i in 1:length(cptlist)) {
            cptc <- cptlist[[ i ]]
            j    <- hosts[ i ]
            ##print(potlist[[j]])
            potlist[[ j ]] <- tableOp( potlist[[ j ]], cptc, "*" )
        }
    .infoPrint(details, 4, {cat("....potlist (after insertion):\n"); print(potlist) })
    potlist
}




## FIXME .findHosts can be replaced by gRbase::get_superset_
.findHosts <- function( xx, yy ){
  unlist(lapply(1:length(xx), function( i ) which(isin(yy, xx[[  i  ]], index=T)>0)[1]))
}




## .insertNA <- function(list.of.tables)
## {
##     lapply(list.of.tables,
##            function(xxx)
##            {
##                xxx[] <- NA
##                xxx
##            } )
## }
##










