####
#### Re-implementation of simulate() function - quite fast...
#### Bristol, March 2008
####

## ##################################################################
##
#' @title Simulate from an independence network
#' @description Simulate data from an independence network.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @name grain-simulate
##
## ##################################################################
#'
#' @param object An independence network
#' @param nsim Number of cases to simulate
#' @param seed An optional integer controlling the random number
#'     generatation
#' @param \dots Not used...
#' @return A data frame

#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords models
#' @examples
#' 
#' tf <- system.file("huginex", "chest_clinic.net", package = "gRain")
#' 
#' chest <- loadHuginNet(tf, details=1)
#' simulate(chest,n=10) 
#' 
#' chest2 <- setFinding(chest, c("VisitToAsia", "Dyspnoea"),
#'                             c("yes", "yes"))
#' simulate(chest2, n=10)
#' 
#' @export simulate.grain
simulate.grain <- function(object, nsim=1, seed=NULL, ...){

    if (!.isComp(object)){
        ##cat("Compiling (and propagating) model ...\n")
        object <- compile(object, propagate=TRUE)
    } else {
        if (!.isProp(object)){
            ## cat("Propagating model...\n")
            object <- propagate(object)
        }
    }
    
    plist  <- pot(object)$pot_equi
    cqlist <- rip(object)$cliques
    splist <- rip(object)$separators
    
    ## Init
    ans           <- matrix(0, nrow=nsim, ncol=length(nodeNames(object)))
    colnames(ans) <- nodeNames(object)
    
    ctab  <- plist[[1]]
    res   <- simulateArray(x=ctab, nsim=nsim)
    ans[,colnames(res)] <- res
    
    ## Iterate
    if (length(cqlist) > 1){
        for (ii in 2:length(cqlist)){
            ctab <- plist[[ii]]
            vn   <- names(dimnames(ctab))
            sp   <- splist[[ii]] ## What we condition on
            if (length(sp) > 0){
                mtab <- tableMargin(ctab, sp)     ## FIXME: Old table-function
                ctab <- tableOp2(ctab, mtab, `/`) ## FIXME: Old table-function
            }
            rr   <- setdiff(vn, sp) ## Variables to be simulated
            ##cat("r:", rr, "s:", sp, "\n")
            if (length(sp)){
                spidx <- match(sp, vn)
                res   <- matrix(0, nrow=nsim, ncol=length(rr))
                colnames(res) <- rr
                un    <- ans[, sp, drop=FALSE]
                ##cat("un:\n"); print(un)
                vals  <- unique(un)
                sc    <- cumprod(apply(vals, 2, max) )
                sc    <- c(1, sc)[1:length(sc)]
                key   <- ((un - 1) %*% sc) + 1
                ##cat(sprintf("key=%s\n", toString(key)))  #browser()
        for(kk in unique(key)){
            nn   <- sum(kk == key)
            idx  <- un[match(kk, key),]
            res[kk==key,] <- simulateArray(ctab, nsim=nn, margin=spidx, value.margin=idx)
        }
            } else {
                res <- simulateArray(x=ctab, nsim=nsim)
            }
            ans[, colnames(res)] <- res
        }
    }
    
    ns <- nodeStates(object)
    vn <- colnames(ans)
    out <- vector("list", ncol(ans))
    names(out) <- vn
    for (jj in 1:ncol(ans)){
        out[[jj]] <- factor(ans[,jj], levels=seq(ns[[jj]]))
        levels(out[[jj]]) <- ns[[jj]]
    }
    out <- as.data.frame(out)
    names(out) <- vn
    out
}



##   ans <- as.data.frame(ans)
##   vn <- names(ans)

##   for (jj in 1:ncol(ans)){
##     #match(vn[jj], names(ns))
##     ans[,jj] <- factor(ans[,jj], levels=seq(ns[[jj]]))
##     levels(ans[,jj]) <- ns[[jj]]
##   }

  #return(ans)


      ##cat(sprintf("vn=%s sp=%s\n", toString(vn), toString(sp)))
      ##cat("ctab:\n");  print(ctab)
      ##cat("mtab:\n"); print(mtab)
      ##cat("ctab (updated):\n"); print(ctab)
