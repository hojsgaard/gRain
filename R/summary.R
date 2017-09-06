summary.grain <- function(object, type='std', ...){

    type <- match.arg(type, c("std", "cliques", "rip", "configurations"))
    
    cat("Independence network: Compiled:", .isComp(object), "Propagated:", .isProp(object), "\n")
    
    if (length(object$evidence)) getEvidence(object)
    
    cat(" Nodes :")
    utils::str(nodeNames(object)) ## $universe$nodes)
    
    if (.isComp(object)){
        
        cq.length <- sapply(rip(object)$cliques, length)
        
        cat(sprintf(" Number of cliques:              %4d \n",  length(cq.length)))
        cat(sprintf(" Maximal clique size:            %4d \n",  max(cq.length)))
        cat(sprintf(" Maximal state space in cliques: %4d \n",
                    max(unlist(lapply(pot(object)$pot_equi, length)) )))
      
        if(length(e <- getEvidence(object))){
            print(e)
        }
        
        switch(type,
               "rip"={
                   cat("\nRIP ordering:\n")
                   print(rip(object))
               },
               "cliques"={
                   cat("\nCliques:\n")
                   .printList(rip(object)$cliques)
               },
               "configurations"={
                   cat("\nConfigurations:\n")
                   nc <- seq_along(rip(object)$cliques)
                   for (i in nc){
                       cat(" clique ", i, ":", length(pot(object)$pot_equi[[i]]), "\n")
                   }
               })
    }
    invisible(object)
}

