##
## Reading / writing Bayesian networks from / to HUGIN net files
##

#' @title Load and save Hugin net files
#' 
#' @description These functions can load a net file saved in the
#'     'Hugin format' into R and save a network in R as a file in the
#'     'Hugin format'.
#'
#' @name load-save-hugin
#' 
#' @aliases loadHuginNet saveHuginNet
#' @param gin An independence network
#' @param file Name of HUGIN net file. Convenient to give the file the
#'     extension '.net'
#' @param description A text describing the network, defaults to
#'     \code{file}
#' @param details Debugging information.
#'
#' @note
#'
#' * In Hugin, it is possible to specify de potential of a node as a
#'   functional relation between other nodes.  In a .net file, such a
#'   specification will appear as 'function' rather than as
#'   'node'. Such a specification is not recognized by `loadHuginNet`.
#'
#' * It is recommended to avoid the text `node` as part of the name of
#'   a node.
#' 
#' @return An object of class `grain`.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{grain}}
#' @references Søren Højsgaard (2012). Graphical Independence
#'     Networks with the gRain Package for R. Journal of Statistical
#'     Software, 46(10), 1-26.
#'     \url{https://www.jstatsoft.org/v46/i10/}.
#' @keywords utilities
#' @examples
#' 
#' 
#' ## Load HUGIN net file
#' tf <- system.file("huginex", "chest_clinic.net", package = "gRain")
#' chest <- loadHuginNet(tf, details=1)
#' chest 
#' 
#' ## Save a copy
#' td <- tempdir()
#' saveHuginNet(chest, paste(td,"/chest.net",sep=''))
#' 
#' ## Load the copy
#' chest2 <- loadHuginNet(paste(td,"/chest.net",sep=''))
#' 
#' tf <- system.file("huginex", "golf.net", package = "gRain")
#' golf <- loadHuginNet(tf, details=1)
#' 
#' saveHuginNet(golf, paste(td,"/golf.net",sep=''))
#' golf2 <- loadHuginNet(paste(td,"/golf.net",sep=''))
#' 
#' @export loadHuginNet
loadHuginNet <- function(file, description=NULL, details=0){

    ## if (is.null(description))
        ## description <- rev(unlist(strsplit(file, "/")))[1]
    ## xxx      <- read_hugin_net(file, details)
    ## yyy      <- hugin_net2internal(xxx)
    ## universe <- as_universe(yyy)
    ## plist    <- lapply(yyy$pot_list, hugin_pot2cptable, universe)

    plist <- loadHuginNet2(file, description=description, details=details)    
    value    <- grain(compileCPT(plist))    
    return(value)
}

#' @export loadHuginNet
loadHuginNet2 <- function(file, description=NULL, details=0){

    if (is.null(description))
        description <- rev(unlist(strsplit(file, "/")))[1]
    xxx      <- read_hugin_net(file, details)
    yyy      <- hugin_net2internal(xxx)
    universe <- as_universe(yyy)
    plist    <- lapply(yyy$pot_list, hugin_pot2cptable, universe)
    plist
}

read_hugin_net <- function(file, details=0){

    .infoPrint(details, 1, cat(".HUGIN netfile:", file,"\n"))
    con <- file(file, "rb")
    node_count <- 0
    
    repeat{
        curr_line <- get_line(con);  #print(curr_line)
        if (!length(curr_line))
            break()
        
        if (.hasToken("node", curr_line)) ## Fragile if 'node' is the name of a variable...
            node_count <- node_count + 1
    }
    close(con)
    
    .infoPrint(details, 3, cat("...there are around", node_count, "nodes \n"))
    
    ## Data structure for holding specification (possibly too long)
    ##
    node_list <- pot_list <- as.list(rep(NA, node_count))
    
    con <- file(file, "rb")
    curr_node <- curr_pot <- 1
    state <- "start"
    repeat{
        curr_line <- get_line(con);  #print(curr_line)
        if (!length(curr_line))
            break()
        switch(state,
               "start"={
                   if (.hasToken("net", curr_line)){
                       state="net"
                       .infoPrint(details, 2, cat("..NET action\n"))
                       wline <- curr_line
                   }
               },
               "net"={
                   wline <- c(wline, curr_line)
                   if (.hasToken("}", curr_line)){
                       state="run1"
                       .infoPrint(details,2,cat("..end NET action\n"))
                   }
               },
               "run1"={
                   if (.hasToken("node", curr_line)){
                       state="node"
                       .infoPrint(details, 2, cat("..NODE action\n"))
                   } else {
                       if (.hasToken("potential", curr_line)){
                           state="potential";
                           .infoPrint(details,2, cat("..POTENTIAL action\n"))
                       } else {
                           if (.hasToken("function", curr_line)){
                               close(con)
                               stop("specification of distribution via 'function' is not permitted\n")
                           }
                           
                       }
                   }
                   wline <- curr_line
               },
               "node"={
                   wline <- c(wline, curr_line)
                   if (.hasToken("}", curr_line)){
                       state="run1";
                       .infoPrint(details, 2, cat("..end NODE action\n"))
                       node_list[[curr_node]] <- wline;
                       curr_node <- curr_node + 1
                   }
               },
               "potential"={
                   wline <- c(wline, curr_line)
                   if (.hasToken("}",curr_line)){
                       state="run1";
                       .infoPrint(details,2, cat("..end POTENTIAL action\n"))
                       pot_list[[curr_pot]] <- wline;
                       curr_pot <- curr_pot + 1
                   }
               }
               )
    }
    close(con)
    
    node_list <- node_list[!sapply(lapply(node_list, is.na), all)]
    pot_list  <- pot_list [!sapply(lapply(pot_list,  is.na), all)]
    
    value <- structure(list(node_list=node_list, pot_list=pot_list))
    return(value)
}




hugin_net2internal <- function(x){

    node_list2 <- lapply(x$node_list, .getNodeSpec) 
    pot_list2 <- lapply(x$pot_list, .getPotentialSpec)
    
    nl <- .makeNodeNamesUnique(node_list2)
    
    repeat{
        if (length(nl$nonunique)==0)
            break()
        nl <- .makeNodeNamesUnique(nl$node_list)
    }
    
    node_list2 <- nl$node_list
    
    value <- structure(list(node_list=node_list2, pot_list=pot_list2))
    class(value)<- "huginnet"
    
    return(value)
}




as_universe <- function(from){
  ccshort   <-sapply(from$node_list, function(x)x$nodeVar)
  ccnames   <-sapply(from$node_list, function(x)x$nodeLabel)
  cclabels  <-lapply(from$node_list, function(x)x$nodeStates)
  names(cclabels) <- ccnames
  di <- c(lapply(cclabels, length),recursive=TRUE)
  list(nodes=ccnames, short=ccshort, levels=cclabels, nlev=di)
}


hugin_pot2cptable <- function(cpot, universe){
  idx <- match(c(cpot[c("nodeVar","parentVar")],recursive=TRUE), universe$short)
  vpa <- universe$nodes[idx]
  v   <- vpa[1]
  cptable(vpa, values=cpot$potential, levels=universe$levels[[v]])
}

get_line   <- function(con) {
  readLines(con, n=1)
}

.hasToken  <- function(token, curr_line) {
  ##print(curr_line)
  curr_line <- gsub("^ +", "", curr_line)
  a <- unlist(strsplit(curr_line, " "))[1]

  if (!is.na(a))
    a == token
  else
    FALSE
}

.tokenIdx <- function(token, x){
  idx <- which(as.logical(lapply(x, function(d) grep(token,d))))
  idx
}


.capWords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s,1,1)),
                           {s <- substring(s,2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

.toCamel <- function(s){
  s <- gsub(" +", " ", s)
  s <- unlist(strsplit(s, " "))
  paste(c(s[1], sapply(s[-1], .capWords)), collapse='')
}


.getNodeSpec <- function(nodeSpec){

    tmp <- nodeSpec[.tokenIdx("node", nodeSpec)]
    nodeVar <- gsub("node +", "", tmp)[1]
    nodeVar <- gsub(" +", "", nodeVar)
    
    tmp <- nodeSpec[.tokenIdx("label", nodeSpec)]
    nodeLabel <- gsub(" +label += +", "", tmp);
    nodeLabel <- gsub(";" , "" , nodeLabel)
    nodeLabel <- gsub('"' , "" , nodeLabel)
    nodeLabel <- gsub(" +", " ", nodeLabel)
    nodeLabel <- gsub(" ", "", nodeLabel)  ## FIXME IS THIS OK?
    
    if (length(nodeLabel) && nchar(nodeLabel)>0){
        
        nodeLabel <- .toCamel(nodeLabel)
        
        nl <- gsub("[^[:alnum:]]","",nodeLabel)
        nodeLabel <- gsub("[^[:alnum:]|\\.]","",nodeLabel)
        
        base<-as.character(0:9)
        if(subsetof(unlist(strsplit(nl,"")), base)){
            nodeLabel <- paste("X",nodeLabel,sep='')
        }
    } else {
        ##if (nchar(nodeLabel)==0)
        nodeLabel <- nodeVar
    }

    tmp <- nodeSpec[.tokenIdx("states", nodeSpec)]
    nodeStates <- gsub(" +states += +", "", tmp);
    nodeStates <- gsub("[\\(,\\);]", "", nodeStates);
    nodeStates <- unlist(strsplit(nodeStates, '\\"'))
    nodeStates <- sapply(nodeStates, function(d) gsub("^ +", "", d))
    nodeStates <- nodeStates[sapply(nodeStates, nchar) > 0]

    nodeStates <- sapply(nodeStates, .toCamel)
    
    nodeStates <- gsub(" +", ".", nodeStates)
    names(nodeStates) <- NULL
    
    value <- list(nodeVar=nodeVar, nodeLabel=nodeLabel, nodeStates=nodeStates)
    value
}

.getPotentialSpec <- function(potSpec){
  tmp <- potSpec[.tokenIdx("potential", potSpec)]
  tmp <- gsub("potential +","", tmp)
  tmp <- gsub("[\\(,\\),|]","", tmp)
  tmp <- gsub(" +"," ", tmp)
  tmp <- unlist(strsplit(tmp," "))
  tmp <- tmp[sapply(tmp, nchar)>0]

  nodeVar <- tmp[1]
  parentVar <- tmp[-1]

  sss  <- paste(potSpec,collapse="") ##; ss <<- sss
  sss2 <- gsub("^.*data[[:space:]]*=([^;]*);(.*)", "\\1", sss) ##; ss2<<-sss2

  ##sss3: ((( 0.5 1.2E-5 ) ( 3E3 0.5 )) ( 0.5 0.5 ) ( 0.5 0.5 )))
  sss3 <- gsub("\\)[^\\)]*\\(", ") (", sss2) ##; ss3<<-sss3

  ## sss4: "  0.5 1.2E-5   3E3 0.5   0.5 0.5   0.5 0.5 "s
  sss4 <- gsub("[\\(,\\),\\}]","", sss3)

  ## sss5: remove leading white space: "0.5 1.2E-5   3E3 0.5   0.5 0.5   0.5 0.5 "
  sss5 <- gsub("^[[:space:]]*","",sss4)
  ## sss6: remove trailing white space: "0.5 1.2E-5   3E3 0.5   0.5 0.5   0.5 0.5"
  sss6 <- gsub("[[:space:]]$*","",sss5)
  ## sss7: split to atoms
  sss7 <- strsplit(sss6, " +")[[1]]

  ###: Now create numerical values
  pot <- as.numeric( sss7 )

  value <- list(nodeVar=nodeVar, parentVar=rev(parentVar), potential=pot)
  ## str(value)
  value
}




.makeNodeNamesUnique <- function(node_list2){
  nl<-t(sapply(node_list2, function(d)unlist(d[1:2])))

  nonunique <- names(which(table(nl[,2])>1))

  if (length(nonunique) > 0){
    cat ("Label(s): {", nonunique, "} appears mode than once in NET file\n")
    for (i in 1:length(nonunique)){
      cnu <- nonunique[i]
      idx <- which(cnu == nl[, 2])
      for (j in idx){
        a <- node_list2[[j]]$nodeVar
        cat("  Replacing label", cnu, " with node name", a, "\n")
        node_list2[[j]]$nodeLabel <- a
      }
    }
  }

  return(list(node_list=node_list2, nonunique=nonunique))
}



#' @export
#' @rdname load-save-hugin
saveHuginNet <- function(gin, file, details=0){

    if (!inherits(gin, "grain"))
        stop("Not a grain object")
    
    if (is.null(gmd <- getgin(gin, "universe")))
        stop("Strange error: no universe in network")

    if (is.null(cptlist <- getgin(gin, "cptlist"))){
        cat("Object does not have 'cptlist' component; creating one for you...\n")
        cptlist <- make_cptlist(gin)
    }
        
    vlab <- gmd$levels
    vnam <- gmd$nodes
    nn   <- length(vlab)
    
    th     <- cumsum(c(0, rep(2 * pi / nn, nn - 1)))
    r      <- 100
    coords <- lapply(th, function(d) round(r + r * c(cos(d), sin(d))))
    
    con <- file(file, "wb")

    ## Write (trivial) net specification
    ##
    writeLines("net\n{", con)
    writeLines("  node_size = (100 30);", con)
    writeLines("\n}\n\n", con)
    
    ## Write node specification
    ##
    for (ii in 1:length(vlab)){
        st <-paste("node ", vnam[ii],"\n","{","\n",sep='')
        writeLines(st, con, sep="")
        ## cat(st)
        st <- paste("   label = \"\";","\n")
        writeLines(st, con, sep="")
        ## cat(st)
        st <- paste("   position = (", paste(coords[[ii]], collapse=' '), ");\n")
        writeLines(st, con, sep="")
        ## cat(st)
        
        st2 <- sapply(vlab[[ii]], function(d) paste('"',d,'"',sep=''))
        st  <- paste("   states = (", paste(st2, collapse=' '), ");\n")
        writeLines(st, con, sep="")
        ## cat(st)
        st <- paste("}\n")
        writeLines(st, con, sep="")
        ## cat(st)
    }
    

    for (ii in 1:length(cptlist)){

        cpot <- cptlist[[ii]]
        nam <- varNames(cpot)    ## BRIS
        lev <- valueLabels(cpot) ## BRIS
        val <- cpot              ## BRIS
        
        v  <- nam[1]
        pa <- nam[-1]
        
        lev   <- rev(lev[-1])
        wval  <- val
        if (length(lev)>0){
            for (kk in 1:length(lev)){
                ##print("splitVec:"); print(wval); print(class(wval))
                wval<-splitVec(wval,length(lev[[kk]]))
            }
        }
        ##print(wval); print(class(wval))
        plx <- printlist(wval)
        
        if (length(pa)){
            st <- paste("potential (",v, "|", paste(rev(pa), collapse=' '),")\n")
            writeLines(st,con,sep="")
            ## cat(st)
            st <- "{\n";
            writeLines(st,con,sep="")
            ## cat(st)
            st <- paste("   data = \n")
            writeLines(st,con,sep="")
            ## cat(st)
            ##a<-lapply(plx, cat, "\n")
            a<-lapply(plx, writeLines, con, sep="\n")
            st <- paste(";\n}\n")
            writeLines(st,con,sep="")
            ## cat(st)
            
        } else {
            st <- paste("potential (", v, ")\n")
            writeLines(st,con,sep="")
            ## cat(st)
            st <- "{\n";
            writeLines(st,con,sep="")
            ## cat(st)
            st <- paste("   data = ", plx, ";\n")
            writeLines(st,con,sep="")
            ## cat(st)
            st <- "}\n\n";
            writeLines(st,con,sep="")
            ## cat(st)
        }
    }
    
  close(con)
}









































