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
#' @param details Debugging information
#' @return An object (a list) of class "huginNet".
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{grain}}
#' @references Søren Højsgaard (2012). Graphical Independence
#'     Networks with the gRain Package for R. Journal of Statistical
#'     Software, 46(10), 1-26.
#'     \url{http://www.jstatsoft.org/v46/i10/}.
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
#' @export loadNet

loadHuginNet <- function(file, description=rev(unlist(strsplit(file, "/")))[1],
                         details=0){

  xxx      <-.readHuginNet(file,details)
  yyy      <-.transformHuginNet2internal(xxx)
  universe <- .asUniverse(yyy)
  plist    <- lapply(yyy$potentialList, .hpot2cptable, universe)
  value    <- grain(compileCPT(plist))
  return(value)
}

loadNet <- function(file, description = rev(unlist(strsplit(file, "/")))[1], details = 0) {
    if(endsWith(file, '.xdsl')){
        plist <- .read_xdsl(file)
        grain(compileCPT(plist))
    } else if (endsWith(file, '.dne')){
        plist <- .read_dne(file)
        grain(compileCPT(plist))
    } else {
        loadHuginNet(file, description=rev(unlist(strsplit(file, "/")))[1],
                     details=0)
    }
}

.extract_dne <- function(lines, pattern, nodes=FALSE){
    presence <- grep(pattern, lines)
    v <- rep(NA, length(lines))
    nms <- sapply(strsplit(lines[presence], ifelse(nodes,' ','=')), '[[', 2)
    if(nodes){
        v <- c(rep(NA, presence[1]-1), nms[cumsum(seq_along(v) %in% presence)])
    } else {
        v[presence] <- gsub(';', '', trimws(gsub("[()]", "", nms)))
    }
    return(v)
}

# CHECK what happens if no decimal is provided (i.e. prob set to zero/one) 
.read_dne <- function(dne){
    f <- readLines(dne)
    nodes <- .extract_dne(f, "node (.*?) \\{", TRUE)
    sts <- .extract_dne(f, "states = (.*?);") # get all statuses
    w_sts <- which(!is.na(sts))
    prt <- .extract_dne(f, "parents = (.*?);")
    chc <- .extract_dne(f, "chance = (.*?);")
    knd <- .extract_dne(f, "kind = (.*?);")
    dsc <- .extract_dne(f, "discrete = (.*?);")
    tb <- data.frame(nodes=nodes[w_sts], 
                     status=sts[w_sts], 
                     parents=prt[!is.na(prt)][1:length(w_sts)],
                     chance=chc[!is.na(chc)][1:length(w_sts)],
                     kind=knd[!is.na(knd)][1:length(w_sts)],
                     discrete=dsc[!is.na(dsc)][1:length(w_sts)],
                     stringsAsFactors=FALSE)
    if(any(tb$chance != 'CHANCE')) {
        stop('Check node types of ', paste(tb$nodes[tb$chance != 'CHANCE'],collapse=','),
             ' from Netica file: only chance node types are allowed.')
    }
    if(any(tb$kind != 'NATURE')) {
        stop('Check node types of ',paste(tb$nodes[tb$kind != 'NATURE'],collapse=','),
             ' from Netica file, only nature node types are allowed.')
    }
    if(any(tb$discrete != 'TRUE')) {
        stop('Check node types of ',paste(tb$nodes[tb$discrete != 'TRUE'],collapse=','),
             ' from Netica file, only discrete nodes are allowed.')
        names(lst) <- tb$nodes
    }
    lst <- vector(mode = "list", length = nrow(tb))
    names(lst) <- nodes[w_sts]
    spb <- str_detect(f, "probs =") # get start of cpt
    stt <- str_detect(f, "title =") # get end of cpt
    ctrl <- FALSE
    for(i in 1:length(spb)){ # isolate between the two 
        if(stt[i]){ # Keep order of if!!
            ctrl <- FALSE
        }
        if(ctrl){
            lst[[nodes[i]]] <- c(lst[[nodes[i]]], trimws(gsub('\t\t','',f[i])))
        }
        if(spb[i]){
            ctrl <- TRUE
        }
    }
    states <- strsplit(tb$status, ', ')
    parents <- strsplit(tb$parents, ', ')
    # na.omit(as.numeric(unlist(strsplit(lst[[i]],"[^[:digit:].]"))))
    lst <- lapply(1:length(lst), function(i){
        # if(!is.null(lst[[i]])){ # case for node different from chance/nature
        a <- unlist(strsplit(lst[[i]], ','))
        probs <- as.numeric(regmatches(a,regexpr("(^|\\d+)\\.\\d+",a)))
        gRain::cptable(c(tb$nodes[i], parents[[i]]), states[[i]], matrix(probs))
        # }
    })
    names(lst) <- tb$nodes
    return(lst)
}

.read_xdsl <- function(xdsl){
    x <- xml2::read_xml(xdsl)    
    recs <- xml2::xml_find_all(x, "//cpt")
    nodes <- xml2::xml_attr(recs, "id")
    probs <- lapply(strsplit(xml2::xml_text(xml2::xml_find_all(x, "//probabilities")), ' '), as.numeric)
    lst <- lapply(1:length(recs), function(i){
        states <- xml2::xml_attr(xml2::xml_find_all(recs[i], './/state'), 'id')
        parents <- unlist(strsplit(xml2::xml_text(xml2::xml_find_all(recs[i], './/parents')), ' '))
        gRain::cptable(c(nodes[i], rev(parents)), states, probs[[i]])
    })
    names(lst) <- nodes
    return(lst)
}

.transformHuginNet2internal <- function(x){
  nodeList2 <- lapply(x$nodeList, .getNodeSpec)
  potentialList2 <- lapply(x$potentialList, .getPotentialSpec)

  nl <- .makeNodeNamesUnique(nodeList2)

  repeat{
    if (length(nl$nonunique)==0)
      break()
    nl <- .makeNodeNamesUnique(nl$nodeList)
  }

  nodeList2 <- nl$nodeList

  value <- structure(list(nodeList=nodeList2, potentialList=potentialList2))
  class(value)<- "huginnet"

  return(value)
}


.readHuginNet <- function(file, details=0){

  .infoPrint(details, 1, cat(".HUGIN netfile:", file,"\n"))
  nodeCount <- 0
  con <- file(file, "rb")
  repeat{
    cline <- .getLine(con);  #print(cline)
    if (!length(cline))
      break()

    if (.hasToken("node", cline)) ## Fragile if 'node' is the name of a variable...
      nodeCount <- nodeCount + 1
  }
  close(con)

  .infoPrint(details, 3, cat("...there are around", nodeCount, "nodes \n"))

  ## Data structure for holding specification (possibly too long)
  ##
  nodeList <- potentialList <- as.list(rep(NA, nodeCount))

  con <- file(file, "rb")
  currNode <- currPotential <- 1
  state<-"start"
  repeat{
    cline <- .getLine(con);  #print(cline)
    if (!length(cline))
      break()
    switch(state,
           "start"={
             if (.hasToken("net",cline)){
               state="net"
               .infoPrint(details, 2, cat("..NET action\n"))
               wline <- cline
             }
           },
           "net"={
             wline <- c(wline, cline)
             if (.hasToken("}",cline)){
               state="run1"
               .infoPrint(details,2,cat("..end NET action\n"))
             }
           },
           "run1"={
             if (.hasToken("node", cline)){
               state="node"
               .infoPrint(details, 2, cat("..NODE action\n"))
             } else {
               if (.hasToken("potential", cline)){
                 state="potential";
                 .infoPrint(details,2, cat("..POTENTIAL action\n"))
               }
             }
             wline <- cline
           },
           "node"={
             wline <- c(wline, cline)
             if (.hasToken("}",cline)){
               state="run1";
               .infoPrint(details,2,cat("..end NODE action\n"))
               nodeList[[currNode]] <- wline;
               currNode <- currNode + 1
             }
           },
           "potential"={
             wline <- c(wline, cline)
             if (.hasToken("}",cline)){
               state="run1";
               .infoPrint(details,2, cat("..end POTENTIAL action\n"))
               potentialList[[currPotential]] <- wline;
               currPotential <- currPotential + 1
             }
           }
           )
  }
  close(con)

  nodeList <- nodeList[!sapply(lapply(nodeList, is.na),all)]
  potentialList <- potentialList[!sapply(lapply(potentialList, is.na),all)]


  value <- structure(list(nodeList=nodeList, potentialList=potentialList))
  return(value)
}



.asUniverse <- function(from){
  ccshort   <-sapply(from$nodeList, function(x)x$nodeVar)
  ccnames   <-sapply(from$nodeList, function(x)x$nodeLabel)
  cclabels  <-lapply(from$nodeList, function(x)x$nodeStates)
  names(cclabels) <- ccnames
  di <- c(lapply(cclabels, length),recursive=TRUE)
  list(nodes=ccnames, short=ccshort, levels=cclabels, nlev=di)
}


.hpot2cptable <- function(cpot, universe){
  idx <- match(c(cpot[c("nodeVar","parentVar")],recursive=TRUE), universe$short)
  vpa <- universe$nodes[idx]
  v   <- vpa[1]
  cptable(vpa, values=cpot$potential, levels=universe$levels[[v]])
}




.getLine   <- function(con) {
  readLines(con, n=1)
}

.hasToken  <- function(token, cline) {
  ##print(cline)
  cline <- gsub("^ +","",cline)
  a <- unlist(strsplit(cline," "))[1]

  if (!is.na(a))
    a==token
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

## .toCamel <- function(s){
##   s<-gsub(" +"," ",s)
##   s<-unlist(strsplit(s, " "))
##   paste(sapply(s, .capWords),collapse='')
## }

.toCamel <- function(s){
  s<-gsub(" +"," ",s)
  s<-unlist(strsplit(s, " "))
  paste(c(s[1],sapply(s[-1], .capWords)),collapse='')
}



.getNodeSpec <- function(nodeSpec){

  tmp <- nodeSpec[.tokenIdx("node", nodeSpec)]
  nodeVar <- gsub("node +","",tmp)[1]
  nodeVar <- gsub(" +","",nodeVar)


  tmp <- nodeSpec[.tokenIdx("label", nodeSpec)]
  nodeLabel <- gsub(" +label += +","",tmp);
  nodeLabel <- gsub(";", "", nodeLabel)
  nodeLabel <- gsub('"',"", nodeLabel)

  nodeLabel <- gsub(" +"," ",nodeLabel)

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
  nodeStates <- gsub(" +states += +","",tmp);
  nodeStates <- gsub("[\\(,\\);]","",nodeStates);
  nodeStates <- unlist(strsplit(nodeStates, '\\"'))
  nodeStates <- sapply(nodeStates, function(d) gsub("^ +","",d))
  nodeStates <- nodeStates[sapply(nodeStates, nchar)>0]

  nodeStates <- sapply(nodeStates, .toCamel)

  nodeStates <- gsub(" +",".", nodeStates)
  names(nodeStates)<-NULL


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
  value
}




.makeNodeNamesUnique <- function(nodeList2){
  nl<-t(sapply(nodeList2, function(d)unlist(d[1:2])))

  nonunique <- names(which(table(nl[,2])>1))

  if (length(nonunique)){
    cat ("Label(s): {", nonunique, "} appears mode than once in NET file\n")
    for (i in 1:length(nonunique)){
      cnu <- nonunique[i]
      idx<-which(cnu ==nl[,2])
      for (j in idx){
        a <- nodeList2[[j]]$nodeVar
        cat("  Replacing label", cnu, " with node name", a, "\n")
        nodeList2[[j]]$nodeLabel <- a
      }
    }
  }

  return(list(nodeList=nodeList2, nonunique=nonunique))
}


#' @rdname load-save-hugin
saveHuginNet <- function(gin, file, details=0){

    if (!inherits( gin, "grain"))
        stop("Not a grain object")
    
    if (is.null(gmd <- getgin(gin, "universe")))
        stop("Strange error: no universe in network")

    if (is.null(cptlist <- getgin(gin, "cptlist"))){
        cat("Object does not have 'cptlist' component; creating one for you...\n")
        cptlist <- mkcptlist(gin)
    }
        
    vlab <- gmd$levels
    vnam <- gmd$nodes
    nn   <- length(vlab)
    
    th     <- cumsum(c(0,rep(2*pi/nn, nn-1)))
    r      <- 100
    coords <- lapply(th, function(d) round(r+r*c(cos(d), sin(d))))
    
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









































