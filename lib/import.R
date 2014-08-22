#' wrapper for importing the table output from hmmsearch
import.pfamTable <- function(file) {
  # load pfam annotation file
  pfam <- read.table(file, sep = "", as.is = T)
  colnames(pfam) <- c("target", "accession", "query_name", "pfam",
                      "evalue", "score", "bias", "evalue2", "score2",
                      "bias2", "exp", "reg", "clu", "ov", "env", 
                      "dom", "rep", "inc", "description")
  return(pfam)
}

#'
import.pfam2go <- function(file) {
  data <- readLines(file)
  # remove header
  data <- data[which(!substring(data, 1, nchar('!')) == '!')]
  # convert to tab seperated strings
  data <- gsub(">","\t", data)
  data <- gsub(";","\t", data)
  data <- strsplit(data,"\t")
  # create data.frame
  data <- cbind(do.call(rbind.data.frame, 
                        lapply(data, function(x){
                          unlist(strsplit(substring(x[[1]][1], 6), "\ "))
                        })), 
                do.call(rbind.data.frame, data)[,2:3])
  names(data) <- c("id","name","go_name","go_id")
  return(data)
}

import.cog2go <- function(file) {
  data <- readLines(file)
  # remove header
  data <- data[which(!substring(data, 1, nchar('!')) == '!')]
  # convert to tab seperated strings
  data <- gsub(">","\t", data)
  data <- gsub(";","\t", data)
  data <- strsplit(data,"\t")
  # create data.frame
  data <- cbind(do.call(rbind.data.frame, 
                        lapply(data, function(x){
                          unlist(strsplit(substring(x[[1]][1], 4), "\ "))
                        })), 
                do.call(rbind.data.frame, data)[,2:3])
  names(data) <- c("id","name","go_name","go_id")
  return(data)
}

import.tigrfams2go <- function(file) {
  data <- readLines(file)
  # remove header
  data <- data[which(!substring(data, 1, nchar('!')) == '!')]
  # convert to tab seperated strings
  data <- gsub(">","\t", data)
  data <- gsub(";","\t", data)
  data <- strsplit(data,"\t")
  # create data.frame
  data <- cbind(do.call(rbind.data.frame, 
                        lapply(data, function(x){
                          unlist(strsplit(substring(x[[1]][1], 15), "\ "))
                        })), 
                do.call(rbind.data.frame, data)[,2:3])
  names(data) <- c("id","name","go_name","go_id")
  return(data)
}

read.obo<-function(obo)
{
  obo <- readLines(obo)
  ## FIND created date
  n <- grep("^date:", obo[1:10])
  date <- as.Date(substr(obo[n], 7,16), "%d:%m:%Y")
  ## empty list to store results
  obo1 <- vector("list")
  # counter for list elements
  z <- 0
  ## start of all stanzas
  n <- grep("^\\[", obo)
  if(length(n)==0){stop("No terms to parse")}
  # add last postion in obo file (+2 )
  n <- c(n, length(obo)+2 )
  for(i in 1:(length(n) - 1) ){
    ## ONLY Terms, skip Typedef
    if(obo[n[i]] == "[Term]"){
      # ids always first, name second?
      id <- gsub("(id: )(.*)", "\\2", obo[n[i] +1 ])
      name<- gsub("(name: )(.*)", "\\2", obo[n[i] +2 ])
      ## get all tag-value pairs for a single stanza
      x <- obo[(n[i] + 1):(n[i + 1] - 2)]
      ## Skip obolete terms
      if(any( grepl("^is_obsolete", x)) ){
        # print(paste(id, "is obsolete"))
        ## check parents - either missing for root OR starts with is_a: or relationship:
      }else{
        p <- x[grep("^is_a|^relationship", x)]
        ## check if missing = root
        if(length(p) == 0){
          print(paste("Note:", id, "is a ROOT node"))
          z <- z + 1
          obo1[[z]] <- c(id, name, NA, NA)
        }else{
          ## may have 1 or more parents...
          for(j in 1:length(p)){
            ## ID prefixes are usually upper-case (GO), some lower and upper (FBbt, LiPrO),
            ## and some with underscore (DC_CL)
            # may have commments after ID
            pid <- gsub("(.* )([[:alpha:]|\\_]+:[[:digit:]]+)(.*)", "\\2", p[j])
            ## capture relationship type??
            rship <- gsub("(.*)( [[:alpha:]|\\_]+:[[:digit:]]+)(.*)", "\\1", p[j])
            rship <- gsub("relationship: ", "", rship)
            rship <- gsub("is_a:", "is_a", rship)
            z <- z + 1
            obo1[[z]] <- c(id, name, rship, pid)
          }
        }
      }
    }
  }
  obo2 <- data.frame(do.call("rbind", obo1), stringsAsFactors=FALSE)
  colnames(obo2) <- c("id", "name", "related", "parent_id")
  class(obo2)<-c("obo", "data.frame")
  attr(obo2, "created") <- date
  attr(obo2, "downloaded") <- Sys.Date()
  print(paste("Loaded", length(unique(obo2$id)), "unique terms and", nrow(obo2), "rows"))
  obo2
}
