#' read.pfam
#'
#' wrapper for read a pfam file in table format into R
#' 
#'@param file           pfam table input file
#'
#'@return data.frame
#'@export
#'
read.pfam <- function(file) {
  # load pfam annotation file
  pfam <- read.table(file, sep = "", as.is = T)
  # change colnames of data.frame
  colnames(pfam) <- c("target", "accession", "query_name", "pfam",
                      "evalue", "score", "bias", "evalue2", "score2",
                      "bias2", "exp", "reg", "clu", "ov", "env", 
                      "dom", "rep", "inc", "description")
  return(pfam)
}
#' read.pfam2go
#'
#' wrapper for read pfam2go mapping files provided by GO Consortium
#' 
#'@param file           mapping input file
#'
#'@return data.frame
#'@export
#'
read.pfam2go <- function(file) {
  # read file line by line
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
  # change colnames
  names(data) <- c("id","name","go_name","go_id")
  return(data)
}
#' read.go2slim
#'
#' wrapper for read pfam2go mapping files provided by GO Consortium
#' 
#'@param file           mapping input file
#'
#'@return data.frame
#'@export
#'
read.go2slim <- function(file) {
  
  # read in mapping file in obo format 
  obo <- readLines(file)
  
  ## empty list to store results
  obo1 <- vector("list")
  # counter for list elements
  z <- 0
  
  ## get all terms from lines
  n <- get.terms(obo)
  if(length(n)==0) stop("No terms to parse")
  # add last postion in obo file (+2 )
  n <- c(n, length(obo)+2 )
  
  for(i in 1:(length(n) - 1) ){
    ## iterate through the terms
    if(obo[n[i]] == "[Term]"){
      # extract the complete Term
      x <- obo[(n[i] + 1):(n[i + 1] - 2)]
      # ids always first, name second?
      id <- get.ID(x)
      name <- get.Name(x)
      namespace <- get.Namespace(x)
      ## get all tag-value pairs for a single stanza
      p <- x[grep("^is_a|^relationship", x)]
      ## check if missing = root
      if(length(p) == 0){
        z <- z + 1
        obo1[[z]] <- c(id, name, namespace, "root", id)
      } else {
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
          obo1[[z]] <- c(id, name, namespace, rship, pid)
        }
      }
    }
  }
  obo2 <- data.frame(do.call("rbind", obo1), stringsAsFactors = FALSE)
  colnames(obo2) <- c("id", "name", "ontology", "related", "parent_id")
  
  linage <- lapply(1:nrow(obo2),function(x) get.linage(x,obo2))
  linage <- do.call("rbind",linage)
  
  obo3 <- data.frame("id" = obo2$id,
                     "name" = obo2$name,
                     "ontology" = obo2$ontology,
                     "linage" = linage)
  # remove switching ontologies
  obo3 <- obo3[-grep("NA", obo3$linage), ]
  
  
  class(obo2)<-c("obo", "data.frame")
  attr(obo2, "created") <- date
  attr(obo2, "downloaded") <- Sys.Date()
  print(paste("Loaded", length(unique(obo2$id)), "unique terms and", nrow(obo2), "rows"))
  obo2
}