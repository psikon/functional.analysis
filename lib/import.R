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
  lines <- readLines(file)
  
  ## empty list to store results
  obo.list <- vector("list")
  # counter for list elements
  list.count <- 0
  
  ## get all terms from lines
  term.index <- get.terms(lines)
  if(length(term.index) == 0) stop("No terms to parse")
  
  # add last postion in obo file (+2 )
  term.index <- c(term.index, length(lines) + 2)
  
  for(i in 1:(length(term.index) - 1)) {
    ## iterate through the terms
    if(lines[term.index[i]] == "[Term]"){
      # extract the complete Term
      single.term <- lines[(term.index[i] + 1):(term.index[i + 1] - 2)]
      # find out all go ids with this annotation
      go.ids <- c(get.ID(single.term), get.alternatives(single.term))
      for (k in 1:length(go.ids)) {
        # extract go id
        id <- go.ids[k]
        # extract name
        name <- get.Name(single.term)
        # extract namespace
        namespace <- get.Namespace(single.term)
        # extract parents 
        p <- single.term[grep("^is_a|^relationship", single.term)]
        # check for parent is a root
        if(length(p) == 0){
          list.count <- list.count + 1
          obo.list[[list.count]] <- c(id, name, namespace, "root", id)
        } else {
          # built up relationships for go id
          for(j in 1:length(p)){
            pid <- gsub("(.* )([[:alpha:]|\\_]+:[[:digit:]]+)(.*)", "\\2", p[j])
            # capture relationship type
            rship <- gsub("(.*)( [[:alpha:]|\\_]+:[[:digit:]]+)(.*)", "\\1", p[j])
            rship <- gsub("relationship: ", "", rship)
            rship <- gsub("is_a:", "is_a", rship)
            list.count <- list.count + 1
            obo.list[[list.count]] <- c(id, name, namespace, rship, pid)
          }
        }
      }
    }
  }
  # build data.frame from oboList
  obo.frame <- data.frame(do.call("rbind", obo.list), 
                          stringsAsFactors = FALSE)
  colnames(obo.frame) <- c("id", "name", "ontology", 
                           "related", "parent_id")
  
  linage <- do.call("rbind", lapply(1:nrow(obo.frame), 
                                    function(x) get.linage(x, obo.frame)))
  
  res <- data.frame("id" = obo.frame$id,
                    "name" = obo.frame$name,
                    "ontology" = obo.frame$ontology,
                    "linage" = linage,
                    stringsAsFactors = FALSE)
  # remove switching ontologies
  res <- res[-grep("NA", res$linage), ]
  return(res)
}