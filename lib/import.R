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

