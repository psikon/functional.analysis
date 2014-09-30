#' create.otu_table
#'
#' generate a phyloseq otu_table from a pfam count table
#' 
#'@param pfamList     list of annotated pfam tables
#'
#'@return otu_table
#'@export
#'
create.otu_table <- function(pfamList, names) {
    # get unique identifier for columns
    sample.name <- as.character(get.metadata()[ , "SampleName"])
    # assign sample name to the data.frames
    for(i in 1:length(pfamList)){
      pfamList[[i]][,"sample"] <- sample.name[i]
    }
    # create list of data.frames with abundance values for every sample
    matList <- lapply(pfamList, function(x) {
        data <- aggregate(x$count, by = list(as.character(x$go.id)), sum)
        names(data) <- c("go.id", unique(x$sample))
        data
    })
    # combine the abundance values to one data.frame
    otumat <- join_all(matList, by = "go.id", type = 'full')
    # change NA values to zero
    otumat[is.na(otumat)] <- 0
    # extract GO:ID's
    rows <- str_trim(as.character(otumat$go.id))
    # remove first row with GO:ID's and change type to matrix
    otumat <- as.matrix(otumat[ , -1])
    # change row.names to GO:ID's
    row.names(otumat) <- rows
    # change colnames to unique sample names
    colnames(otumat) <- sample.name
    # create otu_table in phyloseq format
    otu_table(otumat, taxa_are_rows = T)
}

#' create.sample_data
#'
#' convert the metadata to phyloseq sample_data
#' 
#'@param metadata     result of get.metadata function
#'
#'@return sample_table
#'@export
#'
create.sample_data <- function(metadata) sample_data(metadata)

#' create.tax_table
#'
#' generate a tax_table from pfam input data
#' 
#'@param linages     ?
#'
#'@return tax_table
#'@export
#'
create.tax_table <- function(slimList) {
  # get unique identifier for columns
  sample.name <- as.character(get.metadata()[ , "SampleName"])
  # extract only relevant data of the data.frames
  data <- lapply(slimList, function(x) data.frame(x$go.id,x$slim.linage,stringsAsFactors=F))
  # combine to one data.frame with unique rows
  data <- unique(do.call(rbind, data))
  # split up the linage string and add "root" level
  lin <- lapply(data$x.slim.linage, function(x) {
    c("root", unlist(str_split(x,", ")))
  })
  max.level <- max(unlist(lapply(lin,length))) 
  lin <- lapply(lin, function(x){
    c(x,rep(NA,max.level-length(x)))
  })
  taxmat <- as.matrix(do.call("rbind", lin))
  # assign go ids to row
  rownames(taxmat) <- str_trim(data$x.go.id)
  colnames(taxmat) <- c("root","ontology",paste0("level", seq(max.level - 2)))
  return(tax_table(taxmat))
}

#' create.phyloseq
#'
#' build a phyloseq object from the three components
#' 
#'@param pfamList     list of pfam ids
#'
#'@return otu_table
#'@export
#'
create.phyloseq <- function(slimList) {
    otu <- create.otu_table(slimList)
    sample <- create.sample_data(get.metadata())
    tax <- create.tax_table(slimList)
    return(phyloseq(otu,sample,tax))
}