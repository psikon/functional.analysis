#'@title clean.pfamId
#'
#'get a vector pfam id as input and remove entries after '.'
#'
#'@return character
clean.pfamId <- function(pfam) {
  pfam <- unlist(lapply(pfam, function(x){
    # split the string at the '.' and return only the first element
    strsplit(x,"\\.")[[1]][1]
  }))
  return(pfam)
}



#'@title calculate.PfamAbundance
#'
#' count the occurences of pfam ids from an imported pfam file
#' the output consists of two fields:
#'    \item pfam id - clean name of the pfam id
#'    \item count   - occurences of the pfam id
#'@return data.frame
#'
calculate.abundance <- function(pfam.table) {
  # count pfam ids and create a table
  raw_counts = rle(as.vector(pfam.table[, "pfam"]))
  # convert the table in a data.frame
  pfam_counts = data.frame("pfam" = clean.pfamId(raw_counts$values),
                           "count" = raw_counts$lengths)
  # make the fields type safe
  pfam_counts$pfam <- as.character(pfam_counts$pfam)
  pfam_counts$count <- as.integer(pfam_counts$count)
  return(pfam_counts)
}