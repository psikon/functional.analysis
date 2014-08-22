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
    
    # count the occurences of pfam ids
    pfam.count <- as.data.frame(table(pfam.table$pfam))
    
    # adjust names
    names(pfam.count) <- c("pfam", "count")
    
    # make the fields type safe
    pfam.count$pfam <- as.character(pfam.count$pfam)
    pfam.count$count <- as.integer(pfam.count$count)
    pfam.count$pfam <- clean.pfamId(pfam.count$pfam)
    # sort the data.frame descending by occurences
    pfam.count <- arrange(pfam.count,desc(count))
    
    return(pfam.count)
}

