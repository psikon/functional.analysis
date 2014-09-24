#' clean.pfamID
#'
#' remove artifacts from the pfam ids to make them comparable to the
#' mapping file
#' 
#'@param pfam           pfam ids
#'
#'@return pfam ids
#'@export
#'
clean.pfamID <- function(pfam) {
  pfam <- unlist(lapply(pfam, function(x){
    # split the string at the '.' and return only the first element
    strsplit(x,"\\.")[[1]][1]
  }))
  return(pfam)
}



#' calculate.abundance
#'
#' count the abundance value for the pfam ids
#' 
#'@param pfam.table       data.frame with pfam ids
#'
#'@return data.frame
#'@export
#'
calculate.abundance <- function(pfam.table) {
    
    # count the occurences of pfam ids
    pfam.count <- as.data.frame(table(pfam.table$pfam))
    
    # adjust names
    names(pfam.count) <- c("pfam", "count")
    
    # make the fields type safe
    pfam.count$pfam <- as.character(pfam.count$pfam)
    pfam.count$count <- as.integer(pfam.count$count)
    pfam.count$pfam <- clean.pfamID(pfam.count$pfam)
    # sort the data.frame descending by occurences
    pfam.count <- arrange(pfam.count, desc(count))
    
    return(pfam.count)
}

