#' pfam.fileList
#'
#' wrapper for saving pfam table files
#' 
#'
#'@return list of pfam table files
#'@keywords internal
#'
pfam.fileList <- function() {
    list("data/pfam/60.pfam.hmmsearch",
         "data/pfam/64.pfam.hmmsearch",
         "data/pfam/70.pfam.hmmsearch")
}

#' clstr.fileList
#'
#' wrapper for saving cd-hit clstr files
#' 
#'
#'@return list of clstr files
#'@export
#'
clstr.fileList <- function() {
    list("data/clstr_files/clustered.60.clstr",
         "data/clstr_files/clustered.64.clstr",
         "data/clstr_files/clustered.70.clstr")
}

#' remap.fileList
#'
#' wrapper for saving remapped pfam files
#' 
#'
#'@return list of remapped pfam files
#'@export
#'
remap.fileList <- function() {
    list("data/pfam_remap/clustered.60.remap.hmmsearch",
         "data/pfam_remap/clustered.64.remap.hmmsearch",
         "data/pfam_remap/clustered.70.remap.hmmsearch")
}