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
         "data/pfam/62.pfam.hmmsearch",
         "data/pfam/64.pfam.hmmsearch",
         "data/pfam/66.pfam.hmmsearch",
         "data/pfam/68.pfam.hmmsearch",
         "data/pfam/70.pfam.hmmsearch",
         "data/pfam/72.pfam.hmmsearch",
         "data/pfam/74.pfam.hmmsearch",
         "data/pfam/76.pfam.hmmsearch",
         "data/pfam/78.pfam.hmmsearch",
         "data/pfam/80.pfam.hmmsearch",
         "data/pfam/82.pfam.hmmsearch")
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
         "data/clstr_files/clustered.62.clstr",
         "data/clstr_files/clustered.64.clstr",
         "data/clstr_files/clustered.66.clstr",
         "data/clstr_files/clustered.68.clstr",
         "data/clstr_files/clustered.70.clstr",
         "data/clstr_files/clustered.72.clstr",
         "data/clstr_files/clustered.74.clstr",
         "data/clstr_files/clustered.76.clstr",
         "data/clstr_files/clustered.78.clstr",
         "data/clstr_files/clustered.80.clstr",
         "data/clstr_files/clustered.82.clstr")
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
         "data/pfam_remap/clustered.62.remap.hmmsearch",
         "data/pfam_remap/clustered.64.remap.hmmsearch",
         "data/pfam_remap/clustered.66.remap.hmmsearch",
         "data/pfam_remap/clustered.68.remap.hmmsearch",
         "data/pfam_remap/clustered.70.remap.hmmsearch",
         "data/pfam_remap/clustered.72.remap.hmmsearch",
         "data/pfam_remap/clustered.74.remap.hmmsearch",
         "data/pfam_remap/clustered.76.remap.hmmsearch",
         "data/pfam_remap/clustered.78.remap.hmmsearch",
         "data/pfam_remap/clustered.80.remap.hmmsearch",
         "data/pfam_remap/clustered.82.remap.hmmsearch")
}