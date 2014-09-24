#' remap.cluster2pfam
#'
#' wrapper for python script remap_cluster2pfam
#' 
#'@description function calls the python script remap_cluster2pfam.py,
#'             that remaps the reads, that have been clustered together 
#'             during annotation phase to the pfam ids of the 
#'             representive sequences
#' 
#'@param pfam.file          pfam table input file
#'@param clstr.file         clstr file of cd-hit
#'@param output.file        output file
#'
#'@return message
#'@export
#'
remap.cluster2pfam <- function(pfam.file, clstr.file, output.file) {
    # create cmd string for commandline call 
    cmd <- paste("python src/remap_cluster2pfam.py -p",
                 pfam.file, "-c", clstr.file, "-o", 
                 output.file, sep = " ")
    # call the python script
    system(cmd, intern = T, ignore.stdout = T, ignore.stderr = T)
    # return status message
    return(message(paste("remap successfull: ", output.file, sep = " ")))
}