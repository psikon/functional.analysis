remap.Clstr2Pfam <- function(pfam.file, clstr.file, output.file) {
    cmd <- paste("python src/remap_cluster2pfam.py -p",
                 pfam.file, "-c", clstr.file, "-o", output.file, sep = " ")
    system(cmd, intern = T, ignore.stdout = T,ignore.stderr = T)
    return(paste("remap successfull: ", output.file, sep = " "))
}