#' Load all requiered databases and create from the pfam Files 
#' a GO annotation object with normal and slim annotations.
#' The results will be saved in a cache object for later processing.

# load the functional annotation project
library(ProjectTemplate)
load.project()

# load GO database from GO.db package
GO.db <- init.GO.db() 
# load pfam2GO mapping from interpro
pfam2go <- read.pfam2go("data/mapping/pfam2go.obo")
go2slim.map <- read.go2slim("data/mapping/slim.generic.obo")

# get lists of all needed files
pfam.files <- pfam.fileList()
clstr.files <- clstr.fileList()
remap.files <- remap.fileList()

# create remaped pfam files
for(i in 1:length(pfam.files)) {
  remap.cluster2pfam(pfam.files[[i]],
                     clstr.files[[i]],
                     remap.files[[i]])
}

# import pfam tables
pfamList <- lapply(remap.files, function(x) read.pfam(x))
# calculate abundance for tables
countList <- lapply(pfamList, function(x) calculate.abundance(x))
save(pfamList,countList,file="cache/pfam.import.Rdata")

# map go annotation to tables
goList <- lapply(countList, function(x) pfam2GO(pfam.count = x,
                                                pfam2go.mapping = pfam2go, 
                                                GO.db = GO.db))
# add slim annotations
slimList <- lapply(goList, function(x) go2slim(pfam.go = x,
                                               go2slim.mapping = go2slim.map))
# save as R cache Object
save(goList,slimList, file = "cache/go.annotations.Rdata")