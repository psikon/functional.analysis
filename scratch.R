library(ProjectTemplate)
load.project()
# in zeile 46 hinzufügen
# colnames(mdf)[3] <- "Abundance"
fixInNamespace(psmelt, pos = "package:phyloseq")

#source("http://bioconductor.org/biocLite.R")
#biocLite("GSEABase")

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
# map go annotation to tables
goList <- lapply(countList, function(x) pfam2GO(pfam.count = x,
                                            pfam2go.mapping = pfam2go, 
                                            GO.db = GO.db))
slimList <- lapply(goList, function(x) go2slim(pfam.go = x,
                                               go2slim.mapping = go2slim.map))

functional <- create.phyloseq(slimList)

plot_bar(cc, fill="level1")
