library(ProjectTemplate)
load.project()

#source("http://bioconductor.org/biocLite.R")
#biocLite("GSEABase")

# load GO database from GO.db package
GO.db <- init.GO.db() 
# load pfam2GO mapping from interpro 
pfam2go <- read.pfam2go("data/mapping/pfam2go.obo")

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

pfam.60 <- read.pfam("data/pfam_remap/clustered.60.remap.hmmsearch") 
pfam.64 <- read.pfam("data/pfam_remap/clustered.64.remap.hmmsearch")
pfam.70 <- read.pfam("data/pfam_remap/clustered.70.remap.hmmsearch")

# calculate abundance
pfam.60.count <- calculate.abundance(pfam.60)
pfam.64.count <- calculate.abundance(pfam.64)
pfam.70.count <- calculate.abundance(pfam.70)

# create GO annotation table for pfam

pfam.60.annotated <- pfam2GO(pfam.60.count, sample.id = 60,
                             pfam2go, GO.db)
pfam.64.annotated <- pfam2GO(pfam.64.count, sample.id = 64,
                             pfam2go, GO.db)
pfam.70.annotated <- pfam2GO(pfam.70.count, sample.id = 70,
                             pfam2go, GO.db)

load.project()
pfamList <- 
create.sample_data(get.metadata())
create.otu_table(list(pfam.60.annotated,
                      pfam.64.annotated,
                      pfam.70.annotated))


