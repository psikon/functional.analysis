library(ProjectTemplate)
load.project()

#source("http://bioconductor.org/biocLite.R")
#biocLite("GSEABase")

# load pfam2GO mapping from interpro 
pfam2go <- import.pfam2go("data/mappings/pfam2go.obo")


pfam.table <- import.pfamTable("data/classify.orf.60.pfam.hmmsearch")

# calculate abundance
pfam.count <- calculate.abundance(pfam.table)

go.db <- init.GO.db() 

# create GO annotation table for pfam
pfam.GO.annotated <- pfam2GO(pfam.count,
                             sample.id = 60,
                             pfam2go,
                             go.db)
save(pfam.GO.annotated,file = "cache/annotated.pfams.Rdata")

#load(file = "cache/annotated.pfams.Rdata")
write.annotated.Pfams(pfam.GO.annotated,"data/annotated.pfams.txt")
data <- import.FromFile("data/annotated.pfams.txt")
