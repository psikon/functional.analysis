library(ProjectTemplate)
library(ggplot2)
library(grid)
load.project()

# create a copy
func <- functional 

# number of mapped reads 
sum(colSums(otu_table(func)))
summary(colSums(otu_table(func)))

# filtering
func <- prune_taxa(taxa_sums(func) > 1, func)
func <- subset_taxa(func, ontology != "cellular_component")
# number of filtered reads
sum(colSums(otu_table(func)))
# summary
summary(colSums(otu_table(func)))

# invstigat biological process
bp <- subset_taxa(func, ontology == "biological_process")
bp

unique_bp <- unique(tax_table(bp)[, rank_names(bp)[3]])
bp <- tax_glom(bp, taxrank = rank_names(bp)[3])
unique_bp <- as.vector(tax_table(bp)[order(taxa_sums(bp), 
               decreasing = TRUE),rank_names(bp)[3]])
length(unique_bp)
unique_bp

bp_c <- round(sort(taxa_sums(bp), decreasing = TRUE)/sum(taxa_sums(bp))*100, 4)
bp_c

mf <- subset_taxa(func, ontology == "molecular_function")
mf

unique_mf <- unique(tax_table(mf)[, rank_names(mf)[3]])
mf <- tax_glom(mf, taxrank = rank_names(mf)[3])
unique_mf <- as.vector(tax_table(mf)[order(taxa_sums(mf), 
                  decreasing = TRUE),rank_names(mf)[3]])
length(unique_mf)
unique_mf

mf_c <- round(sort(taxa_sums(mf), decreasing = TRUE)/sum(taxa_sums(mf))*100, 4)
mf_c
