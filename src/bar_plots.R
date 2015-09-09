library(rmisc)
require.all("ggplot2", "grid", "gridExtra", "scales")
source("src/load_functional.R")
# in zeile 46 hinzufügen
# colnames(mdf)[3] <- "Abundance"
fixInNamespace(psmelt, pos = "package:phyloseq")

#creat copy
func <- functional

# filter
func <- prune_taxa(taxa_sums(func) > 1, func)
func <- subset_taxa(func, ontology != "cellular_component")

rel_abundance <- function(x, n) (x/n)/(sum(x)/n)

free.bp <- subset_taxa(func, ontology == "biological_process")
free.bp <- subset_samples(free.bp, HoldingCondition == "free living")
free.bp <- transform_sample_counts(free.bp, function(OTU) OTU/sum(OTU)*100)
free.bp <- tax_glom(free.bp, "level1")
free.bp <- prune_taxa(taxa_sums(free.bp)/nsamples(free.bp) >= 1, free.bp)
free.bp

tax_count <- table(tax_table(free.bp)[, "level1"])
tax_count <- tax_count[as.vector(tax_table(free.bp)[, "level1"])]
free.bp.df <- data.frame(
  Condition = "free living",
  Abundance = rel_abundance(taxa_sums(free.bp), nsamples(free.bp)),
  nOTU = tax_count,
  Functions = as.vector(tax_table(free.bp)[, "level1"])
)
free.bp.df[,c("Abundance","Functions")]

mari.bp <- subset_taxa(func, ontology == "biological_process")
mari.bp <- subset_samples(mari.bp, HoldingCondition == "mariculture")
mari.bp <- transform_sample_counts(mari.bp, function(OTU) OTU/sum(OTU)*100)
mari.bp <- tax_glom(mari.bp, "level1")
mari.bp <- prune_taxa(taxa_sums(mari.bp)/nsamples(mari.bp) >= 1, mari.bp)
mari.bp

tax_count <- table(tax_table(mari.bp)[, "level1"])
tax_count <- tax_count[as.vector(tax_table(mari.bp)[, "level1"])]
mari.bp.df <- data.frame(
  Condition = "mariculture",
  Abundance = rel_abundance(taxa_sums(mari.bp), nsamples(mari.bp)),
  nOTU = tax_count,
  Functions = as.vector(tax_table(mari.bp)[, "level1"])
)
mari.bp.df

all.bp <- subset_taxa(func, ontology == "biological_process")
all.bp <- transform_sample_counts(all.bp, function(OTU) OTU/sum(OTU)*100)
all.bp <- tax_glom(all.bp, "level1")
all.bp <- prune_taxa(taxa_sums(all.bp)/nsamples(all.bp) >= 1, all.bp)
all.bp

tax_count <- table(tax_table(all.bp)[, "level1"])
tax_count <- tax_count[as.vector(tax_table(all.bp)[, "level1"])]
all.bp.df <- data.frame(
  Condition = "all",
  Abundance = rel_abundance(taxa_sums(all.bp), nsamples(all.bp)),
  nOTU = tax_count,
  Functions = as.vector(tax_table(all.bp)[, "level1"])
)
#all.bp.df <- all.bp.df[all.bp.df$Abundance > 0.01, ]
#all.bp.df$Abundance <- all.bp.df$Abundance * 100
#all.bp.df[,c("Abundance","Functions")]



df.bp <- rbind(free.bp.df, mari.bp.df, all.bp.df)
df.bp <- df.bp[df.bp$Abundance > 0.01, ]
df.bp$Abundance <- df.bp$Abundance * 100
df.bp$Condition <- c(rep("free-living",7), 
                     rep("mariculture",5),
                     rep("shared",6))
bp <- ggplot(df.bp, aes(x = factor(1), y = Abundance, fill = Functions)) + 
  geom_bar(stat = "identity") +
  facet_grid(facets = . ~ Condition) + 
  coord_polar(theta = "y") +
  xlab('') + ylab('') + ggtheme_core +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")
bp
tiff("graphs/pie_chart.bp.tiff", compression = "lzw", ,
     width = 10, height = 3, units = "in",res = 600)
  bp
dev.off()

jpeg("graphs/pie_chart.bp.jpg", width = 10, height = 3, 
     units = "in", res = 300)
bp
dev.off()

free.mf <- subset_taxa(func, ontology == "molecular_function")
free.mf <- subset_samples(free.mf, HoldingCondition == "free living")
free.mf <- transform_sample_counts(free.mf, function(OTU) OTU/sum(OTU)*100)
free.mf <- tax_glom(free.mf, "level1")
free.mf <- prune_taxa(taxa_sums(free.mf)/nsamples(free.mf) >= 1, free.mf)
free.mf

tax_count <- table(tax_table(free.mf)[, "level1"])
tax_count <- tax_count[as.vector(tax_table(free.mf)[, "level1"])]
free.mf.df <- data.frame(
  Condition = "free living",
  Abundance = rel_abundance(taxa_sums(free.mf), nsamples(free.mf)),
  nOTU = tax_count,
  Functions = as.vector(tax_table(free.mf)[, "level1"])
)
free.mf.df

mari.mf <- subset_taxa(func, ontology == "molecular_function")
mari.mf <- subset_samples(mari.mf, HoldingCondition == "mariculture")
mari.mf <- transform_sample_counts(mari.mf, function(OTU) OTU/sum(OTU)*100)
mari.mf <- tax_glom(mari.mf, "level1")
mari.mf <- prune_taxa(taxa_sums(mari.mf)/nsamples(mari.mf) >= 1, mari.mf)
mari.mf

tax_count <- table(tax_table(mari.mf)[, "level1"])
tax_count <- tax_count[as.vector(tax_table(mari.mf)[, "level1"])]
mari.mf.df <- data.frame(
  Condition = "mariculture",
  Abundance = rel_abundance(taxa_sums(mari.mf), nsamples(mari.mf)),
  nOTU = tax_count,
  Functions = as.vector(tax_table(mari.mf)[, "level1"])
)
mari.mf.df

all.mf <- subset_taxa(func, ontology == "molecular_function")
all.mf <- transform_sample_counts(all.mf, function(OTU) OTU/sum(OTU)*100)
all.mf <- tax_glom(all.mf, "level1")
all.mf <- prune_taxa(taxa_sums(all.mf)/nsamples(all.mf) >= 1, all.mf)
all.mf

tax_count <- table(tax_table(all.mf)[, "level1"])
tax_count <- tax_count[as.vector(tax_table(all.mf)[, "level1"])]
all.mf.df <- data.frame(
  Condition = "all",
  Abundance = rel_abundance(taxa_sums(all.mf), nsamples(all.mf)),
  nOTU = tax_count,
  Functions = as.vector(tax_table(all.mf)[, "level1"])
)
#all.mf.df <- all.mf.df[all.mf.df$Abundance > 0.01, ]
#all.mf.df$Abundance <- all.mf.df$Abundance * 100
#all.mf.df[,c("Abundance","Functions")]


df.mf <- rbind(free.mf.df, mari.mf.df, all.mf.df)
df.mf <- df.mf[df.mf$Abundance > 0.01, ]
df.mf$Abundance <- df.mf$Abundance * 100
mf <- ggplot(df.mf, aes(x = factor(1), y = Abundance, fill = Functions)) + 
  geom_bar(stat = "identity") +
  facet_grid(facets = . ~ Condition) + 
  coord_polar(theta = "y") +
  xlab('') + ylab('') + ggtheme_core +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")
mf
tiff("graphs/pie_chart.mf.tiff", compression = "lzw", ,
     width = 8, height = 3, units = "in",res = 600)
mf
dev.off()

jpeg("graphs/pie_chart.mf.jpg", width = 8, height = 3, 
     units = "in", res = 300)
mf
dev.off()

tiff("graphs/pie_chart.combined.tiff", compression = "lzw", ,
     width = 8, height = 6.5, units = "in",res = 600)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1)))
  print(bp, vp = viewport(layout.pos.row = 1,
                                  layout.pos.col = 1))
  print(mf, vp = viewport(layout.pos.row = 2,
                                   layout.pos.col = 1))
dev.off()

jpeg("graphs/pie_chart.combined.jpg", width = 8, 
     height = 6.5, units = "in", res = 300)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1)))
  print(bp, vp = viewport(layout.pos.row = 1,
                        layout.pos.col = 1))
  print(mf, vp = viewport(layout.pos.row = 2,
                        layout.pos.col = 1))
dev.off()

#########################
### detailed Overview ###
#########################

# biological process
symbiosis <- func
symbiosis <- transform_sample_counts(symbiosis, function(OTU) OTU/sum(OTU)*100)
symbiosis <- subset_taxa(symbiosis, ontology == "biological_process")
symbiosis <- subset_taxa(symbiosis, level1 == "symbiosis")
symbiosis <- prune_taxa(taxa_sums(symbiosis)/nsamples(symbiosis) >= 0.1, symbiosis)
sym_list <- otu_table(symbiosis)
rowSums(sym_list)/12
bp <- as.list(GOBPANCESTOR)
sym_list <- sym_list[1:4]
round(sort(taxa_sums(symbiosis), decreasing = TRUE)/sum(taxa_sums(symbiosis))*100, 4)

nitrogen <- func
nitrogen <- subset_taxa(nitrogen, ontology == "biological_process")
nitrogen <- transform_sample_counts(nitrogen, function(OTU) OTU/sum(OTU)*100)
nitrogen <- subset_taxa(nitrogen, level1 == "cellular nitrogen compound metabolic process")
nitrogen <- prune_taxa(taxa_sums(nitrogen)/nsamples(nitrogen) >= 1, nitrogen)
round(sort(taxa_sums(nitrogen), decreasing = TRUE)/sum(taxa_sums(nitrogen))*100, 4)


bioproc <- func
bioproc <- subset_taxa(bioproc, ontology == "biological_process")
bioproc <- transform_sample_counts(bioproc, function(OTU) OTU/sum(OTU)*100)
bioproc <- subset_taxa(bioproc, level1 == "biosynthetic process")
bioproc <- prune_taxa(taxa_sums(bioproc)/nsamples(bioproc) >= 0.1, bioproc)
otu_table(bioproc)
round(sort(taxa_sums(bioproc), decreasing = TRUE)/sum(taxa_sums(bioproc))*100, 4)

loco <- func
loco <- subset_taxa(loco, ontology == "biological_process")
loco <- transform_sample_counts(loco, function(OTU) OTU/sum(OTU)*100)
loco <- subset_taxa(loco, level1 == "locomotion")
loco <- prune_taxa(taxa_sums(loco)/nsamples(loco) >= 1, loco)
otu_table(loco)
round(sort(taxa_sums(loco), decreasing = TRUE)/sum(taxa_sums(loco))*100, 4)

protein <- func
protein <- subset_taxa(protein, ontology == "biological_process")
protein <- transform_sample_counts(protein, function(OTU) OTU/sum(OTU)*100)
protein <- subset_taxa(protein, level1 == "protein maturation")
protein <- prune_taxa(taxa_sums(protein)/nsamples(protein) >= 1, protein)
otu_table(protein)
round(sort(taxa_sums(protein), decreasing = TRUE)/sum(taxa_sums(protein))*100, 4)

transport <- func
transport <- subset_taxa(transport, ontology == "biological_process")
transport <- transform_sample_counts(transport, function(OTU) OTU/sum(OTU)*100)
transport <- subset_taxa(transport, level1 == "transport")
transport <- prune_taxa(taxa_sums(transport)/nsamples(transport) >= 0.1, transport)
round(sort(taxa_sums(transport), decreasing = TRUE)/sum(taxa_sums(transport))*100, 4)

colSums(otu_table(transport))

small_mol <- func
small_mol <- subset_taxa(small_mol, ontology == "biological_process")
small_mol <- transform_sample_counts(small_mol, function(OTU) OTU/sum(OTU)*100)
small_mol <- subset_taxa(small_mol, level1 == "small molecule metabolic process")
small_mol <- prune_taxa(taxa_sums(small_mol)/nsamples(small_mol) >= 0.1, small_mol)
small_mol
round(sort(taxa_sums(small_mol), decreasing = TRUE)/sum(taxa_sums(small_mol))*100, 4)

# molecular function

ion <- func
ion <- subset_taxa(ion, ontology == "molecular_function")
ion <- transform_sample_counts(ion, function(OTU) OTU/sum(OTU)*100)
ion <- subset_taxa(ion, level1 == "ion binding")
ion <- prune_taxa(taxa_sums(ion)/nsamples(ion) >= 1, ion)
otu_table(ion)

oxi <- func
oxi <- subset_taxa(oxi, ontology == "molecular_function")
oxi <- transform_sample_counts(oxi, function(OTU) OTU/sum(OTU)*100)
oxi <- subset_taxa(oxi, level1 == "oxidoreductase activity")
oxi <- prune_taxa(taxa_sums(oxi)/nsamples(oxi) >= 1, oxi)
oxi

trans <- func
trans <- subset_taxa(trans, ontology == "molecular_function")
trans <- transform_sample_counts(trans, function(OTU) OTU/sum(OTU)*100)
trans <- subset_taxa(trans, level1 == "transmembrane transporter activity")
trans <- prune_taxa(taxa_sums(trans)/nsamples(trans) >= 0.1, trans)
otu_table(trans)

kinase <- func
kinase <- subset_taxa(kinase, ontology == "molecular_function")
kinase <- transform_sample_counts(kinase, function(OTU) OTU/sum(OTU)*100)
kinase <- subset_taxa(kinase, level1 == "kinase activity")
kinase <- prune_taxa(taxa_sums(kinase)/nsamples(kinase) >= 0.1, kinase)
otu_table(kinase)

pep <- func
pep <- subset_taxa(pep, ontology == "molecular_function")
pep <- transform_sample_counts(pep, function(OTU) OTU/sum(OTU)*100)
pep <- subset_taxa(pep, level1 == "peptidase activity")
pep <- prune_taxa(taxa_sums(pep)/nsamples(pep) >= 0.1, pep)
otu_table(pep)

enz <- func
enz <- subset_taxa(enz, ontology == "molecular_function")
enz <- transform_sample_counts(enz, function(OTU) OTU/sum(OTU)*100)
enz <- subset_taxa(enz, level1 == "enzyme regulator activity")
enz <- prune_taxa(taxa_sums(enz)/nsamples(enz) >= 0.1, enz)
otu_table(enz)

sig <- func
sig <- subset_taxa(sig, ontology == "molecular_function")
sig <- transform_sample_counts(sig, function(OTU) OTU/sum(OTU)*100)
sig <- subset_taxa(sig, level1 == "signal transducer activity")
sig <- prune_taxa(taxa_sums(sig)/nsamples(sig) >= 0.1, sig)
otu_table(sig)

trans <- func
trans <- subset_taxa(trans, ontology == "molecular_function")
trans <- transform_sample_counts(trans, function(OTU) OTU/sum(OTU)*100)
trans <- subset_taxa(trans, level1 == "transferase activity")
trans <- prune_taxa(taxa_sums(trans)/nsamples(trans) >= 0.1, trans)
otu_table(trans)
