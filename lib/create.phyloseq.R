

create.otu_table <- function(pfamList) {
    
    # get unique identifier for columns
    sample.name <- as.character(get.metadata()[,"SampleName"])
    # create list of data.frames with abundance values for every sample
    matList <- lapply(pfamList, function(x) {
        data <- aggregate(x$count, by = list(x$go.id), sum)
        names(data) <- c("go.id", as.character(unique(x$sample)))
        data
    })
    # combine the abundance values to one data.frame
    otumat <- join_all(matList, by = "go.id", type = 'full')
    # change NA values to zero
    otumat[is.na(otumat)] <- 0
    # extract GO:ID's
    rows <- str_trim(as.character(otumat$go.id))
    # remove first row with GO:ID's and change type to matrix
    otumat <- as.matrix(otumat[, -1])
    # change row.names to GO:ID's
    row.names(otumat) <- rows
    # change colnames to unique sample names
    colnames(otumat) <- sample.name
    # create otu_table in phyloseq format
    otu_table(otumat, taxa_are_rows = T)
}

create.sample_data <- function(metadata) sample_data(metadata)

create.tax_table <- function(linages) {
    taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(otumat), ncol = 7)
    rownames(taxmat) <- rownames(otumat)
    colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", 
                          "Species")
    taxmat
}
generate.phyloseq <- function() {
    otu_tbl <- create.otu_table(list(pfam.60.annotated,pfam.64.annotated, pfam.70.annotated))
    sample_tbl <- create.sample_data(get.metadata())
    tax_tbl <- create.tax_table()
    phyloseq(otu_tbl,sample_tbl,tax_tbl)
}