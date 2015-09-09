#' Load the cached GO Annotation object an6t init the project itself
#'

library(ProjectTemplate)
load.project()

# load the cached object
load("cache/go.annotations.Rdata")
# create the phyloseq object
functional <- create.phyloseq(slimList)