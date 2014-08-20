#' load the GO.db package and return entries as list
init.GO.db <- function() {
  require(GO.db)
  as.list(GOTERM)
}

#' return GO id for given GO id
get.GOid <- function(id, GO.db = NULL) {
  if(is.null(GO.db)) GO.db <- init.GO.db()
  GOID(GO.db[[grep(str_trim(id), names(GO.db))]])
}

#' return GO Term for given GO id
get.Term <- function(id, GO.db = NULL) {
  if(is.null(GO.db)) GO.db <- init.GO.db()
  Term(GO.db[[grep(str_trim(id), names(GO.db))]])
}

#' return Synonyms for given GO id
get.Synonym <- function(id, GO.db = NULL) {
  if(is.null(GO.db)) GO.db <- init.GO.db()
  Term(GO.db[[grep(str_trim(id), names(GO.db))]])
}

#' return definition for given GO id
get.Definition <- function(id, GO.db = NULL) {
  if(is.null(GO.db)) GO.db <- init.GO.db()
  Definition(GO.db[[grep(str_trim(id), names(GO.db))]])
}
#' return Ontology for given GO id
get.Ontology <- function(id, GO.db = NULL) {
  if(is.null(GO.db)) GO.db <- init.GO.db()
  return(Ontology(GO.db[[grep(str_trim(id), names(GO.db))]]))
}

# return all ancestors for a given GO id concatenated as one string 
get.Ancestor <- function(id) {
  # load the three GO ontologies
  mf <- as.list(GOMFANCESTOR)
  bp <- as.list(GOBPANCESTOR)
  cc <- as.list(GOCCANCESTOR)
  # create a list of ancestors in case of multiple id's are given
  ancestor_list <- lapply(id, function(x) {
    # trim whitspace and convert id in character 
    x <- str_trim(as.character(x))
    # check go id for ontology
    if (x %in% names(mf)) {
      ancestor <- mf[[grep(x, names(mf))]]
      ancestor <- paste(ancestor[-length(ancestor)], collapse = ",")
    } else if (x %in% names(bp)) {
      ancestor <- bp[[grep(x, names(bp))]]
      ancestor <- paste(ancestor[-length(ancestor)], collapse = ",")
    } else if (x %in% names(cc)) {
      ancestor <- cc[[grep(x, names(cc))]]
      ancestor <- paste(ancestor[-length(ancestor)],collapse = ",")
    }
    ancestor
  })
  # return the list as a vector of ancestor strings
  as.vector(do.call(rbind, ancestor_list))
}