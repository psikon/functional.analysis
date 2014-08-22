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
  index <- grep(str_trim(id), names(GO.db))
  if (length(index) > 0) {
      ontology <- Ontology(GO.db[[index]])
  } else {
      ontology <- NA_character_
  }
  return(ontology)
}

# return all ancestors for a given GO id concatenated as one string 
get.Ancestor <- function(id, ancestor.db = NULL) {
    if (is.null(ancestor.db)) ancestor.db <- list(mf = as.list(GOMFANCESTOR),
                                                 bp = as.list(GOBPANCESTOR),
                                                 cc = as.list(GOCCANCESTOR))
  # create a list of ancestors in case of multiple id's are given
  ancestor_list <- lapply(id, function(x) {
    # trim whitspace and convert id in character 
    x <- str_trim(as.character(x))
    # check go id for ontology
    if (x %in% names(ancestor.db[["mf"]])) {
      ancestor <- ancestor.db[["mf"]][[grep(x, names(ancestor.db[["mf"]]))]]
      ancestor <- paste(ancestor[-length(ancestor)], collapse = ",")
    } else if (x %in% names(ancestor.db[["bp"]])) {
        ancestor <- ancestor.db[["bp"]][[grep(x, names(ancestor.db[["bp"]]))]]
      ancestor <- paste(ancestor[-length(ancestor)], collapse = ",")
    } else if (x %in% names(ancestor.db[["cc"]])) {
        ancestor <- ancestor.db[["cc"]][[grep(x, names(ancestor.db[["cc"]]))]]
      ancestor <- paste(ancestor[-length(ancestor)],collapse = ",")
    } else {
        ancestor <- x
    }
    ancestor
  })
  # return the list as a vector of ancestor strings
  as.vector(do.call(rbind, ancestor_list))
}