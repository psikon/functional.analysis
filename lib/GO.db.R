#' init.GO.db
#'
#' load the GO.db package and return entries as list
#' 
#'@export
#'
init.GO.db <- function() {
  require(GO.db)
  as.list(GOTERM)
}
#' get.GOid
#'
#' return GO id for a given GO id
#' 
#'@param id     go id
#'@param GO.db  GO.db object for faster access
#'
#'@return go id
#'@export
#'
get.GOid <- function(id, GO.db = NULL) {
  if(is.null(GO.db)) GO.db <- init.GO.db()
  GOID(GO.db[[grep(str_trim(id), names(GO.db))]])
}
#' get.GOTerm
#'
#' return GO Term for a given GO id
#' 
#'@param id     go id
#'@param GO.db  GO.db object for faster access
#'
#'@return go term
#'@export
#'
get.GOterm <- function(id, GO.db = NULL) {
  if(is.null(GO.db)) GO.db <- init.GO.db()
  Term(GO.db[[grep(str_trim(id), names(GO.db))]])
}

#' get.GOsynonym
#'
#' return GO synonyms for a given GO id
#' 
#'@param id     go id
#'@param GO.db  GO.db object for faster access
#'
#'@return go synonym
#'@export
#'
get.GOsynonym <- function(id, GO.db = NULL) {
  if(is.null(GO.db)) GO.db <- init.GO.db()
  Term(GO.db[[grep(str_trim(id), names(GO.db))]])
}
#' get.GOdefinition
#'
#' return GO definition for a given GO id
#' 
#'@param id     go id
#'@param GO.db  GO.db object for faster access
#'
#'@return go definition
#'@export
#'
get.GO.definition <- function(id, GO.db = NULL) {
  if(is.null(GO.db)) GO.db <- init.GO.db()
  Definition(GO.db[[grep(str_trim(id), names(GO.db))]])
}

#' get.GOontology
#'
#' return GO ontology for a given GO id
#' 
#'@param id     go id
#'@param GO.db  GO.db object for faster access
#'
#'@return go ontology
#'@export
#'
get.GOontology <- function(id, GO.db = NULL) {
  if(is.null(GO.db)) GO.db <- init.GO.db()
  index <- grep(str_trim(id), names(GO.db))
  if (length(index) > 0) {
      ontology <- Ontology(GO.db[[index]])
  } else {
      ontology <- NA_character_
  }
  return(ontology)
}

#' get.GOancestor
#'
#' return the GO ids for all ancestors up to the ontology level 
#' for a given GO id
#' 
#'@param id           go id
#'@param ancestor.db  GO.db ancestor object for faster access
#'
#'@return vector of ancestor strings
#'@export
#'
get.GOancestor <- function(id, ancestor.db = NULL) {
    if (is.null(ancestor.db)) {
      ancestor.db <- list(mf = as.list(GOMFANCESTOR),
                          bp = as.list(GOBPANCESTOR),
                          cc = as.list(GOCCANCESTOR))
    }
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
  return(as.vector(do.call(rbind, ancestor_list)))
}