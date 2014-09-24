#' pfam2go
#'
#' assign for pfam ids the corresponding go ids, go ontologies and 
#' ancestor strings from the mapping object
#' 
#' @description convert a pfam count table into a new data.frame 
#'              with the following fields:
#'                sample.id - identifier of sample
#'                db        - origin of annotation
#'                pfam      - pfam id of the protein
#'                name      - name of the found pfam id
#'                go.id     - GO:ID from resulting from mapping
#'                count     - number of occurences of the pfam id
#'                ontology  - ontology of the GO:Id (one of BP, MF or CC)
#'                ancestor  - string with all ancestors of the GO:Id
#'                
#'@param pfam.table       pfam count table
#'@param sample.id        id of the current sample
#'@param pfam2go.mapping  pfam2go mapping object
#'@param GO.db            GO.db object for faster access
#'
#'@return data.frame
#'@export
pfam2GO <- function(pfam.table, 
                    sample.id, 
                    pfam2go.mapping, 
                    GO.db = NULL) {
  
  # load requiered Databases from GO.db package - 
  # for speedup the look up process
  if (is.null(GO.db)) GO.db <- init.GO.db()
  # create the ancestor db
  ancestor.db <- list(mf = as.list(GOMFANCESTOR),
                      bp = as.list(GOBPANCESTOR),
                      cc = as.list(GOCCANCESTOR))
  
  data <- pblapply(pfam.table$pfam, function(x) {
    # find mapped GO:Ids for pfam id
    go <- pfam2go.mapping[grep(x, as.character(pfam2go.mapping$id)), ]
    # if go is empty -> list element is NULL else:
    if(nrow(go) > 0) {
      # create data.frame with supplementary data 
      # ontology and ancestors 
      supplemental <- do.call(rbind.data.frame,
                              lapply(go$go_id, function(x) {
                                data.frame("ontology" = get.GOontology(x, GO.db),
                                           "ancestor" = get.GOancestor(x, ancestor.db))
                                }))
      # create data.frame from go mapping and combine it 
      # with supplementary data.frame and count information
      res <- cbind(data.frame(
        "sample" = as.character(sample.id),
        "db" = as.character("pfam"),
        "pfam" = as.character(go$id),
        "name" = as.character(go$name),       
        "go.id" = as.character(go$go_id),
        "count" = pfam.table[grep(as.character(x)[1],
                                  pfam.table$pfam),
                             "count"]),
        supplemental)
      }
  })
  # remove null elements from list and combine all results 
  # to one data.frame
  data <- do.call(rbind.data.frame, data[!sapply(data, is.null)])
  return(data)
}

# mapGO2slim <- function(annotated.pfam, slim.annotation) {
#     slim <- lapply(annotated.pfam$go_id, function(x) {
#                       x <- str_trim(as.character(x))
#                       ancestors <- get.Ancestor(x)
#                       slim.ids <- ancestors[which(ancestors %in% slim.annotation$id)]     
#                       if (is.null(slim.ids)) NA else slim.ids
#                       })
#     cbind(annotated.pfam,do.call(rbind,slim))
# }

