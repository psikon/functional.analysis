#' pfam2go
#'
#' assign for pfam ids the corresponding go ids, go ontologies and 
#' ancestor strings from the mapping object
#' 
#' @description convert a pfam count table into a new data.frame 
#'              with the following fields:
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
pfam2GO <- function(pfam.count, 
                    pfam2go.mapping, 
                    GO.db = NULL) {
  
  # load requiered Databases from GO.db package - 
  # for speedup the look up process
  if (is.null(GO.db)) GO.db <- init.GO.db()
  # create the ancestor db
  ancestor.db <- list(mf = as.list(GOMFANCESTOR),
                      bp = as.list(GOBPANCESTOR),
                      cc = as.list(GOCCANCESTOR))
  
  data <- pblapply(pfam.count$pfam, function(x) {
    # find mapped GO:Ids for pfam id
    go <- pfam2go.mapping[grep(x, as.character(pfam2go.mapping$id)), ]
    # if go is empty -> list element is NULL else:
    if(nrow(go) > 0) {
      # create data.frame with supplementary data 
      # ontology and ancestors 
      supplemental <- do.call(rbind.data.frame,lapply(go$go_id, 
                function(x) {
                  data.frame("ontology" = get.GOontology(x, GO.db),
                  "ancestor" = get.GOancestor(x, ancestor.db))
                                }))
      # create data.frame from go mapping and combine it 
      # with supplementary data.frame and count information
      res <- cbind(data.frame(
        "db" = as.character("pfam"),
        "pfam" = as.character(go$id),
        "name" = as.character(go$name),       
        "go.id" = as.character(go$go_id),
        "count" = pfam.count[grep(as.character(x)[1],
                                  pfam.count$pfam),
                             "count"], 
        stringsAsFactors = F),
        supplemental)
      }
  })
  # remove null elements from list and combine all results 
  # to one data.frame
  data <- do.call(rbind.data.frame, data[!sapply(data, is.null)])
  return(data)
}
#' go2slim
#'
#' assign a linage basing of the slim annotation of gene ontology to
#' the go.ids 
#' 
#' @description convert the annotated go id table to a new data.frame
#'              with following fields
#'                pfam      - pfam id of the protein
#'                go.id     - GO:ID from resulting from mapping
#'                count     - number of occurences of the pfam id
#'                ontology  - ontology of the GO:Id (one of BP, MF or CC)
#'                linage    - linage from the go slim annotation
#'                
#'@param pfam.go          go annotated pfam table
#'@param go2slim.mapping  go2slim mapping object
#'
#'@return data.frame
#'@export
go2slim <- function(pfam.go, go2slim.mapping) {
    # assign slim linage to every row
    slim.linage <- apply(pfam.go, 1, function(x) {
      # get ancestor
      ancestor <- x["ancestor"]
      # convert string to vector
      ancestor <- unlist(str_split(ancestor,","))
      # find index of last hit 
      idx <- max(which(ancestor %in% go2slim.mapping$id == TRUE))
      # get annotation mapping
      slim <- go2slim.mapping[which(go2slim.mapping$id %in% ancestor[idx]),]
      # only assign one linage to go id
      slim$linage[1]
    })
    # build new data.frame with linage
    slim.annotation <- data.frame(pfam.go$pfam,
                                  pfam.go$go.id, 
                                  pfam.go$count,
                                  pfam.go$ontology,
                                  as.vector(slim.linage),
                                  stringsAsFactors = F)
    # adjust column names
    colnames(slim.annotation) <- c("pfam.id", "go.id",
                                   "count", "ontology",
                                   "slim.linage")
    return(slim.annotation)    
}

getGOfromSlim <- function(functions, slim_list) {
  index <- unlist(lapply(functions, FUN = function(x) {
    grep(x, slim_list[[1]]$slim.linage)
  }))
  data <- slim_list[[1]][index, ]
  return(data)
}

aggregateGObySlim <- function(functions) {
  data <- lapply(functions$go.id, function(x){
    data <- functions[functions$go.id == x,]
    data.frame(aggregate(count ~ go.id, data = data, FUN = sum),
               slim.linage = unique(functions$slim.linage))
  })
  data <- unique(do.call(rbind, data))
}



