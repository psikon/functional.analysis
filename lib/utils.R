get.terms <- function(lines) {
  # find the structure of "[TERM]" in lines object
  return(grep("^\\[", lines))
}

get.ID <- function(term) {
  return(gsub("(id: )(.*)", "\\2", term[grep("^id:", term)]))
}

get.Name <- function(term) {
  return(gsub("(name: )(.*)", "\\2",  term[grep("*name:", term)]))
}

get.Namespace <- function(term) {
  return(gsub("(namespace: )(.*)", "\\2",  term[grep("*namespace:", term)]))
}

get.alternatives <- function(term) {
  return(gsub("(alt_id: )(.*)", "\\2",  term[grep("alt_id:", term)]))
}

is.root <- function(parent) {
  return(if(length(parent) == 0) TRUE else FALSE)
}

get.linage <- function(index, complete.obo) {
  # extract actual row from data.frame
  data <- complete.obo[index, ]
  # add first element to linage
  linage <- data$name
  # get realtionship
  rship <- data$related
  # get ontology
  onto <- data$ontology
  # build up linage until element is root
  while(rship != "root") {
    # find the parent in data.frame
    data <- complete.obo[grep(data$parent_id, complete.obo$id), ]
    # if parent has different ontology abort and give empty linage
    if (all(unique(data$ontology) != onto)) {
      linage <- NA
      break
    }
    # check if only one relationship in parent
    if (length(unique(data$related)) == 1) rship <- unique(data$related)
    # change 'occurs_in' and 'has_part' to normal 'is_a'
    if(rship == "occurs_in") rship <- "is_a"
    if(rship == "has_part") rship <- "is_a"
    # if parent is root compress it without relationship
    if (any(data$related == "root")) {
      data <- unique(data)
    } else {
      # if not follow path of relationship
      data <- unique(data[grep(rship, data$related), ])
    }
    # add parent to linage
    linage <- append(linage, data$name)
    # adjust relationship
    rship <- data$related
  }
  # return linage as character string in reverse order
  return(toString(rev(linage)))    
}

equalize.linage <- function(linage) {
  data <- str_split(linage, ", ")
  max(do.call(rbind,lapply(data, function(x) length(x))))
  
}
