#' get.metadata
#'
#' wrapper for saving metadata informations for the actual project
#' 
#'
#'@return data.frame
#'@keywords internal
#'
get.metadata <- function() {
    # number of samples
    num.samples = 12
    # unique identifier for sample
    sample.id <- c(60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82)
    # name of the sample
    sample.name <- c("ef_free_60", "ef_free_62", "ef_free_64", 
                     "ef_free_66", "ef_free_68", "ef_cage_70",
                     "ef_cage_72", "ef_cage_74", "ef_cage_76",
                     "ef_cage_78", "ef_cage_80", "ef_cage_82")
    # location of collection
    location <- rep("Jakarta", num.samples)
    # host
    host <- rep("Epinephelus fuscoguttatus", num.samples)
    # Holding Condition of sample
    holding.condition = c(rep("free living", 5), rep("mariculture", 7))
    # description of sample
    description <- rep("fecal sample", num.samples)
    # quality of DNA Isolation
    dna.quality <- rep("good", num.samples)
    # create the data.frame
    data <- data.frame(SampleId = sample.id, SampleName = sample.name, 
                       Location = location, Host = host, 
                       HoldingCondition = holding.condition, 
                       Description = description, DNAQuality = dna.quality)
    # change row.names
    row.names(data)<- sample.name
    # return data
    return(data)
}