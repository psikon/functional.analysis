get.metadata <- function() {
    
    # unique identifier for sample
    sample.id <- c(60,64,70)
    # name of the sample
    sample.name <- c("ef_free_60", "ef_free_64", "ef_cage_70")
    # location of collection
    location <- rep("Jakarta", 3)
    # host
    host <- rep("Epinephelus fuscoguttatus", 3)
    # Holding Condition of sample
    holding.condition = c("free living", "free living", "mariculture")
    # description of sample
    description <- rep("fecal sample", 3)
    # quality of DNA Isolation
    dna.quality <- rep("good", 3)
    # create the data.frame
    data <- data.frame(SampleId = sample.id, SampleName = sample.name, 
                       Location = location, Host = host, 
                       HoldingCondition = holding.condition, 
                       Description = description, DNAQuality = dna.quality)
    # change row.names
    row.names(data)<- sample.name
    # return data
    data
}