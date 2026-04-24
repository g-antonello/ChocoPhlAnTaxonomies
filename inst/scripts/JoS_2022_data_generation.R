library(parkinsonsMetagenomicData)
library(biobakeryUtils)
library(magrittr)
suppressMessages(library(TreeSummarizedExperiment))

data(sampleMetadata)

# select a study
JoS_metadata.df <- sampleMetadata[sampleMetadata$study_name == "JoS_2022",]

# download metaphlan profiles
JoS_2022.tse <- returnSamples(JoS_metadata.df, data_type = "relative_abundance")

# the enhance function: 
# 1- formats the taxonomy table
# 2- shortens rownames of the TSE object
# 3- calculates some other commonly used/requested data transformations: 
#    relative abundance and counts. MetaPhlAn's default is percent

JoS_2022.tse <- pMD_enhance_MetaPhlAn(JoS_2022.tse, addPhyloTree = FALSE)

# save data
usethis::use_data(JoS_2022.tse)
