# WallenZD_2022 data generation

#' this dataset comes from the Wallen Nat. Comms paper in 2022 
#' (https://www.nature.com/articles/s41467-022-34667-x)
#' 
#' The processed source data for the analysis is published on Zenodo
#' (https://zenodo.org/records/7246185/), in a .xlsx file called 
#' `Source_Data_24Oct2022.xlsx`
#' 
#' The following code shows the steps needed to generate the source data
local_fileName <- file.path(tempdir(), "wallenZD_zenodo.xlsx")
download.file("https://zenodo.org/records/7246185/files/Source_Data_24Oct2022.xlsx?download=1", destfile = local_fileName)

# destdir
destdir <- "inst/extdata/"

#' 01 - generation of `WallenZD_2022_metaphlan3_profiles.tsv`
#' md5sum: b76a7852a67fbb74d3b54ade931aa0a2
# UNKNOWN/UNCLASSIFIED row gets lost with this method
mpa_profles <- openxlsx2::read_xlsx(local_fileName, sheet = "metaphlan_rel_ab")
write.table(mpa_profles, file = bzfile(file.path(destdir, "WallenZD_2022_metaphlan3_profiles.tsv.bz2")), quote= FALSE, row.names = FALSE, sep = "\t")

#' 02 - generation of `WallenZD_2022_subjMetadata.tsv.bz2`
#' md5sum: 4f702d810a93785c1a61f8b6ccb433f3
subj_metadata <- openxlsx2::read_xlsx(local_fileName, sheet = "subject_metadata") 
write.table(subj_metadata, file = bzfile(file.path(destdir, "WallenZD_2022_subjMetadata.tsv.bz2")), quote= FALSE, row.names = FALSE, sep = "\t")

# create ready-to-use .rda file too
WallenZD_2022.tse <- mia::importMetaPhlAn(
  file = file.path(destdir, "WallenZD_2022_metaphlan3_profiles.tsv.bz2"), 
  col.data = file.path(destdir, "WallenZD_2022_subjMetadata.tsv.bz2"))

# write data
usethis::use_data(WallenZD_2022.tse)
