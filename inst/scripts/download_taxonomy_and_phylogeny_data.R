## ---- Initial Setup ----
library(zen4R) 
library(tools) 
library(ape)   

# Paths are relative to the package base directory
mpa_doi <- "https://doi.org/10.5281/zenodo.3965462"
outdir  <- "../tmp/metaphlan_downloads"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## ---- Download latest standard taxonomy from cmprod1 ----
cmprod_url <- "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases"
mpa_latest <- readLines(sprintf("%s/mpa_latest", cmprod_url), warn = FALSE)

tax_dest <- file.path("..", "inst", "extdata", paste0(mpa_latest, "_speciesTaxonomy.tsv.bz2"))
if (!dir.exists(dirname(tax_dest))) dir.create(dirname(tax_dest), recursive = TRUE)

download.file(
  url      = sprintf("%s/%s_species.txt.bz2", cmprod_url, mpa_latest),
  destfile = tax_dest
)

latest_mpa_taxonomy <- read.delim(tax_dest, header = FALSE, col.names = c("SGB", "Taxonomy"))

## ---- List all relevant Zenodo entries of MetaPhlAn ----
zenodo_versions_with_doi <- get_versions(mpa_doi)
zenodo_versions_with_doi <- zenodo_versions_with_doi[, c("date", "version", "doi")]
zenodo_versions_with_doi <- zenodo_versions_with_doi[!grepl("beta|test", zenodo_versions_with_doi$version), ]
zenodo_versions_with_doi <- zenodo_versions_with_doi[!duplicated(zenodo_versions_with_doi$doi), ]
zenodo_versions_with_doi <- zenodo_versions_with_doi[rev(order(zenodo_versions_with_doi$date)), ]
rownames(zenodo_versions_with_doi) <- seq_len(nrow(zenodo_versions_with_doi))

## ---- Download and unzip all versions (~3 min) ----
for (i in seq_len(nrow(zenodo_versions_with_doi))) {
  if (i %% 8 == 0) Sys.sleep(30) # manually tested that a cooldown every 8 downloads is useful
  download_zenodo(doi = zenodo_versions_with_doi$doi[i], path = outdir)
}

invisible(sapply(
  list.files(outdir, pattern = "\\.zip$", full.names = TRUE),
  function(x) unzip(zipfile = x, overwrite = TRUE, exdir = gsub("\\.zip$", "", x))
))

## ---- Download tree exceptions from cmprod ----
# v3.1 corrected tree
v31_dir <- file.path(outdir, "MetaPhlAn-3.1.0")
if (!dir.exists(v31_dir)) dir.create(v31_dir, recursive = TRUE)
download.file(
  "http://cmprod1.cibio.unitn.it/biobakery3/mpa_v31_CHOCOPhlAn_201901_species_tree.nwk",
  destfile = file.path(v31_dir, "mpa_v31_CHOCOPhlAn_201901_species_tree.nwk")
)

# v4.2+ latest tree
v42_dir <- file.path(outdir, "MetaPhlAn-4.2.0")
if (!dir.exists(v42_dir)) dir.create(v42_dir, recursive = TRUE)
download.file(
  "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan25_CHOCOPhlAnSGB_202503.nwk",
  destfile = file.path(v42_dir, "mpa_vJan25_CHOCOPhlAnSGB_202503.nwk")
)

## ---- Build Tree Directory Table ----
trees.vec         <- list.files(outdir, pattern = "\\.nwk$", recursive = TRUE, full.names = TRUE)
trees.md5sum      <- md5sum(trees.vec)
trees.split.list  <- strsplit(trees.vec, "/", fixed = TRUE)
trees.list        <- lapply(trees.vec, read.tree)
names(trees.list) <- basename(trees.vec)

trees.df <- data.frame(
  filePath          = trees.vec,
  tree_filename     = basename(trees.vec),
  MetaPhlAn_version = sapply(trees.split.list, function(x) gsub("MetaPhlAn-", "", x[4], fixed = TRUE)),
  tree_Ntips        = sapply(trees.list, Ntip),
  md5sum            = trees.md5sum,
  stringsAsFactors  = FALSE
)

# Remove the allegedly incorrect v3.1.0 tree (The Zenodo archive ships the v3.0 one)
trees.df <- trees.df[
  !(trees.df$MetaPhlAn_version == "3.1.0" &
      trees.df$tree_filename == "mpa_v30_CHOCOPhlAn_201901_species_tree.nwk"), ]

## ---- Save Processed Trees ----
tree_destdir <- "inst/extdata/trees"
if (!dir.exists(tree_destdir)) dir.create(tree_destdir, recursive = TRUE)

invisible(lapply(trees.df$tree_filename, function(filename) {
  # Get the specific tree object from our list
  tr <- trees.list[[filename]]
  out_path <- file.path(tree_destdir, paste0(filename, ".bz2"))
  write.tree(phy = tr, file = bzfile(out_path))
}))

# Write the manifest
write.table(
  trees.df[!duplicated(trees.df$md5sum),][, -1], # drop filePath
  file      = file.path(tree_destdir, "trees_files_directory.tsv"),
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

## ---- Download MetaPhlAn taxonomies from cmprod1 ----
mpa_tax_dir <- "inst/extdata/mpa_taxonomy/"
if (!dir.exists(mpa_tax_dir)) dir.create(mpa_tax_dir, recursive = TRUE)

tax_urls <- c(
  "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan25_CHOCOPhlAnSGB_202503_species.txt.bz2",
  "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202403_species.txt.bz2",
  "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202307_species.txt.bz2",
  "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212_species.txt.bz2",
  "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202403_species.txt.bz2"
)

for (url in tax_urls) {
  download.file(url, destfile = file.path(mpa_tax_dir, basename(url)))
}

## ---- Format and write GTDB taxonomy tables ----
# Search recursively in the tmp download folder for GTDB files
GTDBs.vec         <- list.files(outdir, pattern = ".*GTDB.*", recursive = TRUE, full.names = TRUE)
GTDBs.md5sum      <- md5sum(GTDBs.vec)
GTDBs.split.list  <- strsplit(GTDBs.vec, "/")

GTDBs.df <- data.frame(
  filePath          = GTDBs.vec,
  GTDB_filename     = basename(GTDBs.vec),
  MetaPhlAn_version = sapply(GTDBs.split.list, function(x) gsub("MetaPhlAn-", "", x[4], fixed = TRUE)),
  md5sum            = GTDBs.md5sum,
  stringsAsFactors  = FALSE
)

GTDBs.df_unique <- GTDBs.df[!duplicated(GTDBs.df$md5sum), ]

# Process GTDB Taxonomy
files_raw <- lapply(GTDBs.df_unique$filePath, read.delim, 
                    header = FALSE, col.names = c("SGB", "taxonomy_long"))
names(files_raw) <- GTDBs.df_unique$GTDB_filename

files_clean <- lapply(files_raw, function(x) {
  ranks <- do.call(rbind, strsplit(x$taxonomy_long, ";"))
  df    <- cbind.data.frame(ranks, SGB = x$SGB, stringsAsFactors = FALSE)
  colnames(df) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB")
  df
})

## ---- Save GTDB files ----
gtdb_destdir <- "inst/extdata/GTDB"
if (!dir.exists(gtdb_destdir)) dir.create(gtdb_destdir, recursive = TRUE)

invisible(lapply(names(files_clean), function(filename) {
  out_path <- file.path(gtdb_destdir, sub("\\.tsv$", ".tsv.bz2", filename))
  if (!grepl("\\.bz2$", out_path)) out_path <- paste0(out_path, ".bz2") # safety
  write.table(
    files_clean[[filename]],
    file      = bzfile(out_path),
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
  )
}))

write.table(
  GTDBs.df_unique[, c("GTDB_filename", "md5sum")],
  file      = file.path(gtdb_destdir, "GTDB_files_directory.tsv"),
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

## ---- Cleanup ----
unlink("tmp/", recursive = TRUE)
