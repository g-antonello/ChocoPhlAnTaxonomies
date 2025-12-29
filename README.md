# CHOCOPhlAn tree abd Genome Taxonomy Database files

This repository aims to collect CHOCOPhlAn versions' trees and genome taxonomy database files 
with some manual curation in case it is needed. 
In brief, the procedure is the following:

  - Download all downloadable MetaPhlAn versions from Zenodo
  - Search for files with the `.nwk` extension (phylogentic tree)
  - Annotate which version of MetaPhlAn they are found in into a .tsv file that goes into `inst/extdata/trees`\
  - rename their tips as SGB[...] with [...] being the original tip label
  - Copy them into inst/extdata
  - Repeat the same with files that have the `GTDB` pattern in its name
  - GTDB files are not just copied, but also formatted as a more typical taxonomy table

