# Load Tree

#' Load mpa tree
#' 
#' The function not only reads the wanted tree, but also renames tips accordingly
#'
#' @param chocophlan_datestamp 
#'
#' @importFrom tidytree read.tree
#' 
#' @returns \code{phylotree} file loaded from extdata
#' 
#' @export
#'
#' @examples
#' 

load_mpa_tree <- function(chocophlan_datestamp){
  
  # validate the ChocoPhlAn version
  chocophlan_datestamp <- validate_mpa_version(chocophlan_datestamp)
   
  # Find the tree
  tree_path <- grep(chocophlan_datestamp, list.files(system.file("extdata/trees/", package = "ChocoPhlAnTaxonomies"), full.names = TRUE), value = TRUE)
  
  # Fix conflicts in cases of old metaphlan versions
  if(length(tree_path) > 1){
    tree_path_chice <- menu(tree_path, title = "Choose correct tree to load")
    tree_path <- tree_path[tree_path_chice]
  }
  
  # Read the tree directly from the connection
  tree <- read.tree(tree_path)
  
  # fix basic tip names to include SGB
  # pre-SGB trees shouldn't be changed, that are versions prior to 2021
  if(as.integer(chocophlan_datestamp) >= 202000){
    tree$tip.label <- paste0("t__SGB", tree$tip.label)  
  } else {
    # shorten names in remaining cases
    tree$tip.label <- sapply(strsplit(tree$tip.label, "|", fixed = TRUE), function(x) x[8])
  }
  return(tree)
}


#' This function downloads a cleaned metagenomic phylogenetic tree from the 
#' CHOCOPhlAn database and incorporates it as the `rowTree` in a 
#' `TreeSummarizedExperiment` object that currently lacks one. It supports 
#' adding a special "UNCLASSIFIED" tip to the root of the tree with branch 
#' length 1 if such reads are present in the data.
#'
#' @param data.tse A \code{TreeSummarizedExperiment} object without an existing
#' `rowTree`.
#' @param CHOCOPhlAn_version A \code{character} string specifying the version 
#' of the CHOCOPhlAn tree to download.
#' Supported versions include "latest" (default), "202503", "202403", and 
#' "202307". These correspond to dated releases available at:
#' \url{http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/}.
#'
#' @return A \code{TreeSummarizedExperiment} object identical to the input but 
#' with the `rowTree` slot populated with the downloaded and pruned phylogenetic
#' tree matching the taxa in `data.tse`.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Checks if `rowTree` is already present; if so, returns the input 
#'   unchanged.
#'   \item Downloads the specified CHOCOPhlAn tree version.
#'   \item Renames tree tips to the format "t__SGB<tip.label>" as per Biobakery
#'   community guidelines.
#'   \item Adds a tip labeled "UNCLASSIFIED" at the root if the dataset contains
#'   such reads.
#'   \item Prunes the tree to keep only tips present in the dataset.
#'   \item Reorders the `TreeSummarizedExperiment` rows to match the tree tip
#'   order.
#'   \item Assigns the pruned tree to the `rowTree` slot.
#' }
#'
#' @importFrom tidytree read.tree keep.tip
#' @importFrom ape node.depth.edgelength
#' @importFrom phytools bind.tip
#' @importFrom TreeSummarizedExperiment rowTree
#' @export
#'
#' @examples
#' \dontrun{
#' WallenZD_2022.tse <- mia::importMetaPhlAn(
#'   file = system.file("extdata", "WallenZD_2022_metaphlan3_profiles.tsv.bz2", package = "biobakeryUtils"),
#'   col.data = system.file("extdata", "WallenZD_2022_subjMetadata.tsv.bz2", package = "biobakeryUtils")
#' )
#' # Assuming tse is a TreeSummarizedExperiment without a rowTree
#' tse_with_tree <- AddPhyloTree_to_mpa_tse(WallenZD_2022.tse, CHOCOPhlAn_version = "201901")
#' }


AddPhyloTree_to_mpa_tse <- function(data.tse, CHOCOPhlAn_version = "latest") {
  if(!is.null(rowTree(data.tse))) {
    message("Tree already present, returning untouched input")
    return(data.tse)
  }
  
  if(CHOCOPhlAn_version == "latest"){
    CHOCOPhlAn_version <- get_mpa_latest()$datestamp
  }
  
  mpa.tre <- load_mpa_tree(chocophlan_datestamp = CHOCOPhlAn_version)
  
  # if there is some UNCLASSIFIED, usually the default
  if(any(grepl("UNCLASSIFIED", rownames(data.tse)))) {
    unclassified_col <- grep("UNCLASSIFIED", rownames(data.tse), value = TRUE)
    new_tip_length <- max(ape::node.depth.edgelength(mpa.tre))
    mpa.tre <- phytools::bind.tip(mpa.tre, where = tidytree::rootnode(mpa.tre), edge.length = new_tip_length, tip.label = unclassified_col)
  }
  
  relevant_tips <- intersect(mpa.tre$tip.label, rownames(data.tse))
  tree_subset <- tidytree::keep.tip(mpa.tre, relevant_tips)
  data_reordered.tse <- data.tse[tree_subset$tip.label,]
  rowTree(data_reordered.tse) <- tree_subset
  return(data_reordered.tse)
  
}




