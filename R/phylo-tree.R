# The first database version to use SGB-style tip labels (YYYYMM format)
SGB_TIP_LABEL_VERSION <- 202001L

#' Load MetaPhlAn phylogenetic tree
#'
#' Reads the CHOCOPhlAn phylogenetic tree for a given database version and
#' standardises tip labels. For versions >= \code{202001} (January 2020,
#' the first release to adopt SGB nomenclature), tips are prefixed with
#' \code{"t__SGB"}; for older versions the eighth pipe-delimited field is
#' extracted as the short label.
#'
#' @param chocophlan_datestamp A \code{character} string specifying the
#'   CHOCOPhlAn database version to load (e.g. \code{"202403"},
#'   \code{"202307"}, \code{"201901"}). The value is validated internally
#'   via \code{validate_chocophlan_version()} before the tree is located.
#' @param multipleTrees \code{integer} representing the wanted tree in case of
#'   multiple trees available for the same timestamp. To find out which are 
#'   available, run the function interactively or read the error message. Default 
#'   is `NULL` to prevent unexpected/unwanted trees being loaded.
#' 
#' @importFrom phytools bind.tip
#' @importFrom tidytree read.tree keep.tip
#' @importFrom ape node.depth.edgelength
#'
#' @returns A \code{list} object with two elements: `Phylo` and `Version`. The 
#' former contains the phylogenetic tree, the latter contains the name of the 
#' metaphlan tree loaded without file extensions
#'
#' @export
#'
#' @examples
#' mpaTree <- load_mpa_tree("202403")
#' mpaTree
load_mpa_tree <- function(chocophlan_datestamp, multipleTrees = NULL) {
  
  # validate chocoplan version
  chocophlan_datestamp <- validate_chocophlan_version(chocophlan_datestamp)
  
  # Find the tree
  tree_path <- grep(
    chocophlan_datestamp,
    list.files(system.file("extdata/trees/", package = "ChocoPhlAnTaxonomies"), full.names = TRUE),
    value = TRUE
  )
  
  # Fail clearly if no tree file is found for this version
  if (length(tree_path) == 0) {
    stop("No tree file found for CHOCOPhlAn version: ", chocophlan_datestamp)
  }
  
  # Fix conflicts in cases of old metaphlan versions
  if (length(tree_path) > 1) {
    formatted_trees <- paste(paste0(seq_along(tree_path), " - ", basename(tree_path)), collapse = "\n")
    
    if(is.null(multipleTrees)){
      stop(
        "Multiple tree files found for datestamp '", chocophlan_datestamp, "'.\n",
        "Choose a tree below and pass its index to the `multipleTrees` parameter:\n\n",
        formatted_trees,
        call. = FALSE
      )
    } 
    
    # validation of multipleTrees values
    if ( multipleTrees < 1 || multipleTrees > length(tree_path)){
      stop("Tree index not valid, reassign index to one of these below:\n\n", formatted_trees,  call. = FALSE)
    }
    
    # if all edge cases are satisfied, subset the tree vector
    tree_path <- tree_path[multipleTrees]    
    
    }
  
  # Read the tree directly from the connection
  tree <- read.tree(tree_path)
  
  # Versions >= 202001 use SGB nomenclature; older versions use the 8th
  # pipe-delimited field of the original tip label as the short name.
  if (as.integer(chocophlan_datestamp) >= SGB_TIP_LABEL_VERSION) {
    tree$tip.label <- paste0("t__SGB", tree$tip.label)
  } else {
    tree$tip.label <- sapply(strsplit(tree$tip.label, "|", fixed = TRUE), function(x) x[8])
  }
  
  return(list("Phylo" = tree, "Version" = gsub(".nwk.bz2", "", basename(tree_path), fixed = TRUE)))
}


#' Add a phylogenetic tree to a MetaPhlAn TreeSummarizedExperiment object
#'
#' Attaches a CHOCOPhlAn phylogenetic tree to the \code{rowTree} slot of a
#' \code{TreeSummarizedExperiment} that currently lacks one. The tree is either
#' supplied directly via \code{mpa.tre} or loaded for the requested
#' \code{chocophlan_datestamp}. An artificial \emph{UNCLASSIFIED} tip is added
#' at the root when such reads are present in the data, and the tree is then
#' pruned to the taxa actually observed in the dataset.
#'
#' @param data.tse A \code{TreeSummarizedExperiment} object whose
#'   \code{rowTree} slot is \code{NULL}. If a tree is already present the
#'   object is returned unchanged with a warning.
#' @param mpa.tre An optional pre-loaded \code{phylo} tree (e.g. from
#'   \code{\link{load_mpa_tree}}). When \code{NULL} (default) the function
#'   loads the tree internally using \code{chocophlan_datestamp}. This
#'   parameter is ignored when a tree is already present in \code{data.tse}.
#' @param chocophlan_datestamp A \code{character} string specifying the
#'   CHOCOPhlAn database version to use when \code{mpa.tre = NULL}. Supported
#'   values: \code{"latest"} (default), \code{"202503"}, \code{"202403"},
#'   \code{"202307"}. \code{"latest"} queries the current release from
#'   \url{http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/}.
#'   Ignored when \code{mpa.tre} is supplied.
#' @param multipleTrees \code{integer} representing the wanted tree in case of
#'   multiple trees available for the same timestamp. To find out which are 
#'   available, run the function interactively or read the error message. Default 
#'   is `NULL` to prevent unexpected/unwanted trees being loaded.
#'
#' @return A \code{TreeSummarizedExperiment} identical to \code{data.tse} but
#'   with the \code{rowTree} slot populated with the pruned phylogenetic tree.
#'   Row order is adjusted to match the tip order of the tree.
#'
#' @details
#' Processing steps:
#' \itemize{
#'   \item Returns \code{data.tse} unchanged (with a warning) if
#'     \code{rowTree} is already populated.
#'   \item Validates that \code{data.tse} is a
#'     \code{TreeSummarizedExperiment}.
#'   \item If \code{mpa.tre} is \code{NULL}, resolves the database version
#'     (\code{"latest"} → \code{get_mpa_latest()}; otherwise
#'     \code{validate_chocophlan_version()}) and calls \code{\link{load_mpa_tree}}.
#'   \item If exactly one row name matches \code{"UNCLASSIFIED"}, appends a
#'     tip at the root with an edge length equal to the maximum node depth
#'     of the tree.
#'   \item Prunes the tree to keep only tips present in
#'     \code{rownames(data.tse)}, erroring if there is no overlap.
#'   \item Reorders \code{data.tse} rows to match the pruned tip order.
#'   \item Assigns the pruned tree to \code{rowTree(data.tse)}.
#' }
#'
#' @importFrom tidytree read.tree keep.tip rootnode
#' @importFrom ape node.depth.edgelength
#' @importFrom phytools bind.tip
#' @importFrom TreeSummarizedExperiment rowTree rowTree<-
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom utils read.delim head
#' @importFrom methods is
#' 
#' @export
#'
#' @examples
#' suppressMessages(library(TreeSummarizedExperiment))
#' data("WallenZD_2022.tse")
#'
#' # Attach the tree from a specific database version
#' tse_with_tree <- AddPhyloTree_to_mpa_tse(
#'   WallenZD_2022.tse,
#'   chocophlan_datestamp = "201901", 
#'   # there are 2 2019 trees, the first one is annotated with
#'   # metaphlan 3.0. I will use that
#'   multipleTrees = 1 
#' )
#'
#' # Or supply a pre-loaded tree to avoid redundant I/O in loops
#' mpa_tree.list <- load_mpa_tree("201901", multipleTrees = 1)
#' mpa_tree_201901 <- mpa_tree.list$Phylo
#' 
#' tse_with_tree <- AddPhyloTree_to_mpa_tse(
#'   WallenZD_2022.tse,
#'   mpa.tre = mpa_tree_201901
#' )
AddPhyloTree_to_mpa_tse <- function(data.tse, mpa.tre = NULL, chocophlan_datestamp = "latest", multipleTrees = NULL) {
  
  # Validate input type up front for a clear error message
  if (!is(data.tse, "TreeSummarizedExperiment")) {
    stop("`data.tse` must be a TreeSummarizedExperiment object.")
  }
  
  if (!is.null(rowTree(data.tse))) {
    warning("Tree already present, returning untouched input")
    return(data.tse)
  }
  
  # Load the tree if one was not provided
  if (is.null(mpa.tre)) {
    mpa.tre.list <- load_mpa_tree(chocophlan_datestamp = chocophlan_datestamp, multipleTrees = multipleTrees)
    mpa.tre <- mpa.tre.list[["Phylo"]]
    mpa.tre.name <- mpa.tre.list[["Version"]]
  } else{
    mpa.tre.name <- deparse(substitute(mpa.tre))
  }
  
  # If there is an UNCLASSIFIED row, attach an artificial tip at the root.
  # We expect at most one such row; if multiple match, error clearly.
  unclassified_rows <- grep("UNCLASSIFIED", rownames(data.tse), value = TRUE)
  if (length(unclassified_rows) > 1) {
    stop(
      "Multiple row names match 'UNCLASSIFIED': ",
      paste(unclassified_rows, collapse = ", "),
      ". Expected at most one."
    )
  }
  if (length(unclassified_rows) == 1) {
    new_tip_length <- max(ape::node.depth.edgelength(mpa.tre))
    mpa.tre <- phytools::bind.tip(
      mpa.tre,
      where      = tidytree::rootnode(mpa.tre),
      edge.length = new_tip_length,
      tip.label  = unclassified_rows
    )
  }
  
  # Restrict the tree to the SGBs found in the dataset
  relevant_tips <- intersect(mpa.tre$tip.label, rownames(data.tse))
  if (length(relevant_tips) == 0) {
    stop(
      "No overlap between tree tips and dataset row names. ",
      "Check that the correct CHOCOPhlAn version is being used."
    )
  }
  tree_subset <- tidytree::keep.tip(mpa.tre, relevant_tips)
  
  # Reorder TSE rows to match tree tip order (good practice for downstream use)
  data_reordered.tse <- data.tse[tree_subset$tip.label, ]
  
  # Assign tree to TSE and return the enriched object
  rowTree(data_reordered.tse) <- tree_subset
  metadata(data_reordered.tse)[["PhyloTree"]] <- mpa.tre.name
  
  return(data_reordered.tse)
}