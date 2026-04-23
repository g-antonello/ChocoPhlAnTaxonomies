
#' Find phylogenetically closest taxa to a target
#'
#' This function ranks SGBs in relation to phylogenetic similarity to a target
#' SGB. Important: the function does not check whether the target SGB is present
#' in a dataset: the user should be mindful of what they are searching.
#'   
#' @param SGB_target \code{character} with SGB name with or without the "t__" 
#' prefix. 
#' @param phyloTree \code{phylo} class tree. Ideally the user should consider 
#' loading trees with the package on GitHub 'g-antonello/ChocoPhlAnTaxonomies'.
#' This ensures that tip labels are identical to taxonomies' terminal ID 
#' (`Species` for chocophlan 2019, `SGB` for all later releases)
#'
#' @importFrom castor get_all_distances_to_tip
#' 
#' @returns \code{double} vector with all taxa distances sorted in increasing
#' order. The vector is named after the SGB considered. All distances are 
#' relative to the \code{SGB_target} from the \code{phyloTree} object. The 
#' self-distance, always 0, is removed from the vector. The vector returned 
#' therefore has 1 element less than the number of tips in the \code{phyloTree}.
#' 
#' @export
#'
#' @examples
#' library(ChocoPhlAnTaxonomies)
#' suppressMessages(library(TreeSummarizedExperiment))
#' phyloTree.list <- load_mpa_tree("202403")
#' dist.vec <- search_nearest_taxa_phylo("t__SGB12137", phyloTree.list$Phylo)
#' 
#' SGBs_distances <- search_nearest_taxa_phylo("t__SGB12137", phyloTree.list$Phylo)
#' head(SGBs_distances)

search_nearest_taxa_phylo <- function(SGB_target, phyloTree){
  
  if(!startsWith(SGB_target, "t__")){
    SGB_target <- paste0("t__", SGB_target)
  }

  if(!(SGB_target %in% phyloTree$tip.label)){
    stop("SGB not found in `phyloTree` tip labels")
  }

  if(!(SGB_target %in% phyloTree$tip.label)){
    stop("SGB not found in `phyloTree` tip labels")
  }
  
  # get all distances from target as a named numeric vector
  ape_coph_distances <- castor::get_all_distances_to_tip(phyloTree, SGB_target)
  ape_coph_distances <- ape_coph_distances[seq_along(phyloTree$tip.label)]
  names(ape_coph_distances) <- phyloTree$tip.label
  
  # remove the self-comparison
  ape_coph_distances_no_target <- sort(ape_coph_distances[!grepl(SGB_target, names(ape_coph_distances))])
  
  # return the vector
  return(ape_coph_distances_no_target)
}

#' Enrich distance-to-SGB vector results
#'
#' After running `search_nearest_taxa_phylo()`, the output can be passed
#' to this function to format is as complete taxonomy data.frame.
#'  
#' @param dist.vec \code{double} vector of distances. the vector's names are the 
#' SGBs.
#' @param taxonomy \code{character} vector, either `mpa` (default) or `GTDB`.
#' NB: if you are interested in eukaryotes, use only `mpa`, GTDB lacks eukaryotes
#' taxonomy
#' @param chocophlan_datestamp A datestamp string or value identifying the
#'   ChocoPhlAn database version (e.g. \code{"202307"}). Passed to
#'   \code{validate_chocophlan_version()} before use.
#' @param tse \code{TreeSummarizedExperiment}. If passed, it is used to add one
#'   extra column to the output saying whether that SGB is present in it.
#' @returns \code{data.frame} of full taxonomies and their distance to the 
#'  `SGB_target` provided in \code{search_nearest_taxa_phylo()}. if tse is provided,
#'  an extra column with presence in its rowData$SGB is provided (TRUE/FALSE).
#' 
#' @seealso \code{\link{search_nearest_taxa_phylo}}, \code{\link{search_nearest_taxa_phylo}}
#' 
#' @examples
#' library(ChocoPhlAnTaxonomies)
#' suppressMessages(library(TreeSummarizedExperiment))
#' data(JoS_2022.tse)
#' JoS_2022.tse <- AddPhyloTree_to_mpa_tse(JoS_2022.tse, chocophlan_datestamp = "202403")
#' 
#' #############################################################################
#' # example with rowTree of the tse, all will be present   ####################
#' 
#' dist.vec <- search_nearest_taxa_phylo("t__SGB12137", rowTree(JoS_2022.tse))
#' 
#' ## example without tse as argument
#' head(enrich_search_nearest_taxa_phylo(dist.vec, taxonomy = "GTDB", 
#'  chocophlan_datestamp = "202403"))
#' 
#' ## example with tse as argument
#' head(enrich_search_nearest_taxa_phylo(dist.vec, taxonomy = "GTDB", 
#' chocophlan_datestamp = "202403", tse = JoS_2022.tse))
#' 
#' #############################################################################
#' # example with full phylogenetic tree, most won't be present ################
#' phyloTree.list <- load_mpa_tree("202403")
#' dist.vec <- search_nearest_taxa_phylo("t__SGB12137", phyloTree.list$Phylo)
#' 
#' ## example without tse as argument
#' head(enrich_search_nearest_taxa_phylo(dist.vec, taxonomy = "GTDB", 
#'    chocophlan_datestamp = "202403"))
#' 
#' ## example with tse as argument
#' head(enrich_search_nearest_taxa_phylo(dist.vec, taxonomy = "GTDB", 
#'    chocophlan_datestamp = "202403", tse = JoS_2022.tse))
#' 
#' @export
enrich_search_nearest_taxa_phylo <- function(dist.vec, taxonomy = "mta", 
                                             chocophlan_datestamp, tse = NULL){
  # first validate the version
  chocophlan_datestamp <- validate_chocophlan_version(chocophlan_datestamp)
  
  # transform dist.vec into data.frame for merging
  dist.df <- data.frame(
    SGB = names(dist.vec),
    phylo_dist = as.double(dist.vec)
  )
  
  # load the taxonomy
  taxonomy.df <- load_taxonomy(taxonomy, chocophlan_datestamp)
  
  # warn about what cannot be found with the select taxonomy, but proceed
  missingSGBs <- setdiff(dist.df$SGB, taxonomy.df$SGB)
  if(length(missingSGBs) > 0) {
    warning(
      length(missingSGBs), " SGBs are in the phylogenetic tree but not in the ",
      taxonomy,
      " taxonomy (see below). If using GTDB, consider `taxonomy = 'mpa'`.:\n",
      paste(missingSGBs, collapse = "; "),
      call. = FALSE
    )
  }
  
  # merge taxonomies and distances into one data.frame
  taxonomy_dist.df <- merge(taxonomy.df, dist.df, by = "SGB")
  
  # put SGB back at the end of the data.frame but before the distances
  taxonomy_dist.df <- taxonomy_dist.df[,c(2:7,1,9)]
  
  # if the tse is provided, validate it and add column "present_in_data
  if (!is.null(tse)) {
    if (!inherits(tse, "TreeSummarizedExperiment")) {
      stop("`tse` must be a TreeSummarizedExperiment")
    }
    taxonomy_dist.df[[paste0("present_in_", deparse(substitute(tse)))]] <- 
      taxonomy_dist.df$SGB %in% rowData(tse)$SGB
  }
  
  # order rows by increasing phylogenetic distance
  taxonomy_dist.df <- taxonomy_dist.df[order(taxonomy_dist.df$phylo_dist),]
  
  return(taxonomy_dist.df)
}
