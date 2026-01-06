
#' Find phylogenetically closest taxa to a target
#'
#' @param phyloTree \code{phylo} tree that can be handled by `ape` and 
#' `tidytree` CRAN packages. IMPORTANT: So far, this function was thought to
#'  work on full taxonomy as tip labels. So before extracting the phyloTree,
#'  please run `rename_rownames.tse(your_data.tse, collapse_taxonomy(as.data.frame(rowData(your_data.tse))))`
#' @param target_SGB \code{character} with SGB name either as "t__SGB12345" or 
#' as full taxonomy.
#' @param taxonomy_sep \code{character} specifying the separation character used
#' in the tree's tip labels. Default is "|"
#'
#' @importFrom ape cophenetic.phylo
#' @importFrom dplyr mutate arrange
#' @returns \code{data.frame} with phylogenetically closest taxa based on the 
#' input phylogenetic tree.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' WallenZD_2022.tse <- mia::importMetaPhlAn(
#'  file = system.file("extdata",
#'                    "WallenZD_2022_metaphlan3_profiles.tsv.bz2",
#'                    package = "biobakeryUtils"),
#'  col.data = system.file("extdata", "WallenZD_2022_subjMetadata.tsv.bz2", 
#'  package = "biobakeryUtils")
#'  )
#'
#' # does not work because WallenZD reported data at species level
#' WallenZD_2022.tse <- AddPhyloTree_to_mpa_tse(
#'   data.tse = WallenZD_2022.tse, 
#'   CHOCOPhlAn_version = "202403"
#'   )
#' }

search_nearest_taxa <- function(phyloTree, target_SGB, taxonomy_sep = "|"){

  # get pariwise distances
  ape_coph_distances <- cophenetic.phylo(phyloTree)
  
  # remove the self-comparison
  diag(ape_coph_distances) <- NA
  
  # get only the comparisons against the target taxon
  all_dists_vs_target <- sort(ape_coph_distances[grepl(target_SGB, rownames(ape_coph_distances)),])
  
  intermediate.df <- expand_taxonomy(names(all_dists_vs_target), sep = taxonomy_sep)
  
  intermediate.df <- mutate(intermediate.df,
    phylo_dist = all_dists_vs_target,
    known_species = !grepl("GGB|SGB|_bacterium|UNCLASSIFIED", Species) | grepl("Isolate|isolate", Species)
  ) 
  
  
  return(arrange(intermediate.df, phylo_dist))
}
