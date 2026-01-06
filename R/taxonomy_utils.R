#' Taxonomy from vector to data.frame
#' 
#' Expand a character vector of taxonomies into a full taxonomy table based on 
#' the separator provided and the prefixes at each level. Tested on metaphlan4
#' output tables
#'
#' @param taxonomy.chr A \code{character} vector with the taxonomy. They should 
#' come from a pre-cleaned taxonomy table, such as after `mia::importMetaphlan`
#' @param sep A \code{character} specifying the separator to use to split 
#' `taxonomy.chr`
#' @param row.names Either `'short'`, `'long'`, or `'none'`. Default is 
#' `'short'`, which returns the highest taxonomic resolution the input 
#' `taxonomy.chr` has for that row. `'long'` returns `taxonomy.chr` as column 
#' names. `'none'` returns a data.frame without rown ames
#
#' @returns A \code{data.frame} of taxonomies
#' @export
#'
#' @examples
#' example_taxonomies_SGB <- c(
#'   "k__Bacteria|p__Firmicutes|c__CFGB8118|o__OFGB8118|f__FGB8118|g__GGB29005|s__GGB29005_SGB41723|t__SGB41723",
#'   "k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Lachnospiraceae_unclassified|s__Lachnospiraceae_unclassified_SGB41492|t__SGB41492",
#'   "k__Bacteria|p__Bacteria_unclassified|c__Bacteria_unclassified|o__Bacteria_unclassified|f__Bacteria_unclassified|g__GGB28904|s__GGB28904_SGB41595|t__SGB41595"
#' )
#' 
#' expand_taxonomy(example_taxonomies_SGB)
#' 
#' example_taxonomies_Species <- c(
#'   "k__Bacteria|p__Firmicutes|c__CFGB8118|o__OFGB8118|f__FGB8118|g__GGB29005|s__GGB29005_SGB41723",
#'   "k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Lachnospiraceae_unclassified|s__Lachnospiraceae_unclassified_SGB41492",
#'   "k__Bacteria|p__Bacteria_unclassified|c__Bacteria_unclassified|o__Bacteria_unclassified|f__Bacteria_unclassified|g__GGB28904|s__GGB28904_SGB41595"
#' )
#' 
#' expand_taxonomy(example_taxonomies_Species)

expand_taxonomy <- function(taxonomy.chr, sep = "|", row.names = "short"){
  
  taxonomy.list <- strsplit(taxonomy.chr, sep, fixed = TRUE)
  taxonomy.df <- as.data.frame(Reduce(rbind, taxonomy.list))
  taxonomies_per_column <- apply(taxonomy.df, 2, function(x) unique(substr(x[!grepl("^UN*", x)], 1, 3)))
  
  final_colnames <- all_taxonomy_levels[taxonomies_per_column]
  
  colnames(taxonomy.df) <- final_colnames
  
  if(row.names == "short"){
    rownames(taxonomy.df) <- sapply(taxonomy.list, function(x) x[length(x)])
  } 
  
  if(row.names == "long"){
    rownames(taxonomy.df) <- taxonomy.chr
  }
  
  if(row.names == "none"){
    rownames(taxonomy.df) <- NULL
  }
  
  return(taxonomy.df)
}

#' Taxonomy from data.frame to vector
#'
#' @param taxonomy.df \code{data.frame} with taxonomic levels on columns and 
#' individual taxa on rows. Taxa levels must be ordered from Kingdom/Domain to 
#' make sense, but the function does not check that
#' @param sep \code{character} specifying the separator avoid 
#' `c("_", "-", " ", ",", "\\t")` separators, because they will likely be 
#' present in the taxonomy table already or because it would conflict with 
#' saving the taxonomy table in .csv or .tsv format
#'
#' @returns A \code{character} of complete taxonomies separated by `sep`
#' @export
#'
#' @examples
#' example_taxonomies_Species.chr <- c(
#' "k__Bacteria|p__Firmicutes|c__CFGB8118|o__OFGB8118|f__FGB8118|g__GGB29005|s__GGB29005_SGB41723",
#' "k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Lachnospiraceae_unclassified|s__Lachnospiraceae_unclassified_SGB41492",
#' "k__Bacteria|p__Bacteria_unclassified|c__Bacteria_unclassified|o__Bacteria_unclassified|f__Bacteria_unclassified|g__GGB28904|s__GGB28904_SGB41595"
#' )
#' 
#' expanded_taxonomies.df <- expand_taxonomy(example_taxonomies_Species.chr)
#' 
#' recollapsed_taxonomies.chr <- collapse_taxonomy(expanded_taxonomies.df)
#' 
#' identical(recollapsed_taxonomies.chr, example_taxonomies_Species.chr)

collapse_taxonomy <- function(taxonomy.df, sep = "|"){
  
  invalidSeps <- c("_", "-", " ", ",", "\\t")
  
  if(sep %in% invalidSeps) {
    warning(paste0("Please avoid this separator: ", sep))
  }
  
  collapsed_taxonomy <- apply(taxonomy.df, 1, function(x) paste(x, collapse = "|"))
  names(collapsed_taxonomy) <- NULL
  return(collapsed_taxonomy)
  
}

#' Enhanced row renaming in a TreeSummarizedExperiment 
#' 
#' This includes renaming not only rows of a TSE, but renames tree tips as well
#'
#' @param data.tse \code{(Tree)SummarizedExperiment} with or without a 
#' phylogenetic tree
#' @param new.rownames = A \code{character} vector of new rownames. 
#' IMPORTANT! Make sure the order of old and new row names makes sense before 
#' running this function
#'
#' @importFrom TreeSummarizedExperiment rowTree
#' 
#' @returns The same \code{TreeSummarizedExperiment}, but with reordered 
#' rownames of both `assays` and `rowData`
#' @export
#'
#' @examples
#' library(mia)
#' WallenZD_2022.tse <- mia::importMetaPhlAn(
#' file = system.file("extdata",
#'                    "WallenZD_2022_metaphlan3_profiles.tsv.bz2",
#'                    package = "biobakeryUtils"),
#' col.data = system.file("extdata", "WallenZD_2022_subjMetadata.tsv.bz2", 
#' package = "biobakeryUtils"))
#' 
#' WallenZD_2022_fullLengthNames.tse <- rename_rownames.tse(WallenZD_2022.tse, 
#' new.rownames = collapse_taxonomy(as.data.frame(rowData(WallenZD_2022.tse))))
#' 
#' head(rownames(WallenZD_2022_fullLengthNames.tse))
#' head(rownames(WallenZD_2022.tse))

rename_rownames.tse <- function(data.tse, new.rownames){
  
  rownames(data.tse) <- new.rownames
  
  # check that if there is a tree, rownames should be ordered in the same order
  # as the tip labels of the tree
  if(!is.null(rowTree(data.tse))) {
    if (!identical(rownames(data.tse), rowTree(data.tse)$tip.label)) {
      stop(
        "Your .tse object does not have identical rownames and tip names in 
        identical order. Consider running `reorder_taxa_with_phyloTree_labels` 
        first"
      )
    }
    rowTree(data.tse)$tip.label <- new.rownames
  }
  
  return(data.tse)
}

#' Reorder rows and tip labels of a TreeSummarizedExperiment
#'
#' @param data.tse \code{TreeSummarizedExperiment} object. It must contain a 
#' phylogenetic tree
#'
#' @importFrom TreeSummarizedExperiment rowTree
#' @importFrom SummarizedExperiment rowData
#'
#' @returns \code{TreeSummarizedExperiment} object with identical rownames and 
#' tip labels order
#' @export
#'
#' @examples
#' \dontrun{
#' # load a tree summarized experiment
#'  WallenZD_2022.tse <- mia::importMetaPhlAn(
#'   file = system.file("extdata",
#'                     "WallenZD_2022_metaphlan3_profiles.tsv.bz2",
#'                     package = "biobakeryUtils"),
#'   col.data = system.file("extdata", "WallenZD_2022_subjMetadata.tsv.bz2", 
#'   package = "biobakeryUtils")
#'   )
#'  # Add phylogenetic tree
#'  
#'  # shuffle rownames of the 
#'  
#'  # check if they are identical to tree tips labels
#'  
#'  tse_reordered <- reorder_taxa_with_phyloTree_labels(tse)
#' }

reorder_taxa_with_phyloTree_labels <- function(data.tse) {
  
  if(is.null(rowTree(data.tse))){
    stop("No phyloTree to reorder taxa with")
  }
  
  data_reordered.tse <- data.tse[rowTree(data.tse)$tip.label,]
  return(data_reordered.tse)
}

rename_features.tse <- function(tse, format = "long"){
  if(format == "long")
  rownames(tse) <- collapse_taxonomy(as.data.frame(rowData(tse)))
}

#' Complete Unknown taxonomy
#'  
#' This is a utility function for the `metaphlan_wrangle function`, which needs
#' it in case it finds taxa that are unique at higher taxonomies levels than the 
#' one requested (e.g. a Genus without any species mapped to it. maybe it is still 
#' a valuable genus, so it has to be repeated with the same name at the taxonomic level below).
#' This is also valuable because another approach would discard it, and the total
#' relative abundances estimated would not add up.
#' 
#' @param x \code{character},  the taxonomy as it appears in the first column 
#' of a metaphlan table
#' @param tax_lvl_int \code{character} the level at which the user wants the 
#' taxonomies aggregated. The default goes all the way to the SGB level (`8`). 
#' Species level is `7`, Genus is `6` and so on.
#'
#' @return \code{character}, the same as the input vector, with the lowest 
#' taxonomy repeated down to the `tax_lvl_int` chosen
#' @importFrom tidyr separate_wider_delim
#' @export
#'
#' @examples
#' 
#' example_taxonomies <- c(
#' "k__Bacteria|p__Firmicutes|c__CFGB1347", 
#' "k__Bacteria|p__Actinobacteria", 
#' "k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Eubacteriales_Family_XIII_Incertae_Sedis|g__Lentihominibacter|s__Lentihominibacter_faecis")
#' 
#' # complete taxonomies to species level (level of taxonomic depth)
#' complete_unknown_taxonomy(example_taxonomies, tax_lvl_int = 7)
#' 

complete_unknown_taxonomy <- function(x, tax_lvl_int = 8){
  
  #input checks 
  if (!is.character(x)) {
    stop("Input 'x' must be a character vector.")
  }
  
  if (!is.numeric(tax_lvl_int) || tax_lvl_int < 1 || tax_lvl_int > 8 || tax_lvl_int != round(tax_lvl_int)) {
    stop("'tax_lvl_int' must be an integer between 1 and 8.")
  }
  
  tmp <- sapply(strsplit(x, "\\|"), function(l) l[length(l)])
  possible_levels <- c("k__","p__","c__","o__", "f__", "g__", "s__", "t__")[1:tax_lvl_int]
  latest_taxonomy <- which(grepl(substr(tmp, 1, 3), possible_levels))
  if(length(latest_taxonomy) == 0){
    latest_taxonomy <- 0
    return(paste(rep(tmp, tax_lvl_int), collapse = "|"))
  } else{
    return(paste(c(x, rep(tmp, tax_lvl_int - latest_taxonomy)), collapse = "|")) 
  }
}

#' Upload full taxonomy data from MetaPhlAn database 
#'
#' @param version \code{character} with CHOCOPhlAn database version as found in 
#' unitn cmprod1 
#'
#' @returns \code{data.frame, data.table} object, with columns being `SGB`, 
#' `Number of alterantive columns`, `list of alternative names`
#' 
#' @importFrom data.table fread
#' @importFrom stringr str_count
#' @importFrom dplyr mutate
#' 
#' @export
#'
#' @examples
#' 
#' head(get_mpa_full_taxonomy("202403"))

get_mpa_full_taxonomy <- function(version = "202403"){
  # supported versions cached in inst/extdata
  versions <- c(
    "202403" = "mpa_vJun23_CHOCOPhlAnSGB_202403_species.txt.bz2", 
    "202503" = "mpa_vJan25_CHOCOPhlAnSGB_202503_species.txt.bz2")
  
  taxonomy <- fread(system.file("extdata", versions[version], 
                                            package = "biobakeryUtils"), 
                                header = FALSE, 
                                col.names = c("SGB", "Main_TaxName"))
  
  taxonomy_annot <- mutate(
    taxonomy,
      N_AlternateTaxNames = str_count(taxonomy$Main_TaxName, ","),
      Alternate_TaxNames = sapply(strsplit(Main_TaxName, ","), 
                                  function(x) 
                                    ifelse(length(x) == 1, 
                                           x, 
                                           paste(x[2:length(x)], collapse = ",")
                                           )
                                  )
    )

  return(taxonomy_annot)
  }

#' Find SGB alternate names
#'
#' @param SGB_id \code{character} of lenght 1 with SGB code 
#' @param version \code{character} of date in the format %Y%m that specifies
#' the CHOCOPhlAn database version
#'
#' @returns A \code{character} vector with the alternate names in CHOCOPhlAn
#' @export
#' 
#' @importFrom stringr str_remove
#' @examples
#' SGB_with_no_altNames <- get_alternate_SGB_names("t__SGB51429", version = "202403")
#' SGB_with_no_altNames
#' 
#' SGB_with_altNames <- get_alternate_SGB_names("t__SGB13547", version = "202403")
#' SGB_with_altNames

get_alternate_SGB_names <- function(SGB_id, version = "202403"){
  # remove leading t__
  SGB_id_clean <- str_remove(SGB_id, pattern =  "t__")
    
  taxonomy <- get_mpa_full_taxonomy(version)
  
  altNames <- strsplit(taxonomy$Alternate_TaxNames[taxonomy$SGB == SGB_id_clean], ",")[[1]]
  names(altNames) <- rep(SGB_id, length(altNames))
  return(altNames)
}
