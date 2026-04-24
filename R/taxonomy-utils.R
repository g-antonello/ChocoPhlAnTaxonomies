#' Fill Unknown Taxonomic Ranks in GTDB-style Taxonomy Table
#'
#' For each rank column that contains empty entries (e.g. \code{"g__"}),
#' fills in a derived name based on the deepest known rank above it, combined
#' with the current rank prefix and the suffix \code{"unknown"}. Rows that are
#' genuinely indistinguishable (identical known taxonomy above the missing rank)
#' will receive identical filled values, which is scientifically correct.
#'
#' @param tax.df A \code{data.frame} with 7 columns in GTDB rank order:
#'   Domain, Phylum, Class, Order, Family, Genus, Species. Each cell should
#'   contain the rank prefix followed by the name (e.g. \code{"g__Natrinema"}),
#'   or just the bare prefix (e.g. \code{"g__"}) when the rank is unknown.
#'
#' @return The input \code{data.frame} with empty rank cells replaced by
#'   derived names of the form \code{"<last_known_name>_<prefix>unknown"},
#'   e.g. \code{"g__Natrinema_s__unknown"}.
#'
#' @examples
#' tax <- data.frame(
#'   Domain = c("d__Archaea", "d__Archaea"),
#'   Phylum = c("p__Halobacteriota", "p__Halobacteriota"),
#'   Class  = c("c__Halobacteria", "c__Halobacteria"),
#'   Order  = c("o__Halobacteriales", "o__Halobacteriales"),
#'   Family = c("f__Natrialbaceae", "f__"),
#'   Genus  = c("g__Natrinema", "g__"),
#'   Species = c("s__", "s__")
#' )
#' 
#' fill_GTDB_unknown_taxonomies(tax)
#'
#' @export
fill_GTDB_unknown_taxonomies <- function(tax.df) {
  prefixes <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")
  empty    <- paste0(prefixes, "")
  
  for (col in seq_along(prefixes)) {
    is_empty <- tax.df[[col]] == empty[col]
    if (!any(is_empty)) next
    
    last_known <- apply(tax.df[is_empty, seq_len(col - 1), drop = FALSE], 1, function(row) {
      row[max(which(row != empty[seq_len(col - 1)]))]
    })
    
    tax.df[is_empty, col] <- paste0(last_known, "..", prefixes[col], "unknown")
  }
  
  tax.df
}

#' Format a Raw Taxonomy Table into a Standardised Wide Data Frame
#'
#' Accepts taxonomy tables in either MetaPhlAn (mpa) or GTDB format and
#' returns a consistently structured \code{data.frame} with one column per
#' taxonomic rank plus an \code{SGB} column.
#'
#' For mpa input (\code{ncol == 2}): the second column is expected to contain
#' pipe-delimited taxonomy strings, optionally comma-separated when multiple
#' names exist for one SGB (only the first is used).
#'
#' For GTDB input (\code{ncol > 2}): unknown ranks are filled via
#' \code{\link{fill_GTDB_unknown_taxonomies}}.
#'
#' @param taxonomy.df A \code{data.frame}. For mpa format: 2 columns
#'   (SGB, taxonomy string). For GTDB format: 8 columns (SGB + 7 rank columns
#'   in order Domain through Species).
#'
#' @return A \code{data.frame} with columns: \code{Domain} (or \code{Kingdom}
#'   for mpa), \code{Phylum}, \code{Class}, \code{Order}, \code{Family},
#'   \code{Genus}, \code{Species}, \code{SGB}.
#'
#' @seealso \code{\link{fill_GTDB_unknown_taxonomies}}, \code{\link{load_taxonomy}}
#' 
#' @examples
#' chocophlan_datestamp <- 202307
#' taxonomy <- "mpa"
#' 
#' # load taxonomy raw
#' taxonomy_path <- grep(
#'  chocophlan_datestamp,
#'  list.files(system.file(sprintf("extdata/%s_taxonomy/", taxonomy), 
#'      package = "ChocoPhlAnTaxonomies"), 
#'    full.names = TRUE), value = TRUE)
#'    
#' taxonomy_raw <- read.delim(taxonomy_path)
#' taxonomy_formatted_manually <- format_taxonomy(taxonomy_raw)
#' taxonomy_cleaned_internally <- load_taxonomy(taxonomy, chocophlan_datestamp)
#' # the content is identical
#' identical(taxonomy_formatted_manually, taxonomy_cleaned_internally)
#' 
#' @export
format_taxonomy <- function(taxonomy.df){
  # Case 1 - taxonomy is mpa
  if(ncol(taxonomy.df) == 2){
    
    # The column name is itself a taxonomy string (mpa artifact) — recover it
    all_tax_strings <- c(
      gsub(".", "|", colnames(taxonomy.df)[2], fixed = TRUE),
      as.character(taxonomy.df[[2]])
    )
    all_sgb <- c(
      gsub(".", "|", colnames(taxonomy.df)[1], fixed = TRUE),
      as.character(taxonomy.df[[1]])
    )
    
    # mpa sometimes has multiple names for the same taxonomy. However, in 
    # practice, the first name is returned in the output
    taxonomies_first_occurrence <- sapply(strsplit(all_tax_strings, ",", fixed = TRUE), `[`, 1)
    
    # expand taxonomies with strsplit 
    taxonomies_first_occurrence_expanded <- do.call(rbind, strsplit(taxonomies_first_occurrence, "|", fixed = TRUE))
    taxonomies_first_occurrence_expanded <- as.data.frame(taxonomies_first_occurrence_expanded)
    
    # bind the SGB column at the end
    taxonomy_final.df <- cbind.data.frame(taxonomies_first_occurrence_expanded, all_sgb)
    
    # rename columns with proper taxonomy
    colnames(taxonomy_final.df) <- c("Kingdom",
                                     "Phylum",
                                     "Class",
                                     "Order",
                                     "Family",
                                     "Genus",
                                     "Species",
                                     "SGB")
  }
  #-----      End of Case 1                               ----
  ############################################################
  
  ############################################################
  # Case 2 - taxonomy is GTDB
  
  if(ncol(taxonomy.df) > 2){
    # fill unknown taxonomies in GTDB, separate function as it is rather long
    taxonomy_final.df <- fill_GTDB_unknown_taxonomies(taxonomy.df)  
  }
  
  # add t__ to SGB column
  taxonomy_final.df$SGB <- paste0("t__", taxonomy_final.df$SGB)
  
  # add UNCLASSIFIED row but taking the prefixes of each tax. level
  # (here i take the first letter of each taxonomic rank and make into a taxonomy)
  unclassified.df <- data.frame(t(paste0(tolower(substr(colnames(taxonomy_final.df), 1, 1)), "__", "UNCLASSIFIED")))
  # fix 8th prefix
  unclassified.df[1,8] <- gsub("s__", "t__", unclassified.df[1,8])
  colnames(unclassified.df) <- colnames(taxonomy_final.df)
  
  taxonomy_final.df <- rbind.data.frame(taxonomy_final.df, unclassified.df)
  
  # remove _group from SGB. it is anyway present in the original files
  taxonomy_final.df$SGB <- gsub("_group", "", taxonomy_final.df$SGB, fixed = TRUE)
  
  # return formatted object
  return(taxonomy_final.df)
}

#' Load and Format a ChocoPhlAn Taxonomy Table
#'
#' Locates the taxonomy file for a given ChocoPhlAn datestamp and taxonomy
#' type, reads it, and returns a formatted wide \code{data.frame} via
#' \code{\link{format_taxonomy}}.
#'
#' @param taxonomy Character string specifying the taxonomy type. Must be
#'   one of \code{"mpa"} or \code{"GTDB"}. Defaults to \code{"mpa"}.
#' @param chocophlan_datestamp A datestamp string or value identifying the
#'   ChocoPhlAn database version (e.g. \code{"202307"}). Passed to
#'   \code{validate_chocophlan_version()} before use.
#'
#' @return A formatted \code{data.frame} with one column per taxonomic rank
#'   and an \code{SGB} column. See \code{\link{format_taxonomy}} for details.
#'
#' @examples
#' tax <- load_taxonomy("mpa", chocophlan_datestamp = "202307")
#' head(tax)
#'
#' @seealso \code{\link{format_taxonomy}}, \code{\link{fill_GTDB_unknown_taxonomies}}
#'
#' @export
load_taxonomy <- function(taxonomy = "mpa", chocophlan_datestamp){
  # validate parameter input
  if(!(taxonomy %in% c("mpa", "GTDB"))) stop("allowed parameters for taxonomy are: 'mpa' (MetaPhlAn) and 'GTDB' (Genome Taxonomy DataBase)")

  chocophlan_datestamp <- validate_chocophlan_version(chocophlan_datestamp)
  
  taxonomy_path <- grep(
    chocophlan_datestamp,
    list.files(system.file(sprintf("extdata/%s_taxonomy/", taxonomy), package = "ChocoPhlAnTaxonomies"), full.names = TRUE),
    value = TRUE
  )
    
  if (length(taxonomy_path) == 0) stop("No ", taxonomy, " taxonomy file found for version: ", chocophlan_datestamp)
    
  # once taxonomy path is found, load it only once
  taxonomy.df <- read.delim(taxonomy_path)
  
  taxonomy_formatted.df <- format_taxonomy(taxonomy.df)
  return(taxonomy_formatted.df)
}

#' Get full taxonomy for a SGB
#'
#' Given a SGB identifier, returns its full pipe-delimited taxonomy string
#' from either the standard MetaPhlAn (\code{"mpa"}) or GTDB taxonomy for
#' the requested CHOCOPhlAn database version. The \code{"t__"} prefix is
#' accepted but not required in \code{SGB}.
#'
#' @param SGB A \code{character} string identifying the SGB to look up,
#'   with or without a \code{"t__"} prefix (e.g. \code{"SGB8271"} or
#'   \code{"t__SGB8271"}).
#' @param taxonomy A \code{character} string specifying which taxonomy to
#'   use. Must be one of \code{"mpa"} (default) or \code{"GTDB"}.
#' @param chocophlan_datestamp A \code{character} string specifying the
#'   CHOCOPhlAn database version (e.g. \code{"202403"}, \code{"202307"},
#'   Default: \code{"latest"}). Validated internally via \code{validate_chocophlan_version()}.
#'
#' @returns A \code{character} string containing the full pipe-delimited
#'   taxonomy for the requested SGB, with the terminal field in the form
#'   \code{"t__SGB<number>"}.  Returns \code{NA_character_} with a warning
#'   if the SGB is not found in the taxonomy table.
#' 
#' @importFrom utils read.delim
#' @export
#'
#' @examples
#' # These are equivalent — the "t__" prefix is added automatically
#' get_SGB_full_taxonomy("SGB8271", "mpa", chocophlan_datestamp = "202403")
#' get_SGB_full_taxonomy("t__SGB8271", "mpa", chocophlan_datestamp = "202403")
#'
#' # GTDB taxonomy
#' get_SGB_full_taxonomy("SGB8271", "GTDB", chocophlan_datestamp = "202403")
#' get_SGB_full_taxonomy("t__SGB8271", "GTDB", chocophlan_datestamp = "202403")
get_SGB_full_taxonomy <- function(SGB, taxonomy = "mpa", chocophlan_datestamp = "latest") {
  
  chocophlan_datestamp <- validate_chocophlan_version(chocophlan_datestamp)
  
  SGB_lookup <- SGB
  # make sure t__ is in front of SGB
  if(!startsWith(SGB, "t__")){
    SGB_lookup <- paste0("t__", SGB_lookup)
    }
  
  # load the requested taxonomy 
  taxonomy.df <- load_taxonomy(taxonomy, chocophlan_datestamp)
  
  # find the taxonomy row that matches the SGB
  matched_rows <- which(taxonomy.df$SGB == SGB_lookup)
  
  if (length(matched_rows) == 0) {
    warning("SGB '", SGB_lookup, "' not found in ", taxonomy, " taxonomy for version ", chocophlan_datestamp)
    return(SGB_lookup)
  } else {
    SGB_lookedUp <- unlist(taxonomy.df[matched_rows,])
    return(SGB_lookedUp)
  }
  }

#' Replace MetaPhlAn taxonomy with GTDB taxonomy in a TreeSummarizedExperiment
#'
#' Swaps the \code{rowData} of a \code{TreeSummarizedExperiment} produced by
#' MetaPhlAn for the corresponding GTDB taxonomy from the CHOCOPhlAn database.
#' Missing taxonomiese from GTDB are incorporated as they are from `mpa`. This
#' applies only to eukaryotes ("t__EUK...")
#'
#' @param data.tse A \code{TreeSummarizedExperiment} whose \code{rowData} uses
#'   MetaPhlAn taxonomy and whose row names are SGB identifiers in the form
#'   \code{"t__SGB<number>"} or \code{"UNCLASSIFIED"}.
#' @param chocophlan_datestamp A \code{character} string specifying the
#'   CHOCOPhlAn database version to use (e.g. \code{"202403"}, \code{"202307"}).
#'   Validated internally via \code{validate_chocophlan_version()}.
#'
#' @returns The input \code{TreeSummarizedExperiment} with \code{rowData}
#'   replaced by the GTDB taxonomy table, subset and ordered to match the
#'   rows of \code{tse}.
#' @importFrom SummarizedExperiment colData rowData rowData<-
#' @importFrom S4Vectors DataFrame
#' @importFrom utils read.delim head
#' @export
#'
#' @examples
#' suppressMessages(library(TreeSummarizedExperiment))
#' data("JoS_2022.tse")
#' 
#' tse_gtdb <- tse_replace_mpa_with_GTDB_taxonomy(JoS_2022.tse, chocophlan_datestamp = "202403")
#' head(rowData(tse_gtdb))
tse_replace_mpa_with_GTDB_taxonomy <- function(data.tse, chocophlan_datestamp = NULL) {
  
  # validate chocophlan datestamp
  chocophlan_datestamp <- validate_chocophlan_version(chocophlan_datestamp)
  
  # load taxonomy
  taxonomy.df <- load_taxonomy("GTDB", chocophlan_datestamp)
  
  # make sure the new taxonomy table has all features that the old one had.
  # if not, add them as they are from the old taxonomy
  missing_taxa_in_GTDB.df <- as.data.frame(rowData(data.tse)[rev(setdiff(rownames(data.tse), taxonomy.df$SGB)),])
  # manipulate Kingdom into Domain in 2 steps
  colnames(missing_taxa_in_GTDB.df) <- colnames(taxonomy.df)
  missing_taxa_in_GTDB.df$Domain <- gsub("k__", "d__", missing_taxa_in_GTDB.df$Domain)
  
  # bind new taxonomy with leftovers of the old one (typically Eukarya)
  taxonomy_new.df <- rbind.data.frame(taxonomy.df, missing_taxa_in_GTDB.df)
  # put rownames in the new taxonomy
  rownames(taxonomy_new.df) <- taxonomy_new.df$SGB
  
  # check that all rownames are there, produce error if not
  stopifnot(all(rownames(data.tse) %in% rownames(taxonomy_new.df)))
  
  # replace the old taxonomy with the new one 
  rowData(data.tse) <- DataFrame(taxonomy_new.df[rownames(data.tse), ])
  return(data.tse)
}


#' Rename rows (and tree tips) of a TreeSummarizedExperiment
#'
#' Renames the rows of a \code{(Tree)SummarizedExperiment}, and if a
#' \code{rowTree} is present, renames the matching tree tip labels in the
#' same operation. The new and old row names must be in the same order.
#'
#' @param data.tse A \code{(Tree)SummarizedExperiment} with or without a
#'   phylogenetic tree.
#' @param new.rownames A \code{character} vector of new row names, the same
#'   length as \code{nrow(data.tse)} and in the same order as the current
#'   row names.
#'
#' @importFrom TreeSummarizedExperiment rowTree
#'
#' @returns The input \code{(Tree)SummarizedExperiment} with updated row
#'   names in both the \code{assays}/\code{rowData} slots and, if present,
#'   the \code{rowTree} tip labels.
#'
#' @export
#'
#' @examples
#' library(TreeSummarizedExperiment)
#' data("WallenZD_2022.tse", package = "ChocoPhlAnTaxonomies")
#'
#' WallenZD_2022_fullNames.tse <- tse_rename_rownames(
#'   WallenZD_2022.tse,
#'   new.rownames = tse_rownames_long(WallenZD_2022.tse))
#'
#' head(rownames(WallenZD_2022_fullNames.tse))
#' head(rownames(WallenZD_2022.tse))
tse_rename_rownames <- function(data.tse, new.rownames) {
  
  if (length(new.rownames) != nrow(data.tse)) {
    stop(
      "`new.rownames` has length ", length(new.rownames),
      " but `data.tse` has ", nrow(data.tse), " rows."
    )
  }
  
  rownames(data.tse) <- new.rownames
  
  if (!is.null(rowTree(data.tse))) {
    if (!identical(rownames(data.tse), rowTree(data.tse)$tip.label)) {
      warning(
        "Row names and tree tip labels are not in identical order. ",
        "Consider running `reorder_taxa_with_phyloTree_labels` first."
      )
    }
    rowTree(data.tse)$tip.label <- new.rownames
  }
  
  return(data.tse)
}


#' Expand TSE row names to full taxonomy strings
#'
#' Replaces short terminal-node row names (e.g. \code{"t__SGB8271"}) with
#' full pipe-delimited taxonomy strings built by collapsing all \code{rowData}
#' columns. Useful when downstream functions require sorting by full taxonomy,
#' or when full taxonomy labels are wanted on a phylogenetic tree.
#'
#' @param tse A \code{TreeSummarizedExperiment} whose \code{rowData} columns
#'   contain the taxonomic ranks to be concatenated.
#'
#' @returns  \code{character} of full-taxonomy rownames ready to be reassigned 
#' to the TSE using `tse_rename_rownames()`. IMPORTANT: There are special characters
#' in the rownames, modeling them will likely silently convert them with 
#' make.names().
#'
#' @importFrom SummarizedExperiment rowData
#' @export
#'
#' @examples
#' data("WallenZD_2022.tse")
#' tse_long <- tse_rownames_long(WallenZD_2022.tse)
#' 
#' # inspect old vs new rownames
#' head(rownames(WallenZD_2022.tse))
#' head(rownames(tse_long))
tse_rownames_long <- function(tse) {
  new_rownames <- apply(rowData(tse), 1, function(x) paste(x, collapse = "|"))
  return(new_rownames)
}


#' Shorten TSE row names to terminal taxonomy identifiers
#'
#' Replaces full pipe-delimited taxonomy row names with only the last
#' pipe-delimited field (e.g. \code{"t__SGB8271"}). This is the inverse of
#' \code{\link{tse_rownames_long}}.
#'
#' @param tse A \code{TreeSummarizedExperiment} whose row names are full
#'   pipe-delimited taxonomy strings.
#'
#' @returns  \code{character} of terminal-only taxonomy rownames ready to be 
#' reassigned to the TSE using `tse_rename_rownames()`. 
#'
#' @export 
#'
#' @examples
#' \dontrun{
#' tse_short <- tse_rownames_short(tse)
#' head(rownames(tse_short))
#' }
tse_rownames_short <- function(tse) {
  new_rownames <- sapply(strsplit(rownames(tse), "|", fixed = TRUE), function(x) x[length(x)])
  return(new_rownames)
}
