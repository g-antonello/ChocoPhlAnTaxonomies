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
#'   \code{"latest"}). Validated internally via \code{validate_chocophlan_version()}.
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
#' # These are equivalent — the "t__" prefix is stripped automatically
#' get_SGB_full_taxonomy("SGB8271", "mpa", chocophlan_datestamp = "202403")
#' get_SGB_full_taxonomy("t__SGB8271", "mpa", chocophlan_datestamp = "202403")
#'
#' # GTDB taxonomy
#' get_SGB_full_taxonomy("SGB8271", "GTDB", chocophlan_datestamp = "202403")
#' get_SGB_full_taxonomy("t__SGB8271", "GTDB", chocophlan_datestamp = "202403")
get_SGB_full_taxonomy <- function(SGB, taxonomy = "mpa", chocophlan_datestamp = "202403") {
  
  taxonomy <- match.arg(taxonomy, choices = c("mpa", "GTDB"))
  
  if (chocophlan_datestamp == "latest") {
    chocophlan_datestamp <- get_mpa_latest()$datestamp
  }
  chocophlan_datestamp <- validate_chocophlan_version(chocophlan_datestamp)
  
  # Strip "t__" prefix so lookups work regardless of input format
  SGB_lookup <- gsub("t__", "", SGB)
  
  if (taxonomy == "mpa") {
    taxonomy_path <- grep(
      chocophlan_datestamp,
      list.files(system.file("extdata/mpa_taxonomy/", package = "ChocoPhlAnTaxonomies"), full.names = TRUE),
      value = TRUE
    )
    if (length(taxonomy_path) == 0) stop("No mpa taxonomy file found for version: ", chocophlan_datestamp)
    
    taxonomy.df <- read.delim(taxonomy_path, header = FALSE)
    matched_rows <- taxonomy.df[[2]][taxonomy.df[[1]] == SGB_lookup]
    
    if (length(matched_rows) == 0) {
      warning("SGB '", SGB_lookup, "' not found in mpa taxonomy for version ", chocophlan_datestamp)
      return(NA_character_)
    }
    
    return(paste(matched_rows, paste0("t__", SGB_lookup), sep = "|"))
  }
  
  if (taxonomy == "GTDB") {
    taxonomy_path <- grep(
      chocophlan_datestamp,
      list.files(system.file("extdata/GTDB/", package = "ChocoPhlAnTaxonomies"), full.names = TRUE),
      value = TRUE
    )
    if (length(taxonomy_path) == 0) stop("No GTDB taxonomy file found for version: ", chocophlan_datestamp)
    
    taxonomy.df <- read.delim(taxonomy_path)
    SGB_lookedUp <- taxonomy.df[taxonomy.df$SGB == SGB_lookup, ]
    
    if (nrow(SGB_lookedUp) == 0) {
      warning("SGB '", SGB_lookup, "' not found in GTDB taxonomy for version ", chocophlan_datestamp)
      return(NA_character_)
    }
    
    SGB_lookedUp$SGB <- paste0("t__", SGB_lookedUp$SGB)
    return(paste(unlist(SGB_lookedUp), collapse = "|"))
  }
}


#' Replace MetaPhlAn taxonomy with GTDB taxonomy in a TreeSummarizedExperiment
#'
#' Swaps the \code{rowData} of a \code{TreeSummarizedExperiment} produced by
#' MetaPhlAn for the corresponding GTDB taxonomy from the CHOCOPhlAn database.
#' An \emph{UNCLASSIFIED} row is added automatically if present in the TSE.
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
  # if a chocophlan datestamp is not provided, try to guess it
  if(is.null(chocophlan_datestamp)){
    if("db_version" %in% colnames(colData(data.tse))){
      chocophlan_datestamp <- sapply(strsplit(unique(colData(data.tse)$db_version), "_"), function(x) x[length(x)])
    } else{
      stop("db_version column not found, please provide a value in 
      chocophlan_datestamp (eg: 202401")
    }
  }
  
  chocophlan_datestamp <- validate_chocophlan_version(chocophlan_datestamp)
  
  taxonomy_path <- grep(
    chocophlan_datestamp,
    list.files(system.file("extdata/GTDB/", package = "ChocoPhlAnTaxonomies"), full.names = TRUE),
    value = TRUE
  )
  if (length(taxonomy_path) == 0) stop("No GTDB taxonomy file found for version: ", chocophlan_datestamp)
  
  taxonomy.df <- read.delim(taxonomy_path)
  taxonomy.df$SGB <- paste0("t__", taxonomy.df$SGB)
  
  # keep missing taxa from origin to ensure all are there in target taxonomy
  missing_taxa_in_GTDB.df <- as.data.frame(rowData(data.tse)[rev(setdiff(rownames(data.tse), taxonomy.df$SGB)),])
  # manipulate Kingdom into Domain in 2 steps
  colnames(missing_taxa_in_GTDB.df) <- colnames(taxonomy.df)
  missing_taxa_in_GTDB.df$Domain <- gsub("k__", "d__", missing_taxa_in_GTDB.df$Domain)
  
  # bind new taxonomy with leftovers of the old one (UNCLASSIFIED and Eukarya)
  taxonomy.df <- rbind.data.frame(taxonomy.df, missing_taxa_in_GTDB.df)
  
  rowData(data.tse) <- DataFrame(taxonomy.df[rownames(data.tse), ])
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
      stop(
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