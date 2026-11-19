#' Format taxonomic names as italicised plot/table labels
#'
#' Turns raw taxonomic names, as found in GTDB or MetaPhlAn (mpa) taxonomy
#' tables, into labels where Latin names are italicised with Markdown
#' (\code{*name*}) or HTML (\code{<i>name</i>}) markup. Markdown labels can be
#' rendered in ggplot2 with \code{ggtext::element_markdown()}.
#'
#' The rules are purely structural, so no list of known names is needed:
#' \itemize{
#'   \item Rank prefixes (\code{d__}, \code{p__}, ..., \code{s__},
#'     \code{t__}) are removed.
#'   \item A Latin word is a capitalised word (\code{^[A-Z][a-z]+$}, genus or
#'     higher rank) or a lowercase word directly following a Latin genus
#'     (the specific epithet). Latin words are italicised; adjacent ones share
#'     one italic span (\code{*Agathobacter rectalis*}).
#'   \item GTDB lineage suffixes (\code{_A}, \code{_BM}, ...) are shown in
#'     parentheses after the word they belong to:
#'     \code{Blautia_A faecis} becomes \code{*Blautia* (A) *faecis*} and
#'     \code{Faecalibacterium prausnitzii_C} becomes
#'     \code{*Faecalibacterium prausnitzii* (C)}.
#'   \item Placeholder names and identifiers (anything containing a digit,
#'     a hyphen or a run of capitals, such as \code{CAG-170},
#'     \code{sp000432135}, \code{SGB4425}) stay roman and unchanged. The GTDB
#'     \code{sp<accession>} form is kept as is.
#'   \item Remaining underscores (mpa style) become spaces, and a standalone
#'     \code{sp} becomes \code{sp.}: \code{Actinomyces_sp_ICM47} becomes
#'     \code{*Actinomyces* sp. ICM47}.
#'   \item Names filled by \code{\link{fill_GTDB_unknown_taxonomies}} are
#'     rewritten from the deepest known rank: \code{g__Natrinema..s__unknown}
#'     becomes \code{*Natrinema* sp. (unknown)}, and e.g.
#'     \code{f__Natrialbaceae..g__unknown} becomes
#'     \code{*Natrialbaceae* (unknown genus)}.
#'   \item A trailing tag in square brackets (e.g. \code{" [SGB4933]"}) is
#'     kept verbatim.
#' }
#'
#' @param x A \code{character} vector (or factor) of taxonomic names, one name
#'   per element (a single rank, not a full pipe-delimited taxonomy string).
#' @param markup Either \code{"markdown"} (default, \code{*italic*}) or
#'   \code{"html"} (\code{<i>italic</i>}).
#'
#' @returns A \code{character} vector of the same length as \code{x}.
#'   \code{NA} values are returned as \code{NA}.
#'
#' @seealso \code{\link{load_taxonomy}}, \code{\link{fill_GTDB_unknown_taxonomies}}
#'
#' @export
#'
#' @examples
#' enhance_taxonomy_labels(c(
#'   "s__Agathobacter rectalis",
#'   "s__Blautia_A faecis",
#'   "s__Faecalibacterium prausnitzii_C",
#'   "s__CAG-170 sp000432135",
#'   "f__Lachnospiraceae",
#'   "g__Natrinema..s__unknown",
#'   "Dialister_invisus",
#'   "Agathobacter rectalis [SGB4933]"
#' ))
#'
#' # HTML markup, e.g. for HTML tables
#' enhance_taxonomy_labels("s__Blautia_A faecis", markup = "html")
#'
#' \dontrun{
#' # in a ggplot2 figure
#' ggplot(df, aes(x = estimate, y = taxon)) +
#'   geom_point() +
#'   scale_y_discrete(labels = enhance_taxonomy_labels) +
#'   theme(axis.text.y = ggtext::element_markdown())
#' }
enhance_taxonomy_labels <- function(x, markup = c("markdown", "html")) {
  markup <- match.arg(markup)
  italic_open  <- switch(markup, markdown = "*", html = "<i>")
  italic_close <- switch(markup, markdown = "*", html = "</i>")

  x   <- as.character(x)
  out <- x
  ok  <- !is.na(x) & nzchar(x)
  out[ok] <- vapply(x[ok], .format_one_taxon, character(1),
                    italic_open = italic_open, italic_close = italic_close,
                    USE.NAMES = FALSE)
  out
}

# ---- internal helpers --------------------------------------------------------

# rank prefix letter -> rank name, used for filled unknowns
.rank_names <- c(d = "domain", k = "kingdom", p = "phylum", c = "class",
                 o = "order", f = "family", g = "genus", s = "species")

# lowercase words that are never a specific epithet, even after a Latin genus
.non_epithet_words <- c("sp", "spp", "unclassified", "unknown", "phylum",
                        "oral", "taxon", "group", "cf", "aff", "endosymbiont")

# ...as are informal descriptors such as "bacterium", "archaeon",
# "proteobacterium", "actinobacterium" (mpa: Lachnospiraceae_bacterium)
.non_epithet_pattern <- "^[a-z]*(bacterium|archaeon)$"  # lowercase only: not Faecalibacterium

.format_one_taxon <- function(name, italic_open, italic_close) {
  name <- trimws(name)

  # 1. keep a trailing [tag] aside
  tag <- regmatches(name, regexpr("\\s*\\[[^]]*\\]$", name))
  if (length(tag)) name <- sub("\\s*\\[[^]]*\\]$", "", name)

  # 2. filled unknowns from fill_GTDB_unknown_taxonomies():
  #    "<known>..<r>__unknown", possibly chained ("o__X..f__unknown..g__unknown")
  unknown_note <- NULL
  if (grepl("\\.\\.[a-z]__unknown$", name)) {
    pieces       <- strsplit(name, "..", fixed = TRUE)[[1]]
    known        <- pieces[1]
    known_rank   <- substr(known, 1, 1)
    target_rank  <- substr(pieces[length(pieces)], 1, 1)
    unknown_note <- if (known_rank == "g" && target_rank == "s") {
      "sp. (unknown)"
    } else {
      sprintf("(unknown %s)", .rank_names[[target_rank]])
    }
    name <- known
  }

  # 3. strip rank prefixes
  name <- gsub("(^|\\s)[dkpcofgst]__", "\\1", name)

  # 4. tokenise: spaces first; then split underscores, except a GTDB suffix
  #    (_A, _BM) at the end of a token, which is kept as a separate marker
  tokens <- list()
  for (word in strsplit(name, "\\s+")[[1]]) {
    if (!nzchar(word)) next
    suffix <- NULL
    if (grepl(".+_[A-Z]{1,2}$", word)) {
      suffix <- sub(".*_([A-Z]{1,2})$", "\\1", word)
      word   <- sub("_[A-Z]{1,2}$", "", word)
    }
    for (part in strsplit(word, "_", fixed = TRUE)[[1]]) {
      if (nzchar(part)) tokens[[length(tokens) + 1]] <- list(text = part, type = "word")
    }
    if (!is.null(suffix)) tokens[[length(tokens) + 1]] <- list(text = sprintf("(%s)", suffix), type = "suffix")
  }

  # 5. decide which words are Latin (italic)
  prev_latin_capital <- FALSE
  for (i in seq_along(tokens)) {
    tok <- tokens[[i]]
    if (tok$type == "suffix") next  # a suffix does not break genus -> epithet
    is_stop    <- tolower(tok$text) %in% .non_epithet_words ||
      grepl(.non_epithet_pattern, tok$text)
    is_capital <- grepl("^[A-Z][a-z]+$", tok$text) && !is_stop
    is_epithet <- grepl("^[a-z]+$", tok$text) && !is_stop && prev_latin_capital
    tokens[[i]]$type <- if (is_capital || is_epithet) "latin" else "plain"
    if (tok$text == "sp") tokens[[i]]$text <- "sp."
    prev_latin_capital <- is_capital
  }

  # 6. assemble, merging adjacent Latin words into one italic span
  texts <- vapply(tokens, `[[`, character(1), "text")
  types <- vapply(tokens, `[[`, character(1), "type")
  runs  <- rle(types == "latin")
  ends  <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1
  pieces <- character(length(runs$lengths))
  for (r in seq_along(pieces)) {
    span <- paste(texts[starts[r]:ends[r]], collapse = " ")
    pieces[r] <- if (runs$values[r]) paste0(italic_open, span, italic_close) else span
  }
  label <- paste(pieces, collapse = " ")

  if (!is.null(unknown_note)) label <- paste(label, unknown_note)
  paste0(label, tag)
}
