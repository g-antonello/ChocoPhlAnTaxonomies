rank_cols_gtdb <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB")
rank_cols_mpa  <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB")

# ---- fill_GTDB_unknown_taxonomies -------------------------------------------

test_that("fill_GTDB_unknown_taxonomies fills empty ranks from the deepest known one", {
  tax <- data.frame(
    Domain  = c("d__Archaea", "d__Archaea", "d__Archaea"),
    Phylum  = c("p__Halobacteriota", "p__Halobacteriota", "p__Halobacteriota"),
    Class   = c("c__Halobacteria", "c__Halobacteria", "c__Halobacteria"),
    Order   = c("o__Halobacteriales", "o__Halobacteriales", "o__Halobacteriales"),
    Family  = c("f__Natrialbaceae", "f__", "f__Natrialbaceae"),
    Genus   = c("g__Natrinema", "g__", "g__Natrinema"),
    Species = c("s__", "s__", "s__Natrinema altunense")
  )
  filled <- fill_GTDB_unknown_taxonomies(tax)

  expect_equal(filled$Species[1], "g__Natrinema..s__unknown")
  expect_equal(filled$Family[2], "o__Halobacteriales..f__unknown")
  # deeper ranks chain onto the already-filled ones
  expect_equal(filled$Genus[2], "o__Halobacteriales..f__unknown..g__unknown")
  # known values untouched, no empty prefix left
  expect_equal(filled[3, ], tax[3, ])
  expect_false(any(unlist(filled) %in% c("d__", "p__", "c__", "o__", "f__", "g__", "s__")))
})

# ---- format_taxonomy / load_taxonomy ----------------------------------------

test_that("load_taxonomy returns a formatted mpa table", {
  tax <- load_taxonomy("mpa", "202307")
  expect_s3_class(tax, "data.frame")
  expect_equal(colnames(tax), rank_cols_mpa)
  expect_true(all(startsWith(tax$SGB, "t__")))
  expect_false(any(grepl("_group", tax$SGB, fixed = TRUE)))
  expect_equal(unname(unlist(tax[nrow(tax), ])),
               c(paste0(c("k", "p", "c", "o", "f", "g", "s"), "__UNCLASSIFIED"), "t__UNCLASSIFIED"))
  expect_false(anyNA(tax))
})

test_that("load_taxonomy returns a formatted GTDB table without empty ranks", {
  tax <- load_taxonomy("GTDB", "202307")
  expect_equal(colnames(tax), rank_cols_gtdb)
  expect_true(all(startsWith(tax$SGB, "t__")))
  expect_true("t__UNCLASSIFIED" %in% tax$SGB)
  expect_false(any(grepl("^[a-z]__$", unlist(tax[rank_cols_gtdb[-8]]))))
})

test_that("load_taxonomy maps 202403 to 202307 and matches format_taxonomy on the raw file", {
  expect_identical(load_taxonomy("GTDB", "202403"), load_taxonomy("GTDB", "202307"))

  raw_path <- grep("202307",
                   list.files(system.file("extdata/mpa_taxonomy", package = "ChocoPhlAnTaxonomies"),
                              full.names = TRUE),
                   value = TRUE)
  expect_identical(format_taxonomy(read.delim(raw_path)), load_taxonomy("mpa", "202307"))
})

test_that("load_taxonomy errors on bad inputs", {
  expect_error(load_taxonomy("NCBI", "202307"), "allowed parameters")
  expect_error(load_taxonomy("GTDB", "199001"), "No GTDB taxonomy file found")
})

# ---- get_SGB_full_taxonomy --------------------------------------------------

test_that("get_SGB_full_taxonomy works with and without the t__ prefix", {
  with_prefix    <- get_SGB_full_taxonomy("t__SGB4933", "mpa", "202307")
  without_prefix <- get_SGB_full_taxonomy("SGB4933", "mpa", "202307")
  expect_identical(with_prefix, without_prefix)
  expect_equal(names(with_prefix), rank_cols_mpa)
  expect_equal(unname(with_prefix["SGB"]), "t__SGB4933")

  gtdb <- get_SGB_full_taxonomy("SGB4933", "GTDB", "202403")
  expect_equal(names(gtdb), rank_cols_gtdb)
  expect_equal(unname(gtdb["Species"]), "s__Agathobacter rectalis")
})

test_that("get_SGB_full_taxonomy warns when the SGB is unknown", {
  expect_warning(res <- get_SGB_full_taxonomy("SGB0", "mpa", "202307"), "not found")
  expect_equal(res, "t__SGB0")
})

# ---- TSE helpers ------------------------------------------------------------

test_that("tse_replace_mpa_with_GTDB_taxonomy swaps rowData and keeps rows", {
  data("JoS_2022.tse", package = "ChocoPhlAnTaxonomies", envir = environment())
  gtdb.tse <- tse_replace_mpa_with_GTDB_taxonomy(JoS_2022.tse, chocophlan_datestamp = "202403")

  rd <- SummarizedExperiment::rowData(gtdb.tse)
  expect_equal(dim(gtdb.tse), dim(JoS_2022.tse))
  expect_equal(rownames(gtdb.tse), rownames(JoS_2022.tse))
  expect_equal(colnames(rd), rank_cols_gtdb)
  expect_true(all(startsWith(rd$Domain, "d__")))
  expect_equal(rd["t__SGB4933", "Species"], "s__Agathobacter rectalis")
})

test_that("tse_replace_mpa_with_GTDB_taxonomy matches _group SGBs to GTDB", {
  data("JoS_2022.tse", package = "ChocoPhlAnTaxonomies", envir = environment())
  gtdb.tse <- tse_replace_mpa_with_GTDB_taxonomy(JoS_2022.tse, chocophlan_datestamp = "202403")
  gtdb <- load_taxonomy("GTDB", "202403")
  rd <- SummarizedExperiment::rowData(gtdb.tse)

  expect_equal(rd["t__SGB8059_group", "Species"], gtdb$Species[gtdb$SGB == "t__SGB8059"])
  expect_equal(rd["t__SGB17167_group", "Species"], "s__Pauljensenia sp000278725")

  group_rows <- grep("_group$", rownames(JoS_2022.tse), value = TRUE)
  expect_gt(length(group_rows), 0)
  expect_equal(rd[group_rows, "Species"],
               gtdb$Species[match(sub("_group$", "", group_rows), gtdb$SGB)])
  # identifiers keep the suffix
  expect_equal(rd[group_rows, "SGB"], group_rows)
})

test_that("tse_replace_mpa_with_GTDB_taxonomy keeps mpa taxonomy for features missing from GTDB", {
  data("JoS_2022.tse", package = "ChocoPhlAnTaxonomies", envir = environment())
  gtdb.tse <- tse_replace_mpa_with_GTDB_taxonomy(JoS_2022.tse, chocophlan_datestamp = "202403")
  euk <- grep("^t__EUK", rownames(JoS_2022.tse), value = TRUE)
  expect_gt(length(euk), 0)
  rd_old <- SummarizedExperiment::rowData(JoS_2022.tse)[euk, ]
  rd_new <- SummarizedExperiment::rowData(gtdb.tse)[euk, ]
  expect_equal(rd_new$Species, rd_old$Species)
  expect_equal(rd_new$Domain, sub("k__", "d__", rd_old$Kingdom))
})

test_that("tse_rownames_long and tse_rownames_short round-trip", {
  data("JoS_2022.tse", package = "ChocoPhlAnTaxonomies", envir = environment())
  long_names <- tse_rownames_long(JoS_2022.tse)
  expect_length(long_names, nrow(JoS_2022.tse))
  expect_true(all(grepl("^k__.*\\|t__", long_names)))

  long.tse <- tse_rename_rownames(JoS_2022.tse, long_names)
  expect_equal(rownames(long.tse), unname(long_names))
  expect_equal(tse_rownames_short(long.tse), rownames(JoS_2022.tse))
})

test_that("tse_rename_rownames checks the length and renames tree tips too", {
  data("JoS_2022.tse", package = "ChocoPhlAnTaxonomies", envir = environment())
  expect_error(tse_rename_rownames(JoS_2022.tse, "too_short"), "has length 1")

  tree.tse  <- AddPhyloTree_to_mpa_tse(JoS_2022.tse, chocophlan_datestamp = "202403")
  new_names <- paste0("taxon_", seq_len(nrow(tree.tse)))
  expect_no_warning(renamed <- tse_rename_rownames(tree.tse, new_names))
  expect_equal(rownames(renamed), new_names)
  expect_equal(TreeSummarizedExperiment::rowTree(renamed)$tip.label, new_names)
})

test_that("tse_rename_rownames maps tips by name when rows and tips are in different order", {
  data("JoS_2022.tse", package = "ChocoPhlAnTaxonomies", envir = environment())
  tree.tse <- AddPhyloTree_to_mpa_tse(JoS_2022.tse, chocophlan_datestamp = "202403")
  shuffled <- tree.tse[rev(seq_len(nrow(tree.tse))), ]
  skip_if(identical(rownames(shuffled), TreeSummarizedExperiment::rowTree(shuffled)$tip.label),
          "subsetting reordered the tree tips as well")

  renamed <- tse_rename_rownames(shuffled, paste0("new_", rownames(shuffled)))
  expect_setequal(TreeSummarizedExperiment::rowTree(renamed)$tip.label, rownames(renamed))
  tip <- TreeSummarizedExperiment::rowTree(renamed)$tip.label[1]
  expect_equal(tip, paste0("new_", TreeSummarizedExperiment::rowTree(shuffled)$tip.label[1]))
})
