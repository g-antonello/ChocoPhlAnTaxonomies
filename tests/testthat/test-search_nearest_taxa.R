# small tree whose tip names are prefixes of each other (t__SGB19 vs t__SGB1914)
toy_tree <- ape::read.tree(
  text = "((t__SGB19:1,t__SGB1914:2):1,(t__SGB4933:1,t__SGB4934:3):2);"
)

# ---- search_nearest_taxa_phylo ----------------------------------------------

test_that("search_nearest_taxa_phylo returns sorted distances to all other tips", {
  d <- search_nearest_taxa_phylo("t__SGB19", toy_tree)
  expect_equal(d, c(t__SGB1914 = 3, t__SGB4933 = 5, t__SGB4934 = 7))
})

test_that("search_nearest_taxa_phylo only removes the target itself (exact match)", {
  # regression: grepl("t__SGB19", ...) used to also drop t__SGB1914
  d <- search_nearest_taxa_phylo("SGB19", toy_tree)
  expect_length(d, ape::Ntip(toy_tree) - 1)
  expect_true("t__SGB1914" %in% names(d))
  expect_false("t__SGB19" %in% names(d))
})

test_that("search_nearest_taxa_phylo accepts SGBs without t__ and errors on unknown ones", {
  expect_identical(search_nearest_taxa_phylo("SGB4933", toy_tree),
                   search_nearest_taxa_phylo("t__SGB4933", toy_tree))
  expect_error(search_nearest_taxa_phylo("SGB0", toy_tree), "SGB not found")
})

# ---- enrich_search_nearest_taxa_phylo ---------------------------------------

test_that("enrich_search_nearest_taxa_phylo adds taxonomy, sorted by distance", {
  d <- search_nearest_taxa_phylo("t__SGB19", toy_tree)
  # default taxonomy is "mpa"
  res <- suppressWarnings(enrich_search_nearest_taxa_phylo(d, chocophlan_datestamp = "202307"))
  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res),
               c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB", "phylo_dist"))
  expect_false(is.unsorted(res$phylo_dist))
  expect_equal(res$SGB[res$SGB == "t__SGB4933"], "t__SGB4933")
})

test_that("enrich_search_nearest_taxa_phylo warns about SGBs missing from the taxonomy", {
  d <- c(t__SGB4933 = 1, t__SGB0 = 2)
  expect_warning(res <- enrich_search_nearest_taxa_phylo(d, taxonomy = "GTDB", chocophlan_datestamp = "202403"),
                 "1 SGBs are in the phylogenetic tree but not in the GTDB taxonomy")
  expect_equal(res$SGB, "t__SGB4933")
})

test_that("enrich_search_nearest_taxa_phylo flags presence in a TSE", {
  data("JoS_2022.tse", package = "ChocoPhlAnTaxonomies", envir = environment())
  present <- "t__SGB15291"
  absent  <- setdiff(load_taxonomy("mpa", "202307")$SGB,
                     SummarizedExperiment::rowData(JoS_2022.tse)$SGB)[1]
  d <- setNames(c(1, 2), c(present, absent))
  res <- enrich_search_nearest_taxa_phylo(d, taxonomy = "mpa", chocophlan_datestamp = "202307",
                                          tse = JoS_2022.tse)
  expect_equal(res$SGB, c(present, absent))
  expect_equal(res$present_in_JoS_2022.tse, c(TRUE, FALSE))
  expect_error(enrich_search_nearest_taxa_phylo(d, "mpa", "202307", tse = data.frame()),
               "must be a TreeSummarizedExperiment")
})
