# ---- load_mpa_tree ----------------------------------------------------------

test_that("load_mpa_tree loads SGB trees with t__SGB tip labels", {
  tree.list <- load_mpa_tree("202307")
  expect_named(tree.list, c("Phylo", "Version"))
  expect_s3_class(tree.list$Phylo, "phylo")
  expect_equal(tree.list$Version, "mpa_vJun23_CHOCOPhlAnSGB_202307")
  expect_true(all(startsWith(tree.list$Phylo$tip.label, "t__SGB")))
  expect_equal(ape::Ntip(tree.list$Phylo), 36273)
})

test_that("load_mpa_tree maps 202403 to the 202307 tree", {
  expect_equal(load_mpa_tree("202403")$Version, "mpa_vJun23_CHOCOPhlAnSGB_202307")
})

test_that("load_mpa_tree asks to choose among multiple trees for old versions", {
  expect_error(load_mpa_tree("201901"), "Multiple tree files found")
  expect_error(load_mpa_tree("201901", multipleTrees = 3), "Tree index not valid")

  v31 <- load_mpa_tree("201901", multipleTrees = 2)
  expect_equal(v31$Version, "mpa_v31_CHOCOPhlAn_201901_species_tree")
  # pre-SGB trees use the species field of the original tip label
  expect_true(all(startsWith(v31$Phylo$tip.label, "s__")))
})

test_that("load_mpa_tree errors when no tree exists for a version", {
  expect_error(load_mpa_tree("199001"), "No tree file found")
})

# ---- AddPhyloTree_to_mpa_tse ------------------------------------------------

test_that("AddPhyloTree_to_mpa_tse attaches a tree matching the rows", {
  data("JoS_2022.tse", package = "ChocoPhlAnTaxonomies", envir = environment())
  tree.tse <- AddPhyloTree_to_mpa_tse(JoS_2022.tse, chocophlan_datestamp = "202403")
  tree <- TreeSummarizedExperiment::rowTree(tree.tse)

  expect_s4_class(tree.tse, "TreeSummarizedExperiment")
  expect_equal(rownames(tree.tse), tree$tip.label)
  expect_true(all(rownames(tree.tse) %in% rownames(JoS_2022.tse)))
  expect_equal(S4Vectors::metadata(tree.tse)$PhyloTree, "mpa_vJun23_CHOCOPhlAnSGB_202307")
  # UNCLASSIFIED is kept, attached at the root
  expect_true("t__UNCLASSIFIED" %in% rownames(tree.tse))
})

test_that("AddPhyloTree_to_mpa_tse accepts a pre-loaded tree", {
  data("JoS_2022.tse", package = "ChocoPhlAnTaxonomies", envir = environment())
  my_tree <- load_mpa_tree("202307")$Phylo
  tree.tse <- AddPhyloTree_to_mpa_tse(JoS_2022.tse, mpa.tre = my_tree)
  expect_equal(S4Vectors::metadata(tree.tse)$PhyloTree, "my_tree")
  expect_equal(rownames(tree.tse), TreeSummarizedExperiment::rowTree(tree.tse)$tip.label)
})

test_that("AddPhyloTree_to_mpa_tse validates its input", {
  data("JoS_2022.tse", package = "ChocoPhlAnTaxonomies", envir = environment())
  expect_error(AddPhyloTree_to_mpa_tse(data.frame()), "must be a TreeSummarizedExperiment")

  tree.tse <- AddPhyloTree_to_mpa_tse(JoS_2022.tse, chocophlan_datestamp = "202403")
  expect_warning(again <- AddPhyloTree_to_mpa_tse(tree.tse, chocophlan_datestamp = "202403"),
                 "Tree already present")
  expect_identical(again, tree.tse)

  # mpa3-style rownames (s__...) do not overlap with an SGB tree
  data("WallenZD_2022.tse", package = "ChocoPhlAnTaxonomies", envir = environment())
  expect_error(AddPhyloTree_to_mpa_tse(WallenZD_2022.tse, chocophlan_datestamp = "202307"),
               "No overlap")
})

test_that("AddPhyloTree_to_mpa_tse keeps _group SGBs and drops only features absent from the tree", {
  data("JoS_2022.tse", package = "ChocoPhlAnTaxonomies", envir = environment())
  tree.tse <- AddPhyloTree_to_mpa_tse(JoS_2022.tse, chocophlan_datestamp = "202403")
  tips <- TreeSummarizedExperiment::rowTree(tree.tse)$tip.label

  group_rows <- grep("_group$", rownames(JoS_2022.tse), value = TRUE)
  expect_gt(length(group_rows), 0)
  expect_true(all(group_rows %in% rownames(tree.tse)))
  expect_true(all(group_rows %in% tips))
  expect_equal(rownames(tree.tse), tips)

  # the only rows dropped are eukaryotes, which are not in the SGB tree
  dropped <- setdiff(rownames(JoS_2022.tse), rownames(tree.tse))
  expect_true(all(startsWith(dropped, "t__EUK")))
})
