test_that("plain Latin binomials and higher ranks are italicised", {
  expect_equal(enhance_taxonomy_labels("s__Agathobacter rectalis"), "*Agathobacter rectalis*")
  expect_equal(enhance_taxonomy_labels("Agathobacter rectalis"), "*Agathobacter rectalis*")
  expect_equal(enhance_taxonomy_labels("f__Lachnospiraceae"), "*Lachnospiraceae*")
  expect_equal(enhance_taxonomy_labels("k__Bacteria"), "*Bacteria*")
})

test_that("GTDB suffixes go in parentheses after the rank they belong to", {
  expect_equal(enhance_taxonomy_labels("s__Blautia_A faecis"), "*Blautia* (A) *faecis*")
  expect_equal(enhance_taxonomy_labels("Faecalibacterium prausnitzii_C"),
               "*Faecalibacterium prausnitzii* (C)")
  expect_equal(enhance_taxonomy_labels("Clostridium_Q saccharolyticum_A"),
               "*Clostridium* (Q) *saccharolyticum* (A)")
  expect_equal(enhance_taxonomy_labels("Streptococcus mitis_BM"), "*Streptococcus mitis* (BM)")
  expect_equal(enhance_taxonomy_labels("p__Firmicutes_A"), "*Firmicutes* (A)")
  expect_equal(enhance_taxonomy_labels("f__2-02-FULL-58-16_A"), "2-02-FULL-58-16 (A)")
})

test_that("GTDB placeholders and sp<accession> are kept as they are", {
  expect_equal(enhance_taxonomy_labels("s__CAG-170 sp000432135"), "CAG-170 sp000432135")
  expect_equal(enhance_taxonomy_labels("1XD42-69 sp003477885"), "1XD42-69 sp003477885")
  expect_equal(enhance_taxonomy_labels("An92 sp900199495"), "An92 sp900199495")
  expect_equal(enhance_taxonomy_labels("Blautia sp000432135"), "*Blautia* sp000432135")
  expect_equal(enhance_taxonomy_labels("Blautia_A sp900120195"), "*Blautia* (A) sp900120195")
  expect_equal(enhance_taxonomy_labels("UBA226 sp9596u"), "UBA226 sp9596u")
})

test_that("unknowns filled by fill_GTDB_unknown_taxonomies() are rewritten", {
  expect_equal(enhance_taxonomy_labels("g__Natrinema..s__unknown"), "*Natrinema* sp. (unknown)")
  expect_equal(enhance_taxonomy_labels("f__Natrialbaceae..g__unknown"),
               "*Natrialbaceae* (unknown genus)")
  expect_equal(enhance_taxonomy_labels("o__Foo..f__unknown..g__unknown"),
               "*Foo* (unknown genus)")
  expect_equal(enhance_taxonomy_labels("g__UBA226..s__unknown"), "UBA226 sp. (unknown)")
})

test_that("MetaPhlAn-style names are handled", {
  expect_equal(enhance_taxonomy_labels("Dialister_invisus"), "*Dialister invisus*")
  expect_equal(enhance_taxonomy_labels("Actinomyces_sp_ICM47"), "*Actinomyces* sp. ICM47")
  expect_equal(enhance_taxonomy_labels("Rothia_SGB49305"), "*Rothia* SGB49305")
  expect_equal(enhance_taxonomy_labels("GGB3345_SGB4425"), "GGB3345 SGB4425")
  expect_equal(enhance_taxonomy_labels("TM7_phylum_sp_oral_taxon_352"),
               "TM7 phylum sp. oral taxon 352")
  expect_equal(enhance_taxonomy_labels("Lachnospiraceae_bacterium"), "*Lachnospiraceae* bacterium")
  expect_equal(enhance_taxonomy_labels("Candidatus_Woesearchaeota_archaeon_CG_4_10"),
               "*Candidatus Woesearchaeota* archaeon CG 4 10")
  expect_equal(enhance_taxonomy_labels("Burkholderiales_proteobacterium"),
               "*Burkholderiales* proteobacterium")
  expect_equal(enhance_taxonomy_labels("Wolbachia_endosymbiont"), "*Wolbachia* endosymbiont")
  expect_equal(enhance_taxonomy_labels(c("UNCLASSIFIED", "unknown", "t__SGB4933")),
               c("UNCLASSIFIED", "unknown", "SGB4933"))
})

test_that("bracket tags, NA, empty strings, factors and HTML markup", {
  expect_equal(enhance_taxonomy_labels("Agathobacter rectalis [SGB4933]"),
               "*Agathobacter rectalis* [SGB4933]")
  expect_equal(enhance_taxonomy_labels("Blautia_A faecis [SGB4810_group]"),
               "*Blautia* (A) *faecis* [SGB4810_group]")
  expect_equal(enhance_taxonomy_labels(c("Blautia faecis", NA, "")),
               c("*Blautia faecis*", NA, ""))
  expect_equal(enhance_taxonomy_labels(factor("Blautia faecis")), "*Blautia faecis*")
  expect_equal(enhance_taxonomy_labels("s__Blautia_A faecis", markup = "html"),
               "<i>Blautia</i> (A) <i>faecis</i>")
  expect_error(enhance_taxonomy_labels("Blautia faecis", markup = "latex"))
})

test_that("every name in the GTDB 202403 (= 202307) taxonomy is formatted consistently", {
  tax <- load_taxonomy("GTDB", chocophlan_datestamp = "202403")
  ranks <- setdiff(colnames(tax), "SGB")
  all_names <- unique(unlist(tax[ranks], use.names = FALSE))
  labels <- enhance_taxonomy_labels(all_names)

  # balanced italics everywhere
  n_stars <- lengths(regmatches(labels, gregexpr("*", labels, fixed = TRUE)))
  expect_true(all(n_stars %% 2 == 0))

  # without markup, a label is the input with prefixes dropped and suffixes in ()
  plain <- !grepl("..", all_names, fixed = TRUE)
  expected <- gsub("_([A-Z]{1,2})(\\s|$)", " (\\1)\\2", sub("^[a-z]__", "", all_names[plain]))
  expect_equal(gsub("*", "", labels[plain], fixed = TRUE), expected)

  # species: Latin genus -> label starts italic; placeholder genus -> roman
  species <- sub("^s__", "", unique(tax$Species[!grepl("..", tax$Species, fixed = TRUE)]))
  species <- species[nzchar(species) & species != "UNCLASSIFIED"]
  genus <- sub("_[A-Z]{1,2}$", "", sub(" .*", "", species))
  species_labels <- enhance_taxonomy_labels(species)
  latin_genus <- grepl("^[A-Z][a-z]+$", genus)
  expect_true(all(startsWith(species_labels[latin_genus], "*")))
  expect_false(any(grepl("*", species_labels[!latin_genus], fixed = TRUE)))

  # sp<accession> epithets are never italicised
  expect_false(any(grepl("\\*sp[0-9G]|sp[0-9G][^ ]*\\*", species_labels)))
  expect_true(any(grepl("^\\*Blautia\\* sp[0-9]", species_labels)))
})
