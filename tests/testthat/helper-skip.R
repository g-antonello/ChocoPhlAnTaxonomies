# Tests that need the biobakery server (validate_chocophlan_version("latest"))
skip_if_no_biobakery <- function() {
  skip_on_cran()
  skip_if_offline("cmprod1.cibio.unitn.it")
}
