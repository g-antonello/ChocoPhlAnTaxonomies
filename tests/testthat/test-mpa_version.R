test_that("validate_chocophlan_version maps 202403 to 202307 and passes others through", {
  expect_equal(validate_chocophlan_version("202403"), "202307")
  expect_equal(validate_chocophlan_version("202307"), "202307")
  expect_equal(validate_chocophlan_version("202503"), "202503")
  expect_equal(validate_chocophlan_version(202403), "202307")
})

test_that("'latest' is resolved through the biobakery server", {
  skip_if_no_biobakery()
  latest <- get_mpa_latest()
  expect_named(latest, c("name", "datestamp"))
  expect_match(latest$datestamp, "^[0-9]{6}$")
  expect_match(validate_chocophlan_version("latest"), "^[0-9]{6}$")
})
