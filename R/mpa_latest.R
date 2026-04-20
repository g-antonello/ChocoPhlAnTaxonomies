get_mpa_latest <- function(){
  cmprod_url <- "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases"
  mpa_latest <- readLines(sprintf("%s/mpa_latest", cmprod_url))
  return(mpa_latest)
}

mpa_latest <- get_mpa_latest()