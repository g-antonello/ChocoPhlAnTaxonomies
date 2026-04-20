#' Get latest ChocoPhlAn release version
#'
#' @returns list of full version name and datestamp
#' @export
#'
#' @examples
#' mpa_latest <- get_mpa_latest()
#' 
get_mpa_latest <- function(){
  cmprod_url <- "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases"
  mpa_latest <- readLines(sprintf("%s/mpa_latest", cmprod_url))
  return(list(
    name = mpa_latest, 
    datestamp = sapply(strsplit(mpa_latest, "_"), function(x) x[length(x)])))
}
