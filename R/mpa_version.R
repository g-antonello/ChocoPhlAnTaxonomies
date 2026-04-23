#' Get latest ChocoPhlAn release version
#'
#' Connects to the biobakery repository to fetch the latest available
#' MetaPhlAn/ChocoPhlAn database version name and extracts its datestamp.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{name}: The full release version name (character).
#'   \item \code{datestamp}: The extracted date string (character).
#' }
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # This requires an internet connection
#' mpa_latest <- get_mpa_latest()
#' print(mpa_latest$datestamp)
#' }
get_mpa_latest <- function(){
  cmprod_url <- "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases"
  mpa_latest <- readLines(sprintf("%s/mpa_latest", cmprod_url))
  return(list(
    name = mpa_latest, 
    datestamp = sapply(strsplit(mpa_latest, "_"), function(x) x[length(x)])))
}

#' Validate and map MetaPhlAn database versions
#'
#' Some MetaPhlAn/ChocoPhlAn releases do not correspond to a new taxonomy table.
#' This function validates the requested version against known exceptions and 
#' maps newer versions to their stable, existing counterparts to ensure 
#' compatibility with the available taxonomy data.
#'
#' @param chocophlan_version A string or numeric representing the MetaPhlAn database 
#'   version (e.g., "202403").
#'
#' @return A character string representing the validated taxonomy version 
#'   datestamp to be used.
#' @export
#'
#' @examples
#' # Maps the 202403 release to the stable 202307 taxonomy
#' validate_chocophlan_version("202403") 
#' 
#' # Returns the original version if no mapping exists
#' validate_chocophlan_version("202307")
#'
validate_chocophlan_version <- function(chocophlan_version){
  
  # in case of 'latest', fetch the latest version and make it into a timestamp
  if (chocophlan_version == "latest") {
    chocophlan_version <- get_mpa_latest()$datestamp
  } 
  
  # in case of 202403, the taxonomies are different in some alternate names, but 
  # the core remains the same. the arbitrary choice for now is to convert 202403
  # into 202307
  if(chocophlan_version == "202403"){
    chocophlan_version <- "202307"
  }
  
  return(chocophlan_version)
}