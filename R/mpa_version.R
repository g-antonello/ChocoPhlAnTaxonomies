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

# validate mpa version
# This function is supposed to validate in-between cases where new versions of 
# mpa/ChocoPhlAndo not correspond to a different taxonomy table. 

validate_mpa_version <- function(mpa_version){
  
  # force the version variable to character
  mpa_version <- as.character(mpa_version)
  
  if(mpa_version == "202403"){
    # this mpa version exists but its taxonomy has not changed since 
    # the previous release
    mpa_version <- "202307"
  }
  
  # future exceptions can be written here
  
  return(mpa_version)
}
