#' Wallen ZD 2022 MetaPhlAn 3.0 relative abundance example data
#'
#' Relative abundance, taxonomy and sample metadata (colData) in 
#' `TreeSummarizedExperiment` format. The data were processed with MetaPhlAn
#' 3.0, using the 201901 database. rows are taxa, columns are samples. 
#' Paper: https://doi.org/10.1038/s41467-022-34667-x
#'
#' @format ##
#'    class: TreeSummarizedExperiment 
#'    dim: 719 724 
#'    metadata(0):
#'      assays(1): metaphlan
#'    rownames(719): s__Methanobacterium_sp_MB1 s__Methanobrevibacter_arboriphilus ...
#'    s__Saccharomyces_cerevisiae s__Blastocystis_sp_subtype_1
#'    rowData names(8): kingdom phylum ... species clade_name
#'    colnames(724): DP310 DP309 ... DC169 DC168
#'    colData names(57): sample_name Case_status ... collection_method total_sequences
#'    reducedDimNames(0):
#'      mainExpName: NULL
#'    altExpNames(6): kingdom phylum ... family genus
#'    rowLinks: NULL
#'    rowTree: NULL
#'    colLinks: NULL
#'    colTree: NULL
#'   
#' @source <https://zenodo.org/records/7246185>
"WallenZD_2022.tse"


#' JoS 2022 MetaPhlAn 4.1 relative abundance example data
#'
#' Relative abundance, taxonomy and sample metadata (colData) in 
#' `TreeSummarizedExperiment` format. The data were processed with MetaPhlAn
#' 4.1, using the 202403 database. rows are taxa, columns are samples. 
#' Paper: https://doi.org/10.1038/s41531-022-00351-6
#'
#' @format ##
#'    class: TreeSummarizedExperiment 
#'    dim: 1785 156 
#'    metadata(1): rowData_extraColumns
#'    assays(3): relative_abundance percent counts
#'    rownames(1785): t__SGB15291 t__SGB2328 ... t__SGB53574 t__SGB9634
#'    rowData names(8): Kingdom Phylum ... Species SGB
#'    colnames(156): 367613f5-f3a3-4448-99da-4d204b93b422 36cd8d84-758f-4c16-b0e5-40aa5446c347 ...
#'    ff1bf78d-3711-43b9-98ae-21337d309c9b ffb05e68-a5df-4cb9-9895-0053a33dea3c
#'    colData names(60): uuid db_version ... JoS_2022_uncurated_SRA.Study number_reads
#'    reducedDimNames(0):
#'      mainExpName: NULL
#'    altExpNames(0):
#'      rowLinks: NULL
#'    rowTree: NULL
#'    colLinks: NULL
#'    colTree: NULL
#'   
#' @source <https://github.com/ASAP-MAC/parkinsonsMetagenomicData>
"JoS_2022.tse"


