#' Mapping between cell line IDs, filenames and cancer type
#'
#' A data frame containing the mappings between cell line IDs
#' and their filenames in GDSC, CCLE, and gCSI (GNE). Additionally
#' the mapping of these cell lines to their cancer types as predicted
#' by either metadata annotation, or through the oncotree api
#' 
#' oncocodes genreated from: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
#' Cancer_passport IDs and TCGA codes were obtained from: https://www.cancerrxgene.org/celllines
#' Oncotree and CVCL ids from PharmacoGx mapping came from: https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/cello_to_onco_mapping_final_with_manual_review.csv
#'
#' @format A data frame with 1621 rows and 14 variables:
#' \describe{
#'   \item{simpleid}{cancer_type: TCGA Cancer study code}
#'   \item{ID}{patients: IDs of patients}
#'   \item{CVCL}{Cellosaurus IDs}
#'   \item{GDSC}{Filenames for GDSC dataset}
#'   \item{CCLE}{Filenames for CCLE dataset}
#'   \item{GNE}{Filenames for gCSI/GNE dataset}
#'   \item{PharmacoGX_ID}{Identifier for PharmacoGx}
#'   \item{SRR}{SRR ID}
#'   \item{EGAF}{EGAF file ID}
#'   \item{EGAN}{EGA sample ID}
#'   \item{gCSI_RNA}{gCSI filename RNA ID}
#'   \item{Var.12}{Unknown}
#'   \item{passport}{Cancer passport IDs from cellosaurus or cancerrxgene}
#'   \item{oncocode}{Oncode from oncotree or cancerrxgene census}
#' }
#' @source \url{https://www.cancerrxgene.org/celllines}
"onco_meta_df"

#' Mapping between TCGA IDs, cell line IDs, filenames and cancer type
#'
#' A list containing the onco_meta_df dataframe, or a TCGA metadata
#' dataframe containing the mapping between sample IDs and their
#' oncotree code.
#'
#' @format A list with 2 elements containing dataframes:
#' \describe{
#'   \item{ccl}{data frame of onco_meta_df}
#'   \item{tcga}{2 column dataframe of sampleID to oncocode}
#' }
"onco_meta"