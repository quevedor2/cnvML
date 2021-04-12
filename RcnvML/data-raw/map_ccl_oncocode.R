#### Purpose ####
##################
## Uses the cell line name and its mapped CVCL id (generated from the 
# CCLid paper) to map its corresponding cancer type and oncocode


# oncocodes genreated from: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
# Cancer_passport IDs and TCGA codes were obtained from: https://www.cancerrxgene.org/celllines
# 
# cellosaurus raw is documented in its respective directory
# onco_meta_df.rds generated from scripts/map_ccl_oncocode.R
# 
# Oncotree and CVCL ids from PharmacoGx mapping came from: https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/cello_to_onco_mapping_final_with_manual_review.csv



#### Function ####
##################
## oncocode_api: Takes in the name of a cancer type (parsed from cellosaurus
# dataset), and tries to use the oncotree API to map the 2nd level 
# oncocode (e.g., BRCA).  If that fails, it will return the 1st level 
# oncocode (e.g. BREAST). If no entry can be found, it will return
# an "Other:" tagged name (e.g. Other:Gingival squamous cell carcinoma)
oncocode_api <- function(name){
  #name="mantle cell lymphoma"
  base <- "http://oncotree.mskcc.org/api/tumorTypes/search/name/"
  name_api <- gsub(" ", "\\%20", tolower(name))
  suffix <- "?exactMatch=false&levels=1%2C2%2C3%2C4%2C5"
  
  call1 <- paste0(base, name_api, suffix)
  get_oncocode <- GET(call1)
  if(get_oncocode$status_code == 200){
    if(length(content(get_oncocode)) >=2){
      content(get_oncocode)[[2]]$parent
    } else {
      content(get_oncocode)[[1]]$parent
    }
  } else {
    paste0("Other:", name)
  }
}

#### MAIN ####
##################
library(httr)
library(jsonlite)

pdir <- '~/git/cnvML/RcnvML/data-raw'
setwd(pdir)

## Read in the datasets
# parsed Cellosaurus xml data in json format
load("cellosaurus_raw.rda")  ## cl.dat  
# Mapping of cell line, filenames, and CVCL ids
load('meta.df.rda')  ## meta.df
# TCGA Study abbreviation oncocodes (https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations)
oncocodes <- read.table('oncocodes.txt', header=FALSE, sep="\t",
                        stringsAsFactors = FALSE, check.names = FALSE, 
                        col.names = c('Code', 'Description'))
# Sanger CancerRxGene oncocodes associated with cell lines
cl_codes <- read.csv("cancerrxgene_celllines.csv", header=TRUE, 
                     stringsAsFactors = FALSE, check.names = FALSE)


## Using the meta.df data structure, maps the cvcl to the cancer type
# oncocode using everything available.
ccl_oncocode <- lapply(meta.df$CVCL, function(cvcl, verbose=FALSE){
  if(verbose) print(cvcl)
  ## Search the Cellosaurus database for the CancerRxGene 
  # cell-model-passport associated with the given CVCL id
  passports <- sapply(cl.dat[[cvcl]]$'xref-list', function(xref){
    if(xref$'.attrs'['database'] == 'Cell_Model_Passport'){
      xref$'.attrs'['accession']
    }
  })
  
  passport_idx <- sapply(passports, is.null)
  if(all(passport_idx)){
    ## CVCL ID is not found in CancerRxGene database
    passport <- NA
    ctype <- NULL
  } else {
    ## Use the CancerRxGene database to map oncocodes
    passport <- passports[[which(!passport_idx)[1]]]
    ctype <- cl_codes[match(passport, cl_codes$cell_model_passports),]$TCGA_Classification
    if(is.na(ctype)){
      ctype <- NULL
    } else if(ctype == 'UNCLASSIFIED') {
      ctype <- NULL
    }
    
  }
  
  ## Use the MSKCC oncotree api to get the parent oncocode for cell line
  if(is.null(ctype)) {
    name <- cl.dat[[cvcl]]$'disease-list'$'cv-term'$text
    if(is.null(name)){
      ctype <- paste0("Other:", 'UNCLASSIFIED')
    } else {
      ctype <- oncocode_api(name)
      Sys.sleep(1)
    }
  }
  
  
  return(data.frame("CVCL"=cvcl,
                    "passport"=passport,
                    "oncocode"=ctype))
}, verbose=TRUE)
ccl_oncocode <- as.data.frame(do.call(rbind, ccl_oncocode))




## Append and save data structure
meta.df <- cbind(meta.df, ccl_oncocode[,-1])
onco_meta_df <- meta.df
usethis::use_data(onco_meta_df, overwrite = T)
#saveRDS(meta.df, file=file.path(pdir, "onco_meta_df.rds"))