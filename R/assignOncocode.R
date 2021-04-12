# Description: Uses the onco_meta_df.rds file made from map_ccl_oncocode.R to move 
## the hilbert images of cell line data into their own respective directories.
## First split by dataset (i.e. gCSI, CCLE, GDSC), then by cancer type (i.e. BRCA, 
## CHOL, LUAD, ...)

# Usage:
## Rscript assignOncocode.R --cntype 'ASCN' --sfc sweep --dataset ccl_aggregate
## Rscript assignOncocode.R --cntype 'ASCN' --sfc sweep --dataset TCGA
## Rscript assignOncocode.R --cntype 'ASCN' --sfc hilbert --dataset ccl_aggregate
## Rscript assignOncocode.R --cntype 'ASCN' --sfc hilbert --dataset TCGA
## Rscript assignOncocode.R --cntype 'TCN' --sfc hilbert --dataset ccl_aggregate
## Rscript assignOncocode.R --cntype 'TCN' --sfc hilbert --dataset TCGA

suppressPackageStartupMessages(library(optparse))
####################
#### Parameters ####
option_list <- list( 
  make_option(c("-c", "--cntype"), type="character", default='ASCN',
              help="Copy number to plot: TCN or ASCN [%default]"),
  make_option(c("-p", "--pdir"), type="character", 
              default='/mnt/work1/users/pughlab/projects/cancer_cell_lines',
              help="Path to seg files [dir]/[dataset]/input/[seg_file] [%default]"),
  make_option(c("-s", "--sfc"), type="character", default='hilbert',
              help="Space-filling curve to use, 'sweep' or 'hilbert' [%default]"),
  make_option(c("-d", "--dataset"), type="character", default='ccl_aggregate',
              help="Dataset to use, 'ccl_aggregate' or 'TCGA' [%default]"),
  make_option(c("-v", "--verbose"), action="store_false", default=FALSE,
              help="Print extra output [%default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

## Load libraries
libraries <- c('filesstrings', 'RcnvML', 'assertthat')
null <- lapply(libraries, function(i) suppressMessages(require(i, character.only = TRUE)))
# Set up parameters
PDIR <- opt$pdir          # Path to seg files [file.path(PDIR, analysis, input, seg_i)]
analysis <- opt$dataset   # dataset to analyze (TCGA or ccl_aggregate)
sfc <- opt$sfc            # whether to use a hilbert SFC or sweep
cntype <- opt$cntype      # Plot ASCN or TCN as colours

# Asserts
assert_that(length(cntype)==1, any(c('TCN', 'ASCN') %in% cntype), 
            msg='--cntype must be either "TCN" or"ASCN"')
assert_that(length(sfc)==1, any(c('sweep', 'hilbert') %in% sfc), 
            msg='--sfc must be either "sweep" or "hilbert"')
assert_that(length(analysis)==1, any(c('ccl_aggregate', 'TCGA') %in% analysis), 
            msg='--dataset must be either "ccl_aggregate" or "TCGA"')

#### Main ####
##############
setwd(file.path(PDIR, analysis))
data("onco_meta_df")

## List all hilbert files
fls <- list.files(file.path("output", sfc, cntype), pattern="png$")
assert_that(length(fls) > 1, 
            msg='Could not locate SFC png images outputed from plotSFC.R')

if(analysis == 'ccl_aggregate'){
  ## Assigns the CCL files to different oncocode directories
  datasets <- c('GDSC', 'CCLE', 'GNE')
  onco_meta_df$GDSC <- gsub(".cel$", "", onco_meta_df$GDSC, ignore.case = TRUE)
  
  copy_status <- sapply(fls, function(fl_i){
    # Map file to the onco_meta_df
    print(fl_i)
    i <- gsub(".png$", "", fl_i)
    idx <- which(i == onco_meta_df[,datasets], arr.ind=TRUE)
    
    # Identity dataset and oncocode
    ds <- datasets[idx[,2]]
    oncocode <- as.character(gsub(":.*", "", onco_meta_df[idx[,1],]$oncocode))
    
    # Copy file to new locations
    if(length(ds) > 0 & length(oncocode) > 0){
      target_dir = file.path("output", sfc, cntype, ds, oncocode)
    } else {
      target_dir = file.path("output", sfc, cntype, "unknown")
    }
    if(!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
    file.copy(file.path("output", sfc, cntype, fl_i),
              file.path(target_dir, fl_i), overwrite = TRUE)
  })
} else if(analysis=='TCGA'){
  
}
  