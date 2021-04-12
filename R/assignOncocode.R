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

# Set up parameters
PDIR <- opt$pdir                   # Path to seg files [file.path(PDIR, analysis, input, seg_i)]
analysis <- opt$dataset           # dataset to analyze (TCGA or ccl_aggregate)
sfc <- opt$sfc                    # whether to use a hilbert SFC or sweep
cntype <- opt$cntype              # Plot ASCN or TCN as colours

# Asserts
assert_that(length(cntype)==1, any(c('TCN', 'ASCN') %in% cntype), 
            msg='--cntype must be either "TCN" or"ASCN"')
assert_that(length(sfc)==1, any(c('sweep', 'hilbert') %in% sfc), 
            msg='--sfc must be either "sweep" or "hilbert"')
assert_that(length(analysis)==1, any(c('ccl_aggregate', 'TCGA') %in% analysis), 
            msg='--dataset must be either "ccl_aggregate" or "TCGA"')

#### Main ####
##############
library(filesstrings)
library(RcnvML)

setwd(file.path(pdir, analysis))
data("onco_meta_df")
# meta_df <- readRDS(file.path(PDIR, "ref", "onco_meta_df.rds"))

## List all hilbert files
fls <- list.files(file.path(PDIR, "output", "hilbert"), pattern="png$")
datasets <- c('GDSC', 'CCLE', 'GNE')
meta_df$GDSC <- gsub(".cel$", "", meta_df$GDSC, ignore.case = TRUE)

copy_status <- sapply(fls, function(fl_i){
  # Map file to the meta_df
  print(fl_i)
  i <- gsub(".png$", "", fl_i)
  idx <- which(i == meta_df[,datasets], arr.ind=TRUE)
  
  # Identity dataset and oncocode
  ds <- datasets[idx[,2]]
  oncocode <- as.character(gsub(":.*", "", meta_df[idx[,1],]$oncocode))
  
  # Copy file to new locations
  if(length(ds) > 0 & length(oncocode) > 0){
    target_dir = file.path(PDIR, "output", "hilbert", ds, oncocode)
  } else {
    target_dir = file.path(PDIR, "output", "hilbert", "unknown")
  }
  if(!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
  file.copy(file.path(PDIR, "output", "hilbert", fl_i),
            file.path(target_dir, fl_i), overwrite = TRUE)
})



