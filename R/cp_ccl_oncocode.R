#### Purpose ####
##################
## Uses the onco_meta_df.rds file made from map_ccl_oncocode.R to move the 
# hilbert images of cell line data into their own respective directories.
# First split by dataset (i.e. gCSI, CCLE, GDSC), then by cancer
# type (i.e. BRCA, CHOL, LUAD, ...)

#### Main ####
##############
library(filesstrings)
library(RcnvML)


PDIR <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/ccl_aggregate'
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



