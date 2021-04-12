#' Customized code to move files to their oncocode directory
#'
#' @param x sample id
#' @param meta_df metadata data frame
#' @param datasets vector of datasets to use
#' @param pattern Pattern to gsub file IDs
#' @param sfc SFC type (sweep or hilbert, path related)
#' @param cntype CN type (TCN or ASCN, path related)
#' @param rmfile remove the file after copying or not
#'
#' @export
mvFileToOncocode <- function(x, meta_df, datasets, pattern, 
                             sfc, cntype, rmfile=TRUE){
  # Map file to the onco_meta_df
  print(x)
  i <- gsub(pattern, "", x)
  idx <- which(i == meta_df[,datasets, drop=FALSE], arr.ind=TRUE)
  
  # Identity dataset and oncocode
  ds <- datasets[idx[,2]]
  oncocode <- as.character(gsub(":.*", "", meta_df[idx[,1],]$oncocode))
  
  # Copy file to new locations
  if(length(ds) > 0 & length(oncocode) > 0){
    target_dir = file.path("output", sfc, cntype, ds, oncocode)
  } else {
    target_dir = file.path("output", sfc, cntype, "unknown")
  }
  png_dir <- file.path("output", sfc, cntype, 'png')
  if(!dir.exists(png_dir)) dir.create(png_dir, recursive = TRUE)
  if(!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
  targstat  <- file.copy(file.path("output", sfc, cntype, x),
                         file.path(target_dir, x), overwrite = TRUE)
  pngstat <- file.copy(file.path("output", sfc, cntype, x),
                       png_dir, overwrite = TRUE)
  if(pngstat & targstat) file.remove(file.path("output", sfc, cntype, x))
}