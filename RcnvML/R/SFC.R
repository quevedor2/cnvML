#' Initiailize HilbertCurve object
#'
#' @description Wrapper function to run the HilbertCurve() function
#' using only a set end point and order
#'
#' @param end Typically, the size of the genome [integer]
#' @param mode 'pixel' or 'normal', default is 'pixel'
#' @param sfc_pos POS values for SFC mapping [dataframe]
#' @param order Order to generate hilbert curve for [integer]
#'
#' @return
#' an object from HilbertCurve::HilbertCurve()
#' @export
#' @importFrom HilbertCurve HilbertCurve
#' @importFrom assertthat assert_that
#' @importFrom intervals Intervals
#' @importFrom utils tail
#' @importFrom grid unit
#' 
genHC <- function(end, order, mode='pixel', sfc_pos=NULL){
  assert_that(is.numeric(end), end > 1, length(end)==1,
              msg="'end' must be a single integer greater than 1")
  assert_that(is.integer(order), end > 2, length(end)==1,
              msg="'order' must be a single integer greater than 2")
  
  # SfcCurve is a copy of HilbertCurve::HilbertCurve() function but allows
  # for an import of @POS values for a different SFC mapping
  hc = SfcCurve(1, end, level = order, mode = mode, 
                reference = TRUE, padding=unit(1, "mm"),
                newpage = FALSE, sfc_pos=sfc_pos)
  return(hc)
}

#' Parse the HilbertCurve object
#'
#' @description Wrapper function to run the HilbertCurve() function
#' but also to parse the returning object into a dataframe that 
#' maps the coordinates to genomic positions
#' @param order Order to generate the HilbertCurve object [integer]
#' @param scale divide the end by a given scale to adjust [integer]
#' @return
#'  A list object with 3 elements:
#'    - 'hc' = HilbertCurve object
#'    - 'mat' = Square matrix with unique X_Y coordinates (e.g. 1_2)
#'    - 'ord' = A data frame that maps to the 'mat' matrix. It dictates
#'    the genomic position for each unique X_Y coordinate
#' @export
#' @importFrom assertthat assert_that
#' @importFrom GenomeInfoDb seqlevelsStyle<-
#' @importFrom intervals Intervals
#' @importFrom intervals interval_overlap
#' @import GenomicRanges
setupRefHcMatrix <- function(order=8, scale=1){
  chr.size.dat <- getChrLength()
  seqlevelsStyle(chr.size.dat) <- 'NCBI'
  
  hc <- genHC(end = as.numeric(max(chr.size.dat$cum.end)/scale), 
              order = as.integer(order))
  
  ## Identify the intervals for each HC bin by dividing by the zoom factor
  bin_df <- data.frame(start = round(start(hc@BINS) / hc@ZOOM),
                       end = round(end(hc@BINS) / hc@ZOOM))
  
  ## For each interval, identify the corresponding genomic position by mapping cumulative pos
  ### Create interval objects of bins and chromosomes
  bin_df[,1] <- bin_df[,1]+1
  bin_ir <- Intervals(bin_df)            # bins intervals
  bin_chrsize <- as.data.frame(mcols(chr.size.dat)[,c('cum.start', 'cum.end')])/scale
  chrsize_ir <- Intervals(bin_chrsize)   # chrs intervals
  
  ### overlap the bin and chrs intervals and flag intervals that span two chrs
  ov_idx <- interval_overlap(bin_ir, chrsize_ir)
  spans_chrs <- which(sapply(ov_idx, length) > 1)
  
  ### Find the cumulative start/end pos for each bin according to its chr
  cumstart_pos <- sapply(ov_idx, function(i) i[1])
  adjstart_pos <- as.matrix(chrsize_ir[cumstart_pos])[,1] - 1
  adjend_pos <- adjstart_pos
  adjend_pos[spans_chrs] <- adjend_pos[spans_chrs+1]
  
  ### assemble genomic loci by substracting cumulative pos from chr position
  gbin_df <- bin_df
  gbin_df$loc.start <-  bin_df$start - adjstart_pos
  gbin_df$loc.end <-  bin_df$end - adjend_pos
  names(adjstart_pos) <- gsub("^chr", "", names(adjstart_pos), ignore.case = T)
  gbin_df$end.chr <-  gbin_df$start.chr <- names(adjstart_pos)
  gbin_df$end.chr[spans_chrs] <- names(adjstart_pos)[spans_chrs+1]
  gbin_df$start <- paste0(gbin_df$start.chr, ":", gbin_df$loc.start)
  gbin_df$end <- paste0(gbin_df$end.chr, ":", gbin_df$loc.end)
  
  ## Combine genomic position with HC matrix pos
  gbin_pos <- cbind(gbin_df, hc@POS)
  gbin_pos$gord <- c(1:nrow(gbin_pos))
  gbin_pos_ord <- gbin_pos[order(gbin_pos$y1, decreasing = TRUE),]
  gbin_pos_ord <- gbin_pos_ord[order(gbin_pos_ord$x1, decreasing=FALSE),]
  ## Associate position in matrix (UID) with a genomic position (start, end)
  gbin_pos_ord$uid <- paste(gbin_pos_ord$x1, gbin_pos_ord$y1, sep="_")
  
  ## Add a filler for the last interval missing from the order dataframe
  last_start <- tail(gbin_pos$loc.end,1)
  last_chr <- tail(gbin_pos$end.chr,1)
  chr_idx <- which(as.character(seqnames(chr.size.dat)) == last_chr)
  last_end <- end(chr.size.dat[chr_idx,])/scale
  gbin_pos_ord <- rbind(gbin_pos_ord, 
                        data.frame('start'=paste0(last_chr, ":", last_start+1),
                                   'end'=paste0(last_chr, ":", last_end),
                                   'loc.start'=last_start+1, 'loc.end'=last_end,
                                   'start.chr'=last_chr, 'end.chr'=last_chr, 
                                   'x1'=max(gbin_pos_ord$x1), 'y1'=0, 
                                   'x2'=max(gbin_pos_ord$x1)+1, 'y2'=0, 
                                   'gord'=nrow(gbin_pos_ord)+1, 
                                   'uid'=paste0(max(gbin_pos_ord$x1), '_0')))
  
  ## Reform the HC matrix using UIDs instead of mapping information
  gbin_pos_mat <- matrix(gbin_pos_ord$uid, nrow=max(gbin_pos_ord$y1)+1, 
                         ncol=max(gbin_pos_ord$x1)+1)
  return(list("hc"=hc,
              "mat"=gbin_pos_mat, 
              "ord"=gbin_pos_ord))
}

#' Map a space-filling curve using a HilbertCurve object
#'
#' @description Maps a space-filling curve to the HilbertCurve object 
#' @param sfc space filling curve [character]
#' @param hc_ord dataframe object from setupRefHcMatrix(order=ord)$ord
#' @param order Order to run setupRefHcMatrix if hc_ord not given [integer]
#' @param uids Order to fill in the space-filling curve to avoid regenerating
#' @return
#' a dataframe object of a reordered X_Y positions correspondign to a different
#' space-filling curve for the same square matrix
#' @export
#' @importFrom assertthat assert_that
#' 
mapSFC <- function(sfc='sweep', hc_ord=NULL, order=NULL, uids=NULL){
  assert_that(is.character(sfc), length(sfc)==1,
              msg="space-filling curve me be a single character vector")
  
  if(is.null(hc_ord)){
    ## Run HilbertCurve if no HC object is given
    assert_that(is.integer(order), msg="integer order must be set for hilbert curve")
    hc_ord <- setupRefHcMatrix(order=order)$ord
  }
  
  # assert proper input of hc_ord
  assert_that(is.data.frame(hc_ord), 
              all(c('gord', "x1", "x2", "y1", "y2") %in% colnames(hc_ord)),
              msg="'hc_ord' is malformed. You need to rerun setupRefHcMatrix()")
  ords <- 4^(1:20)
  maxn <- 4^which.min(abs(ords-max(hc_ord$x1))) -1
  
  # If a mapping is given, ensure its proper format
  if(!is.null(uids)) assert_that(is.character(uids), 
                                 length(uids) == (maxn+1)^2, 
                                 all(grepl("[0-9]*_[0-9]*", uids)),
                                 msg="'uids' is malformed")
  
  if(grepl('^sweep$', sfc, ignore.case = T)){
    bins <- hc_ord
    bins <- bins[order(bins$gord),]
    
    # Create the Sweep order
    if(is.null(uids)){
      uids <- apply(expand.grid(c(0:maxn), 
                                c(0:maxn)), 
                    1, paste, collapse="_")
    }
    bins$uid <- uids 
    bins$x1 <- as.integer(gsub("_.*", "", uids))
    bins$y1 <- as.integer(gsub("^.*_", "", uids))
    bins$x2 <- c(bins$x1[-1], -1)
    bins$y2 <- c(bins$y1[-1], -1)
  }
  
  return('hc_ord'=bins)
}