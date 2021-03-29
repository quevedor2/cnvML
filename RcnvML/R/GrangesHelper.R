#' Calculate cumulative genomic positions
#'
#' @param gr0 GenomicRanges object to get the cumulative positions for
#' @param seqstyle seq style from GenomeInfoDb (Default='ENSEMBL')
#'
#' @return
#' GRanges object with the cumulative postiions added to it
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomeInfoDb seqnames
#' @importFrom methods as
#' @export
calcCumulativeGPos <- function(gr0, seqstyle='ENSEMBL'){
  chrs <- getChrLength()
  seqlevelsStyle(chrs) <- seqstyle
  
  gr0 <- unlist(as(lapply(split(gr0, seqnames(gr0)), function(chr_gr){
    # Matches chr position
    chr_id <- unique(as.character(seqnames(chr_gr)))
    ref_idx <- match(chr_id, as.character(seqnames(chrs)))
    # Calculates cumulative position
    chr_gr$cum.start <- start(chr_gr) + chrs[ref_idx,]$cum.start - 1
    chr_gr$cum.end<- end(chr_gr) + chrs[ref_idx,]$cum.start - 1
    return(chr_gr)
  }), "GRangesList"))
  return(gr0)
}

#' Intersect and annotate sample with reference GR
#' @description Intersects the sample GR (s_gr) with the 
#' reference GR (ref_gr) and annotates all of the overlaps
#' in the sample GR with the reference
#' 
#' @param s_gr Sample GenomicRanges object
#' @param ref_gr Reference GenomicRanges object
#' @param cpos cpos
#'
#' @return
#' A GRanges intersected between s_gr and ref_gr with annotations
#' from ref_gr
#' @importFrom assertthat assert_that
#' @importFrom S4Vectors subjectHits
#' @import GenomicRanges
#' @export
intersectAndAnnotate <- function(s_gr, ref_gr, cpos=TRUE){
  ## Intersects s_gr with ref_gr, then annotates all of s_gr component
  ## on the new intersect
  print(unique(s_gr$ID))
  end(s_gr) <- end(s_gr)-1 # Prevents merging of adjacent segments
  int_gr <- intersect(ref_gr, s_gr, ignore.strand=FALSE)
  ov_idx <- findOverlaps(int_gr, s_gr)
  mcols(int_gr) <- mcols(s_gr[subjectHits(ov_idx),])
  if(cpos) int_gr <- calcCumulativeGPos(int_gr)
  
  return(int_gr)
}