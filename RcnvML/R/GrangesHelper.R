#' Calculate cumulative genomic positions
#'
#' @param gr0 GenomicRanges object to get the cumulative positions for
#' @param seqstyle seq style from GenomeInfoDb (Default='ENSEMBL')
#'
#' @return
#' GRanges object with the cumulative postiions added to it
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomeInfoDb seqnames
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
