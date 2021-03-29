#' hg19 chromosome sizes
#' @description Returns a GRanges object of the hg19 chromosome sizes
#' as well as their cumulative starts, ends, and midpoints
#' @return 
#' GenomicRanges object for the hg19 chromosomes and cumulative start/end points
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 Hsapiens
#' @export
getChrLength <- function(){
  #require(BSgenome.Hsapiens.UCSC.hg19)
  chr.lengths = seqlengths(Hsapiens)[1:24]
  chr.len.gr <- makeGRangesFromDataFrame(data.frame("chrom"=names(chr.lengths),
                                                    "loc.start"=rep(1, length(chr.lengths)),
                                                    "loc.end"=chr.lengths))
  chr.len.gr$cum.end <- cumsum(as.numeric(end(chr.len.gr)))
  chr.len.gr$cum.start <- chr.len.gr$cum.end - (end(chr.len.gr) -1)
  chr.len.gr$cum.mid <- (chr.len.gr$cum.start + ((chr.len.gr$cum.end - chr.len.gr$cum.start)/2))
  return(chr.len.gr)
}
