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


#' Add an apha layer to the matrix
#'
#' @param matX A matrix of RGB colours ranging in value from 0 to 255
#' @param col Colour to add to the matrix
#' @param alpha Alpha value of the coloure
#' @importFrom abind abind
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices col2rgb
#' @export
alphaIt <- function(matX, col, alpha=0.5){
  #x <- abs(1-(matX/255))
  x <- (matX/255)
  
  col_ramp <- colorRampPalette(colors = c("white", col))(255)
  feat_cols <- col_ramp[matX+1]
  feat_rgb <- col2rgb(feat_cols)
  
  x2 <- abind(matrix(feat_rgb[1,]/255, ncol=ncol(matX)), 
              matrix(feat_rgb[2,]/255, ncol=ncol(matX)), 
              matrix(feat_rgb[3,]/255, ncol=ncol(matX)), 
              x/(1/alpha), along=3 )
  return(x2)
}

#' Add cumulative genomic position
#'
#' @param dat Input data, either a seg dataframe or data
#' @param ref GRanges object from getChrLength
#' @param dat.type either 'data' or 'seg'
#' @importFrom BiocGenerics match
#' @export
addCumPos <- function(dat, ref, dat.type){
  m.row.idx <- match(as.character(dat$chrom), as.character(seqnames(ref)))
  if(dat.type=='data'){
    dat$cpos <- ref[m.row.idx,]$cum.start +  dat$pos - 1
    dat$chr.stat
  } else if(dat.type == 'seg'){
    dat$cloc.start <- ref[m.row.idx,]$cum.start +  dat$loc.start - 1
    dat$cloc.end <- ref[m.row.idx,]$cum.start +  dat$loc.end - 1
  }
  dat$chr.stat <- (m.row.idx %% 2) + 1
  return(dat)
}
