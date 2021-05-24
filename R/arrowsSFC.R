# Purpose: cnvML-Supplementary Figure 1

# Description:  Map Space-Filling Curves (SFC)
## This script is meant to plot the mapping direction of a
## Hilbert space filling curve at 4th and 6th order, as well
## as a Sweep space filling curve at the same orders. It is
## designed to highlight how the 1D->2D mapping differs between
## the two approaches and not address the issues of locality
## or recursiveness

library(GenomicRanges)
library(RcnvML)
PDIR <-'/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnvML/hilbert_modeling'
setwd(PDIR)

###################################
#### Set up Sweep Genomic Bins ####
scale <- 1000
chr.size.dat <- getChrLength()
seqlevelsStyle(chr.size.dat) <- 'NCBI'

dir.create("arrows", recursive = F, showWarnings = F)
sapply(c(4,6,8), function(ord){
  ord <- as.integer(ord)
  gbin_pos <- setupRefHcMatrix(order=ord)
  na <- lapply(c('sweep', 'scan', 'diagonal', 'morton', 'random'), function(sfc){
    bins <- mapSFC(sfc=sfc, hc_ord=gbin_pos$ord)
    
    maxn <- sqrt(nrow(bins))
    
    
    
    pdf(file.path("arrows", paste0(sfc, "_", ord, ".pdf")))
    hc <- genHC(max(chr.size.dat$cum.end)/scale, ord, mode='normal',
                sfc_pos=bins[-nrow(bins),c('x1', 'y1', 'x2', 'y2')])
    dev.off()
  })
  
  seqinfo <- tryCatch({
    seqlevelsStyle(seqinfo) <- genomeStyle
    seqinfo
  }, error=function(e){
    if(genomeStyle =='NCBI'){
      seqnames(seqinfo) <- gsub("^chr", "", seqnames(seqinfo), ignore.case=TRUE)
      seqnames(seqinfo) <- gsub("^M$", "chrM", seqnames(seqinfo), ignore.case=TRUE)
      genome(seqinfo) <- 'GRCh37.p13'
    }
    seqinfo
  })
  
  pdf(file.path("arrows", paste0("hilbert_", ord, ".pdf")))
  hc <- genHC(max(chr.size.dat$cum.end)/scale, ord, mode='normal')
  dev.off()
})

