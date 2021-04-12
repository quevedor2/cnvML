library(GenomicRanges)
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
  bins <- mapSFC(sfc='sweep', hc_ord=gbin_pos$ord)
  
  pdf(file.path("arrows", paste0("sweep_", ord, ".pdf")))
  hc <- genHC(max(chr.size.dat$cum.end)/scale, ord, mode='normal',
              sfc_pos=bins[-nrow(bins),c('x1', 'y1', 'x2', 'y2')])
  dev.off()
  
  pdf(file.path("arrows", paste0("hilbert_", ord, ".pdf")))
  hc <- genHC(max(chr.size.dat$cum.end)/scale, ord, mode='normal')
  dev.off()
})

