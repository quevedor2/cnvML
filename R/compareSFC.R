# Purpose: cnvML-Supplementary Figure 1

# Description:  Compare Space-Filling Curves (SFC)
## This script is meant to take the cytoband information
## from the hg19 genome and encode the one-dimensional genome
## into a 2D representation using the Hilbert space-filling
## curve and a Sweep space-filling curve.  Additionally, this 
## script will generate the euclidean distance between regions 
## in the 2D representation.

library(RColorBrewer)
library(HilbertCurve)
library(HilbertVis)
library(GenomicRanges)
library(RcnvML)   # devtools::install_github("quevedor2/cnvML/RcnvML")

## 
# This directory just guides the output, there are no inputs
# for this script
PDIR <-'/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnvML/hilbert_modeling'
setwd(PDIR)

################
####  Main  ####
for(order in c(4,6,8)){
  #############################
  #### Set up HilbertCurve ####
  gbin_pos <- setupRefHcMatrix(order=order)
  gbin_pos$ord$chr <- paste0("chr", gsub(":.*", "", gbin_pos$ord$start))
  gbin_pos$ord$loc.start <- as.integer(gsub("^.*:", "", gbin_pos$ord$start))
  gbin_pos$ord$loc.end <- as.integer(gsub("^.*:", "", gbin_pos$ord$end))
  
  ###################################
  #### Set up Sweep Genomic Bins ####
  bins <- mapSFC(sfc='sweep', hc_ord=gbin_pos$ord)
  
  #########################################################
  #### Visualize difference using HilbertCurve package ####
  chr.size.dat <- getChrLength()
  seqlevelsStyle(chr.size.dat) <- 'NCBI'
  scale <- 1*10^4
  chr_ranges <- as.data.frame(mcols(chr.size.dat))[,c('cum.start', 'cum.end')]
  ir <- IRanges::IRanges(start = as.numeric(chr_ranges$cum.start)/scale,
                         end = as.numeric(chr_ranges$cum.end)/scale)
  fgfr3 <- IRanges::IRanges(start=(1795020+chr.size.dat['chr4',]$cum.start)/scale,
                            end=(1810594+chr.size.dat['chr4',]$cum.start)/scale)
  
  # Chrom Colors
  n <- 24
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  # Genomic Tiling
  pdf(paste0("chrom-tile_ord", order, ".pdf"))
  hc <- genHC(max(chr.size.dat$cum.end)/scale, order,
              sfc_pos=bins[-nrow(bins),c('x1', 'y1', 'x2', 'y2')])
  # hc@POS <- bins[-nrow(bins),c('x1', 'y1', 'x2', 'y2')]
  hc_layer(hc, ir, mean_mode = "absolute", col = col_vector[1:24])
  hc_layer(hc, fgfr3, mean_mode = "absolute", col = 'black')
  dev.off()
  
  # Hilbert Curve
  pdf(paste0("chrom-hc_ord", order, ".pdf"))
  hc <- genHC(max(chr.size.dat$cum.end)/scale, order)
  hc_layer(hc, ir, mean_mode = "absolute", col = col_vector[1:24])
  hc_layer(hc, fgfr3, mean_mode = "absolute", col = 'black')
  dev.off()
  
  
  
  ########################################################
  #### Map distance between adjacent regions on image ####
  ## bins 2d distance
  gpos <- gbin_pos$ord
  gpos <- gpos[order(gpos$gord),]
  pos = HilbertVis::hilbertCurve(order)
  
  st_frac  <- 1-round(chr.size.dat$cum.start / max(chr.size.dat$cum.end),3)
  end_frac <- 1-round(chr.size.dat$cum.end / max(chr.size.dat$cum.end),3)
  frac_df <- data.frame("start"=st_frac, "end"=end_frac, 
                        "col"=as.character(col_vector[1:24]), stringsAsFactors = F)
  
  px <- 600
  xdist=as.matrix(dist(bins[,c('x1', 'y1')]))
  png(paste0("tile-dist_ord", order, ".png"), width = px, height = px, units = "px")
  plotEuclidDist(xdist, frac_df)
  dev.off()
  
  xdist=as.matrix(dist(gpos[,c('x1', 'y1')]))
  png(paste0("hc-dist_ord", order, ".png"), width = px, height = px, units = "px")
  plotEuclidDist(xdist, frac_df)
  dev.off()
}
