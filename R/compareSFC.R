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
library(ComplexHeatmap)     # Visualize distance between bins
library(intervals)
library(GenomicRanges)
library(AneuploidyScore)    # devtools::install_github("quevedor2/aneuploidy_score")

## 
# This directory just guides the output, there are no inputs
# for this script
PDIR <-'/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnvML/hilbert_modeling'
setwd(PDIR)

## Load cytoband/chr-arm data
data("ucsc.hg19.cytoband")  # hg19 cytoband from AneuploidyScore


cytoarm <- cytobandToArm(ucsc.hg19.cytoband)
cytoarm <- do.call(rbind, lapply(cytoarm, function(i){
  i[order(factor(rownames(i), c('p', 'cen', 'q'))),]
}))


###################
#### Functions ####
getChrLength <- function(){
  require(BSgenome.Hsapiens.UCSC.hg19)
  chr.lengths = seqlengths(Hsapiens)[1:24]
  chr.len.gr <- makeGRangesFromDataFrame(data.frame("chrom"=names(chr.lengths),
                                                    "loc.start"=rep(1, length(chr.lengths)),
                                                    "loc.end"=chr.lengths))
  chr.len.gr$cum.end <- cumsum(as.numeric(end(chr.len.gr)))
  chr.len.gr$cum.start <- chr.len.gr$cum.end - (end(chr.len.gr) -1)
  chr.len.gr$cum.mid <- (chr.len.gr$cum.start + ((chr.len.gr$cum.end - chr.len.gr$cum.start)/2))
  return(chr.len.gr)
}

addCumPos <- function(dat, ref, dat.type){
  if(class(seg) == 'GRanges'){
    m.row.idx <- match(as.character(seqnames(dat)), as.character(seqnames(ref)))
    
    dat$cloc.start <- ref[m.row.idx,]$cum.start + start(dat) - 1
    dat$cloc.end <- ref[m.row.idx,]$cum.start + end(dat) - 1
  } else {
    m.row.idx <- match(as.character(dat$chrom), as.character(seqnames(ref)))
    if(dat.type=='data'){
      dat$cpos <- ref[m.row.idx,]$cum.start +  dat$pos - 1
      dat$chr.stat
    } else if(dat.type == 'seg'){
      dat$cloc.start <- ref[m.row.idx,]$cum.start +  dat$loc.start - 1
      dat$cloc.end <- ref[m.row.idx,]$cum.start +  dat$loc.end - 1
    }
  }
  dat$chr.stat <- (m.row.idx %% 2) + 1
  return(dat)
}

genHC <- function(end, order){
  hc = HilbertCurve(1, end, level = order, mode = "pixel", 
                    reference = TRUE, padding=unit(1, "mm"),
                    newpage = FALSE)
  return(hc)
}

setupRefHcMatrix <- function(order=8){
  chr.size.dat <- getChrLength()
  seqlevelsStyle(chr.size.dat) <- 'NCBI'
  
  hc <- genHC(max(chr.size.dat$cum.end), order)

  ## Identify the intervals for each HC bin by dividing by the zoom factor
  bin_df <- data.frame(start = round(start(hc@BINS) / hc@ZOOM),
                       end = round(end(hc@BINS) / hc@ZOOM))
  
  ## For each interval, identify the corresponding genomic position by mapping cumulative pos
  ### Create interval objects of bins and chromosomes
  bin_df[,1] <- bin_df[,1]+1
  bin_ir <- intervals::Intervals(bin_df)            # bins intervals
  bin_chrsize <- as.data.frame(mcols(chr.size.dat)[,c('cum.start', 'cum.end')])
  chrsize_ir <- intervals::Intervals(bin_chrsize)   # chrs intervals
  
  ### overlap the bin and chrs intervals and flag intervals that span two chrs
  ov_idx <- intervals::interval_overlap(bin_ir, chrsize_ir)
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
  gbin_df$start <- with(gbin_df, paste0(start.chr, ":", loc.start))
  gbin_df$end <- with(gbin_df, paste0(end.chr, ":", loc.end))
  
  ## Combine genomic position with HC matrix pos
  gbin_pos <- cbind(gbin_df, hc@POS)
  gbin_pos$gord <- c(1:nrow(gbin_pos))
  gbin_pos_ord <- gbin_pos[order(gbin_pos$y1, decreasing = TRUE),]
  gbin_pos_ord <- gbin_pos_ord[order(gbin_pos_ord$x1, decreasing=FALSE),]
  ## Associate position in matrix (UID) with a genomic position (start, end)
  gbin_pos_ord$uid <- with(gbin_pos_ord, paste(x1, y1, sep="_"))
  gbin_pos_ord <- rbind(gbin_pos_ord, 
                        data.frame('start'='Y:59279092', 'end'='Y:59326329',
                                   'loc.start'=59279092, 'loc.end'=59326329,
                                   'start.chr'='Y', 'end.chr'='Y', 
                                   'x1'=max(gbin_pos_ord$x1)+1, 'y1'=0, 
                                   'x2'=max(gbin_pos_ord$x1)+1, 'y2'=0, 
                                   'gord'=nrow(gbin_pos_ord)+1, 
                                   'uid'=paste0(max(gbin_pos_ord$x1)+1, '_0')))
  ## Reform the HC matrix using UIDs instead of mapping information
  gbin_pos_mat <- matrix(gbin_pos_ord$uid, nrow=max(gbin_pos_ord$y1)+1, ncol=max(gbin_pos_ord$x1)+1)
  return(list("hc"=hc,
              "mat"=gbin_pos_mat, 
              "ord"=gbin_pos_ord))
}

hilbertBinSize <- function(order, scale, chr.size.dat){
  hil_n <- 4^order
  scaled_bp_per_n <- (max(chr.size.dat$cum.end)/scale) / hil_n
  bp_per_n <- scaled_bp_per_n * scale
  return(bp_per_n)
}

plotEuclidDist <- function(xdist, frac_df){
  ht = Heatmap(xdist, name = "dist", 
               cluster_rows = F, cluster_columns = F, 
               show_row_names = F, show_column_names = F, 
               heatmap_legend_param = list(title = "distance"))
  draw(ht, padding = unit(c(5, 5, 5, 2), "mm"))
  
  line_seg <- frac_df$start[-1]
  decorate_heatmap_body("dist", {
    grid.segments(x0=c((1-line_seg), rep(0, length(line_seg))), 
                  x1=c((1-line_seg), rep(1, length(line_seg))), 
                  y0=c(rep(0, length(line_seg)), line_seg),
                  y1=c(rep(1, length(line_seg)), line_seg), 
                  gp = gpar(lty = 2))
    apply(frac_df, 1, function(i, xpos, ypos){
      # side along the y-axis
      i1 <- as.numeric(i[1])
      i2 <- as.numeric(i[2])
      grid.polygon(x=c(xpos, xpos, xpos-0.03), y=c(i2, i1, i2+(i1-i2)/2),
                   gp=gpar(fill=i[3]))
      
      # top along the x-axis
      i1 <- 1-i1
      i2 <- 1-i2
      grid.polygon(x=c(i2, i1, i1+(i2-i1)/2), y=c(ypos, ypos, ypos+0.03), 
                   gp=gpar(fill=i[3]))
    }, xpos=0, ypos=1)
  })
}

#############################
#### Set up HilbertCurve ####
#order <- 8
for(order in c(4,6,8)){
  gbin_pos <- setupRefHcMatrix(order=order)
  gbin_pos$ord$chr <- paste0("chr", gsub(":.*", "", gbin_pos$ord$start))
  gbin_pos$ord$loc.start <- as.integer(gsub("^.*:", "", gbin_pos$ord$start))
  gbin_pos$ord$loc.end <- as.integer(gsub("^.*:", "", gbin_pos$ord$end))
  
  #############################
  #### Set up Genomic Bins ####
  # bin_size <- hilbertBinSize(order, scale, chr.size.dat)
  # chr_sizes <- setNames(end(chr.size.dat), as.character(seqnames(chr.size.dat)))
  # bins_gr   <- tileGenome(chr_sizes, tilewidth=bin_size, 
  #                         cut.last.tile.in.chrom=T)
  #
  bins <- gbin_pos$ord
  bins <- bins[order(bins$gord),]
  uids <- apply(expand.grid(c(0:(max(bins$x1)-1)), 
                            c(0:(max(bins$x1)-1))), 
                1, paste, collapse="_")
  bins$uid <- uids 
  bins$x1 <- as.integer(gsub("_.*", "", uids))
  bins$y1 <- as.integer(gsub("^.*_", "", uids))
  bins$x2 <- c(bins$x1[-1], -1)
  bins$y2 <- c(bins$y1[-1], -1)
  
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
  hc <- genHC(max(chr.size.dat$cum.end)/scale, order)
  hc@POS <- bins[-nrow(bins),c('x1', 'y1', 'x2', 'y2')]
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
