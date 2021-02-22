library(HilbertCurve)
library(HilbertVis)
library(circlize)
library(DescTools)

library(AneuploidyScore)

# Load cytoband/chr-arm data
data("ucsc.hg19.cytoband")
data("seg")

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

setupRefHcMatrix <- function(order=8){
  chr.size.dat <- getChrLength()
  seqlevelsStyle(chr.size.dat) <- 'NCBI'
  
  hc = HilbertCurve(1, max(chr.size.dat$cum.end), 
                    level = order, mode = "pixel", 
                    reference = TRUE, padding=unit(1, "mm"),newpage = FALSE)
  
  ## Identify the intervals for each HC bin by dividing by the zoom factor
  bin_df <- data.frame(start = round(start(hc@BINS) / hc@ZOOM),
                       end = round(end(hc@BINS) / hc@ZOOM))
  
  ## For each interval, identify the corresponding genomic position by mapping cumulative pos
  gbin_df <- apply(bin_df, 1, function(i){
    rng <- sapply(c('start', 'end'), function(id){
      idx <- which(((i[id] + 1) >= chr.size.dat$cum.start) &
                     ((i[id] + 1) <= chr.size.dat$cum.end))
      gpos <- start(chr.size.dat)[idx] + ((i[id]+1) - chr.size.dat$cum.start[idx])
      seq <- as.character(seqnames(chr.size.dat))[idx]
      paste(c(seq, gpos), collapse=":")
    })
    rng
  })
  gbin_df <- as.data.frame(t(gbin_df))
  
  ## Combine genomic position with HC matrix pos
  gbin_pos <- cbind(gbin_df, hc@POS)
  gbin_pos_ord <- gbin_pos[order(gbin_pos$y1, decreasing = TRUE),]
  gbin_pos_ord <- gbin_pos_ord[order(gbin_pos_ord$x1, decreasing=FALSE),]
  ## Associate position in matrix (UID) with a genomic position (start, end)
  gbin_pos_ord$uid <- with(gbin_pos_ord, paste(x1, y1, sep="_"))
  gbin_pos_ord <- rbind(gbin_pos_ord, 
                        data.frame('start'='Y:59279092', 'end'='Y:59326329',
                                   'x1'=255, 'y1'=0, 'x2'=256, 'y2'=0, 'uid'='255_0'))
  ## Reform the HC matrix using UIDs instead of mapping information
  gbin_pos_mat <- matrix(gbin_pos_ord$uid, nrow=max(gbin_pos_ord$y1)+1, ncol=max(gbin_pos_ord$x1)+1)
  return(list("mat"=gbin_pos_mat, "ord"=gbin_pos_ord))
}


##############
#### Main ####
chr.size.dat <- getChrLength()
seqlevelsStyle(chr.size.dat) <- 'NCBI'
scale <- 1e4
order=8
a_col=c("black", "red")
b_col=c("black", "blue")
breaks <- seq(0, 5, by=0.1)
#col2rgb(colorRampPalette(a_col)(length(breaks)))
#col2rgb(colorRampPalette(b_col)(length(breaks)))



## Add cumulative loci position
seqlevelsStyle(seg) <- 'ENSEMBL'
chr.data <- addCumPos(seg, chr.size.dat, dat.type='seg')
ir <- IRanges(chr.data$cloc.start/scale, chr.data$cloc.end/scale)

## Fit Hilbert Curve and colorize
hc = HilbertCurve(1, max(chr.size.dat$cum.end)/scale, level = order, mode = "pixel", 
                  reference = TRUE, padding=unit(1, "mm")) # title = basename(segf), 
a1 = circlize::colorRamp2(breaks = breaks,
                          colors = colorRampPalette(a_col)(length(breaks)))
a2 = circlize::colorRamp2(breaks = breaks,
                          colors = colorRampPalette(b_col)(length(breaks)))
#cols <- DescTools::MixColor(col1 = a1(breaks), col2 = a2(breaks))
# Ensure that Red/Blue channels range from 0-255. MixColor halves both values when mixing
cols <- DescTools::MixColor(col1 = a1(chr.data$nMajor), col2 = a2(chr.data$nMinor))
cols <- apply(col2rgb(cols),2,function(x) rgb(x[1]*2, x[2]*2, x[3]*2, maxColorValue = 255))

hc_layer(hc, ir, mean_mode = "absolute", col = cols)
if(annotate){
  hc_polygon(hc, x1=chr.size.dat$cum.start/scale, x2=chr.size.dat$cum.end/scale)
  hc_text(hc, x1=chr.size.dat$cum.start/scale, x2=chr.size.dat$cum.end/scale, 
          labels = c(1:22, "X", "Y"), gp = gpar(fontsize = 10), 
          centered_by = "polygon")
}
dev.off()

gbin_pos <- setupRefHcMatrix(order=8)
gbin_pos$ord$chr <- paste0("chr", gsub(":.*", "", gbin_pos$ord$start))
gbin_pos$ord$loc.start <- as.integer(gsub("^.*:", "", gbin_pos$ord$start))
gbin_pos$ord$loc.end <- as.integer(gsub("^.*:", "", gbin_pos$ord$end))

set.seed(123)
xindices <- sample(1:256, size=1000, replace = T) # c(1:nrow(gbin_pos$mat))
yindices <- sample(1:256, size=1000, replace = T) # c(1:ncol(gbin_pos$mat))
xy_indices <- data.frame(xindices, yindices) # expand.grid(xindices, yindices)
dmat <- apply(xy_indices, 1, function(xy){
  xidx <- xy[1]
  yidx <- xy[2]
  
  # In the matrix, get the reference UIDs for xidx and yidx
  ref <- gbin_pos$mat[xidx, yidx]
  
  # Get all the immediately surrounding UIDs to the ref xidx and yidx
  xalt <- xidx + c(-1,0,1)
  yalt <- yidx + c(-1,0,1)
  altidx <- expand.grid(xalt, yalt)[-5,]
  oob <- (altidx < 1 | altidx > 255)
  if(any(oob)) altidx <- altidx[-which(as.logical(rowSums(oob))),]
  alt <- apply(altidx, 1, function(i) {gbin_pos$mat[i[1], i[2]]})
  
  # Obtain the corresponding genomic positions for reference and comp UIDs
  ref_df <- gbin_pos$ord[match(ref, gbin_pos$ord$uid),]
  alt_df <- gbin_pos$ord[match(alt, gbin_pos$ord$uid),]
  
  # Get the genomic distance between the locis
  cols <- c('chr', 'loc.start', 'loc.end', 'uid')
  comp_bed <- suppressWarnings(cbind(ref_df[,cols],
                                     alt_df[,cols]))
  colnames(comp_bed) <- apply(expand.grid(cols, c("A", "B")), 1, paste, collapse="_")
  comp_bed$d <- abs(with(comp_bed, loc.end_A - loc.end_B))
  diff_chr <- with(comp_bed, chr_A != chr_B)
  if(any(diff_chr)) comp_bed$d[which(diff_chr)] <- -1 #-1 if differnet chromosomes
  
  return(comp_bed)
})
dmat <- as.data.frame(do.call(rbind, dmat))



##CIRCOS DEMO
library(circlize)
dist_threshold <- 50*1e6
gap_dmat <- dmat[which((dmat$d > dist_threshold) | (dmat$d == -1)),]

bed1 <- gap_dmat[,grep("_A$", colnames(dmat))]
bed2 <- gap_dmat[,grep("_B$", colnames(dmat))]
colnames(bed1) <- colnames(bed2) <- c("chr", "start", "end", "UID")
circos.initializeWithIdeogram()
circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.4), 
                   border = TRUE)




