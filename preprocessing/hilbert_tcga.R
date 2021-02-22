library(HilbertCurve)
library(HilbertVis)
library(circlize)
library(DescTools)

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

analysis <- 'ccl_aggregate' # 'TCGA', 'ccl_aggregate'
annotate <- FALSE
PDIR <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines'
if(analysis=='TCGA'){
  seg_files <- "TCGA_mastercalls.abs_segtabs.fixed.txt"
} else if(analysis=='ccl_aggregate'){
  seg_files <- c('CCLE_cna_hg19.seg', 'GDSC_cna_hg19.seg', 'gCSI_cna_hg19.seg')
} else {
  stop("Analysis must be either 'TCGA' or 'ccl_aggregate'")
}

for(seg_i in seg_files){
  segf <- file.path(PDIR, analysis, "input", seg_i)
  segd <- read.table(segf, sep="\t", header=TRUE, stringsAsFactors = FALSE)
  if(analysis=='TCGA'){
    colnames(segd)[2:4] <- c('chrom', 'loc.start', 'loc.end')
    
    normal_segf <- file.path(PDIR, '1000G', 'eacon', 'symlinks', 'seg', '1000g.seg')
    normal_segd <- read.table(normal_segf, sep="\t", header=TRUE, stringsAsFactors = FALSE)
    colnames(normal_segd) <- c('Sample', 'chrom', 'Chrom', 'loc.start', 'loc.end', 
                               'Length', 'Modal_Total_CN', 'Modal_HSCN_1', 'Modal_HSCN_2')
    segd2 <- plyr::rbind.fill(segd, normal_segd)
    #segd <- segd2[grep("^HG", segd2$Sample),]
    segd <- segd2
  } else {
    colnames(segd)[c(1:3,6,11,7,8)] <- c('chrom', 'loc.start', 'loc.end',
                                  'Sample', 'Modal_Total_CN', 
                                  'Modal_HSCN_1', 'Modal_HSCN_2')
    if(any(segd$chrom %in% c('chrX', 'chrY'))){
      segd <- segd[-which(segd$chrom %in% c('chrX', 'chrY')),]
    }
  }
  
  
  segl <- split(segd, segd$Sample)
  
  hilberts <- lapply(segl, function(seg){
    sample_id <- unique(seg$Sample) 
    print(paste0(seg_i, "-Sample: ", sample_id, " (", 
                 grep(sample_id, names(segl)), "/", length(segl), ")"))
    
    if(!all(is.na(seg$chrom))){
      seg$chrom <- gsub("23", "X", seg$chrom) 
      seg$chrom <- gsub("24", "Y", seg$chrom) 
      seg$chrom <- gsub("^chr", "", seg$chrom, ignore.case = TRUE)
      chr.data <- addCumPos(seg, chr.size.dat, dat.type='seg')
      ir <- IRanges(chr.data$cloc.start/scale, chr.data$cloc.end/scale)
      
      ## Plot Allele specific copy-number
      png(file.path(PDIR, analysis, "output", "hilbert", paste0(sample_id, ".png")), 
          width=300, height=300)
      hc = HilbertCurve(1, max(chr.size.dat$cum.end)/scale, level = order, mode = "pixel", 
                        reference = TRUE, padding=unit(1, "mm")) # title = basename(segf), 
      a1 = circlize::colorRamp2(breaks = breaks,
                                colors = colorRampPalette(a_col)(length(breaks)))
      a2 = circlize::colorRamp2(breaks = breaks,
                                colors = colorRampPalette(b_col)(length(breaks)))
      #cols <- DescTools::MixColor(col1 = a1(breaks), col2 = a2(breaks))
      # Ensure that Red/Blue channels range from 0-255. MixColor halves both values when mixing
      cols <- DescTools::MixColor(col1 = a1(chr.data$Modal_HSCN_1), col2 = a2(chr.data$Modal_HSCN_2))
      cols <- apply(col2rgb(cols),2,function(x) rgb(x[1]*2, x[2]*2, x[3]*2, maxColorValue = 255))
      
      hc_layer(hc, ir, mean_mode = "absolute", col = cols)
      if(annotate){
        hc_polygon(hc, x1=chr.size.dat$cum.start/scale, x2=chr.size.dat$cum.end/scale, col='white')
        hc_text(hc, x1=chr.size.dat$cum.start/scale, x2=chr.size.dat$cum.end/scale, 
                labels = c(1:22, "X", "Y"), gp = gpar(fontsize = 10), 
                centered_by = "polygon")
      }
      dev.off()
      
    }
  })
}






