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
PDIR <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/TCGA'
segf <- file.path(PDIR, "input", "TCGA_mastercalls.abs_segtabs.fixed.txt")
segd <- read.table(segf, sep="\t", header=TRUE, stringsAsFactors = FALSE)
colnames(segd)[2:4] <- c('chrom', 'loc.start', 'loc.end')

chr.size.dat <- getChrLength()
seqlevelsStyle(chr.size.dat) <- 'NCBI'
scale <- 1e5
order=7

segl <- split(segd, segd$Sample)
hilberts <- lapply(segl[-c(1:180)], function(seg){
  sample_id <- unique(seg$Sample) 
  print(paste0("Sample: ", sample_id, " (", grep(sample_id, names(segl)), 
               "/", length(segl), ")"))
  
  if(!all(is.na(seg$chrom))){
    seg$chrom <- gsub("23", "X", seg$chrom) 
    seg$chrom <- gsub("24", "Y", seg$chrom) 
    chr.data <- addCumPos(seg, chr.size.dat, dat.type='seg')
    ir <- IRanges(chr.data$cloc.start/scale, chr.data$cloc.end/scale)
    
    ## Plot Allele specific copy-number
    png(file.path(PDIR, "output", "hilbert", paste0(sample_id, ".png")), 
        width=300, height=300)
    hc = HilbertCurve(1, max(chr.size.dat$cum.end)/scale, level = order, mode = "pixel", 
                      reference = TRUE) # title = basename(segf), 
    breaks <- seq(0, 5, by=0.1)
    a1 = circlize::colorRamp2(breaks = breaks,
                              colors = colorRampPalette(c("white", "red", "darkred"))(length(breaks)))
    a2 = circlize::colorRamp2(breaks = breaks,
                              colors = colorRampPalette(c("white", "blue", "darkblue"))(length(breaks)))
    cols <- DescTools::MixColor(col1 = a1(chr.data$Modal_HSCN_1), col2 = a2(chr.data$Modal_HSCN_2))
    hc_layer(hc, ir, mean_mode = "absolute", col = cols)
    hc_polygon(hc, x1=chr.size.dat$cum.start/scale, x2=chr.size.dat$cum.end/scale)
    hc_text(hc, x1=chr.size.dat$cum.start/scale, x2=chr.size.dat$cum.end/scale, 
            labels = c(1:22, "X", "Y"), gp = gpar(fontsize = 10), 
            centered_by = "polygon")
    dev.off()
    
  }
})

