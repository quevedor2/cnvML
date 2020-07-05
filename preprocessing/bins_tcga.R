library(dplyr)
library(GenomicRanges)

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

formatSeg <- function(segd, analysis, PDIR=NULL){
  if(analysis=='TCGA'){
    colnames(segd)[2:4] <- c('chrom', 'loc.start', 'loc.end')
    
    normal_segf <- file.path(PDIR, '1000G', 'eacon', 'symlinks', 'seg', '1000g.seg')
    normal_segd <- read.table(normal_segf, sep="\t", header=TRUE, stringsAsFactors = FALSE)
    colnames(normal_segd) <- c('Sample', 'chrom', 'Chrom', 'loc.start', 'loc.end', 
                               'Length', 'Modal_Total_CN', 'Modal_HSCN_1', 'Modal_HSCN_2')
    segd2 <- plyr::rbind.fill(segd, normal_segd)
    segd <- segd2
  } else {
    colnames(segd)[c(1:3,6,11,7,8)] <- c('chrom', 'loc.start', 'loc.end',
                                         'Sample', 'Modal_Total_CN', 
                                         'Modal_HSCN_1', 'Modal_HSCN_2')
    if(any(segd$chrom %in% c('chrX', 'chrY'))){
      segd <- segd[-which(segd$chrom %in% c('chrX', 'chrY')),]
    }
  }
  
  if(sum(grepl("chrom", colnames(segd), ignore.case=TRUE)) > 1){
    segd <- segd[,-grep("Chrom", colnames(segd))]
  }
  if(any(is.na(segd$chrom))) segd <- segd[-which(is.na(segd$chrom)),]
  segd$chrom <- gsub("23", "X", segd$chrom) %>% gsub("24", "Y", .)
  return(segd)
}

hilbertBinSize <- function(order, scale, chr.size.dat){
  hil_n <- 4^order
  scaled_bp_per_n <- (max(chr.size.dat$cum.end)/scale) / hil_n
  bp_per_n <- scaled_bp_per_n * scale
  return(bp_per_n)
}

##############
#### Main ####
## Load in Metadata
METADIR='/mnt/work1/users/pughlab/projects/cancer_cell_lines/TCGA'
META <- file.path(METADIR, 'input/merged_sample_quality_annotations.trimmed.tsv')
meta <- read.table(META, header=TRUE, check.names = FALSE, 
                   sep="\t", stringsAsFactors = FALSE)
m_idx <- NULL

## Load in SEG files
analysis <- 'TCGA'
PDIR <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines'
if(analysis=='TCGA'){
  seg_files <- "TCGA_mastercalls.abs_segtabs.fixed.txt"
} else if(analysis=='ccl_aggregate'){
  seg_files <- c('CCLE_cna_hg19.seg', 'GDSC_cna_hg19.seg', 'gCSI_cna_hg19.seg')
} else {
  stop("Analysis must be either 'TCGA' or 'ccl_aggregate'")
}
OUTDIR <- file.path(PDIR, analysis, "output", "bin")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

chr.size.dat <- getChrLength()
seqlevelsStyle(chr.size.dat) <- 'NCBI'
scale <- 1e5
order=7

for(seg_i in seg_files){
  segf <- file.path(PDIR, analysis, "input", seg_i)
  segd <- read.table(segf, sep="\t", header=TRUE, stringsAsFactors = FALSE)
  segd <- formatSeg(segd, analysis, PDIR=PDIR)
  seg_gr <- GenomicRanges::makeGRangesFromDataFrame(segd, keep.extra.columns = TRUE)
  
  bin_size <- hilbertBinSize(order, scale, chr.size.dat)
  chr_sizes <- setNames(end(chr.size.dat), as.character(seqnames(chr.size.dat)))
  bins_gr   <- tileGenome(chr_sizes, tilewidth=bin_size, 
                          cut.last.tile.in.chrom=T)
  seqlevelsStyle(bins_gr) <- seqlevelsStyle(seg_gr) <- 'UCSC'
  bins_l <- split(bins_gr, seqnames(bins_gr))
  bins_gr$ID <- as.character(unlist(sapply(bins_l, function(b){
    paste(as.character(seqnames(b)), c(1:length(b)), sep="_")
  })))
  
  blank_l2r <- setNames(rep(NA, length(bins_gr)), bins_gr$ID)
  seg_l <- split(seg_gr, seg_gr$Sample)
  modal <- paste0("Modal_", c('Total_CN', 'HSCN_1')) #'HSCN_2', 
  for(m in modal){
    print(paste0(m, "..."))
    l2r_mat <- sapply(seg_l, function(seg){
      seg_l2r <- blank_l2r
      ov <- findOverlaps(seg, bins_gr)
      seg_l2r[subjectHits(ov)] <- GenomicRanges::mcols(seg[queryHits(ov),])[,m]
      seg_l2r
    })
    l2r_mat <- as.data.frame(t(l2r_mat))
    
    if(is.null(m_idx)){
      m_idx <- sapply(rownames(l2r_mat), function(i) grep(paste0("^", i), 
                                                          x=meta$aliquot_barcode)[1])
    }
    l2r_mat$cancer_type <- meta[m_idx,]$`cancer type`
    l2r_mat$cancer_type[grep("^HG0", rownames(l2r_mat))] <- 'Normal'

    write.table(l2r_mat, file=file.path(OUTDIR, paste0(m, ".csv")), sep=",",
                col.names=TRUE, row.names=TRUE, quote=FALSE)
  }
}
