suppressPackageStartupMessages(library(optparse))

####################
#### Parameters ####
option_list <- list( 
  make_option(c("-p", "--pdir"), type="character", 
              default='/mnt/work1/users/pughlab/projects/cancer_cell_lines',
              help="Path to seg files [dir]/[dataset]/input/[seg_file] [%default]"),
  make_option(c("-d", "--dataset"), type="character", default='ccl_aggregate',
              help="Dataset to use, 'ccl_aggregate' or 'TCGA' [%default]"),
  make_option(c("-s", "--segfile"), type="character", default=NULL,
              help="If --dataset is set to custom, indicate the name of the seg_file [%default]"),
  make_option(c("-m", "--metafile"), type="character", default=NULL,
              help="If --dataset is set to custom, indicate the meta RDS file path [%default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

########################
#### Load Libraries ####
library(dplyr)
library(RcnvML)
library(GenomicRanges)

###################
#### Functions ####
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
## Load in SEG files
analysis <- opt$dataset 
m_idx <- NULL
PDIR <- opt$pdir
if(analysis=='TCGA'){
  seg_files <- "TCGA_mastercalls.abs_segtabs.fixed.txt"
  METADIR='/mnt/work1/users/pughlab/projects/cancer_cell_lines/TCGA'
  META <- file.path(METADIR, 'input/merged_sample_quality_annotations.trimmed.tsv')
  meta <- read.table(META, header=TRUE, check.names = FALSE, 
                     sep="\t", stringsAsFactors = FALSE)
} else if(analysis=='ccl_aggregate'){
  seg_files <- c('CCLE_cna_hg19.seg', 'GDSC_cna_hg19.seg', 'gCSI_cna_hg19.seg')
  META <- file.path(PDIR, analysis, "ref", "onco_meta_df.rds")
  meta <- readRDS(META)
  meta$GDSC <- gsub(".cel$", "", meta$GDSC, ignore.case = TRUE)
  colnames(meta) <- gsub("oncocode", "cancer type", colnames(meta))
  colnames(meta) <- gsub("^GNE$", "gCSI", colnames(meta))
} else if(analysis=='custom'){
  seg_files <- opt$segfile
  meta <- readRDS(opt$metafile)
} else {
  stop("Analysis must be either 'TCGA' or 'ccl_aggregate' or 'custom'")
}
OUTDIR <- file.path(PDIR, analysis, "output", "bin")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

chr.size.dat <- getChrLength()
seqlevelsStyle(chr.size.dat) <- 'NCBI'
scale <- 1e4
order=8

for(seg_i in seg_files){
  ds_id <- gsub("_.*", "", seg_i)
  print(paste0(ds_id, "..."))
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
  modal <- paste0("Modal_", c('Total_CN', 'HSCN_2', 'HSCN_1')) # 
  for(m in modal){
    print(paste0(m, "..."))
    l2r_mat <- sapply(seg_l, function(seg){
      seg_l2r <- blank_l2r
      ov <- findOverlaps(seg, bins_gr)
      seg_l2r[subjectHits(ov)] <- GenomicRanges::mcols(seg[queryHits(ov),])[,m]
      seg_l2r
    })
    l2r_mat <- as.data.frame(t(l2r_mat))
    
    ds <- switch(analysis,
                 TCGA='aliquot_barcode',
                 ccl_aggregate=gsub("_.*", "", seg_i))
    m_idx <- sapply(rownames(l2r_mat), 
                    function(i) grep(paste0("^", i), 
                                     x=meta[,ds])[1])
    l2r_mat$cancer_type <- gsub(":.*", "", as.character(meta[m_idx,]$`cancer type`))
    l2r_mat$cancer_type[grep("^HG0", rownames(l2r_mat))] <- 'Normal'

    dir.create(file.path(OUTDIR, ds_id), recursive = TRUE, showWarnings = FALSE)
    write.table(l2r_mat, file=file.path(OUTDIR, ds_id, paste0(m, "_matrix.csv")), 
                sep=",", col.names=TRUE, row.names=TRUE, quote=FALSE)
  }
  saveRDS(bins_gr, file.path(OUTDIR, "bins_ref.rds"))
}
