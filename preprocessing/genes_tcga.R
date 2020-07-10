###################
#### Functions ####
###################
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

getMapping <- function(in.col='ENTREZID', 
                       out.cols=c("SYMBOL", "ENSEMBL")){
  suppressPackageStartupMessages(require(org.Hs.eg.db))
  gene.map <- as.data.frame(sapply(out.cols, function(oc){
    mapIds(org.Hs.eg.db, keys=keys(org.Hs.eg.db, in.col),
           keytype="ENTREZID", column=oc, multiVals = 'first')
  }))
  gene.map$ENTREZID <- keys(org.Hs.eg.db, in.col)
  # gene.map <- select(org.Hs.eg.db, keys=keys(org.Hs.eg.db, in.col), 
  #                    keytype="ENTREZID", columns=out.cols)
  
  gene.map[,1] <- as.character(gene.map[,1])
  gene.map[,2] <- as.character(gene.map[,2])
  gene.map[,3] <- as.character(gene.map[,3])
  gene.map
}

assignEntrezToSegment <- function(cnv0, anno){
  olaps = findOverlaps(cnv0, anno)
  mcols(olaps)$gene_id = anno$gene_id[subjectHits(olaps)]  # Fixed the code here
  cnv_factor = factor(queryHits(olaps), levels=seq_len(queryLength(olaps)))
  gene.id <- IRanges::splitAsList(mcols(olaps)$gene_id, cnv_factor)
  return(gene.id)
}

splitSegmentByGene <- function(cnv0, cols){
  seg.entrez <- apply(as.data.frame(mcols(cnv0)), 1, function(i){
    ids <- unlist(strsplit(x = as.character(unlist(i[['gene_id']])), split=","))
    segs <- do.call(rbind, replicate(length(ids), round(unlist(i[cols]),3), simplify = FALSE))
    
    as.data.frame(cbind(segs, 'ENTREZ'=ids))
  })
  seg.entrez <- do.call(rbind, seg.entrez)
  if(any(duplicated(seg.entrez$ENTREZ))) seg.entrez <- seg.entrez[-which(duplicated(seg.entrez$ENTREZ)),]
  return(seg.entrez)
}

annotateCNVs <- function(cnv, txdb, anno=NULL,
                         cols=c("seg.mean", "nA", "nB"),
                         l2r.dat=NULL){
  stopifnot(is(cnv, "GRanges"), is(txdb, "TxDb"))
  
  ## Assign EntrezID to each segment 
  if(is.null(anno)) anno = genes(txdb)
  cnv$gene_id = assignEntrezToSegment(cnv, anno)
  if(!is.null(l2r.dat)) l2r.dat$gene_id <-  assignEntrezToSegment(l2r.dat, anno)
  
  ## Split the segments by genes
  seg.entrez <- splitSegmentByGene(cnv, cols)
  if(!is.null(l2r.dat)) {
    l2r.entrez <- splitSegmentByGene(l2r.dat, cols='Log2Ratio')
    seg.l2r.entrez <- merge(seg.entrez, l2r.entrez, by='ENTREZ', all.x=TRUE)
    seg.entrez <- seg.l2r.entrez
  }
  
  ## Map ensembl and HUGO IDs to the ENTREZ ids
  seg.anno <- merge(seg.entrez, getMapping(),
                    by.x="ENTREZ", by.y="ENTREZID", all.x=TRUE)
  if(any(duplicated(seg.anno$ENTREZ))) seg.anno <- seg.anno[-which(duplicated(seg.anno$ENTREZ)),]
  cols <- colnames(seg.entrez)[grep("ENTREZ", colnames(seg.entrez), invert = TRUE)]
  for(each.col in cols){
    seg.anno[,each.col] <- as.numeric(as.character(seg.anno[,each.col]))
  }
  
  list("seg"=cnv, "genes"=seg.anno)  
}

getGenes <- function(genome.build="hg19", make.into.gr=FALSE){
  switch(genome.build,
         hg18={
           suppressPackageStartupMessages(require(TxDb.Hsapiens.UCSC.hg18.knownGene))
           if(!exists("TxDb.Hsapiens.UCSC.hg18.knownGene")) stop("Requires TxDb.Hsapiens.UCSC.hg18.knownGene")
           package <- TxDb.Hsapiens.UCSC.hg18.knownGene
         },
         hg19={
           suppressPackageStartupMessages(require(TxDb.Hsapiens.UCSC.hg19.knownGene))
           if(!exists("TxDb.Hsapiens.UCSC.hg19.knownGene")) stop("Requires TxDb.Hsapiens.UCSC.hg19.knownGene")
           package <- TxDb.Hsapiens.UCSC.hg19.knownGene
         },
         hg38={
           suppressPackageStartupMessages(require(TxDb.Hsapiens.UCSC.hg38.knownGene))
           if(!exists("TxDb.Hsapiens.UCSC.hg38.knownGene")) stop("Requires TxDb.Hsapiens.UCSC.hg38.knownGene")
           package <- TxDb.Hsapiens.UCSC.hg38.knownGene
         },
         stop("genome must be 'hg19' or 'hg38'"))
  
  if(make.into.gr){
    genes0 <- genes(package)
    idx <- rep(seq_along(genes0), elementNROWS(genes0$gene_id))
    genes <- granges(genes0)[idx]
    genes$gene_id = unlist(genes0$gene_id)
    genes
  } else {
    list(txdb=package,
         txdb.genes=genes(package))
  }
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

##############
#### MAIN ####
##############
library(parallel)
library(dplyr)

##################################
#### 1) Load in data and meta ####
analysis <- 'TCGA' # 'TCGA' or 'ccl_aggregate'
PDIR <- file.path('/mnt/work1/users/pughlab/projects/cancer_cell_lines', analysis)
fl_id <- switch(analysis,
                ccl_aggregate='all_cna_hg19.seg',
                TCGA='TCGA_mastercalls.abs_segtabs.fixed.txt')
segf <- file.path(PDIR, "input", fl_id)
segd <- read.table(segf, sep="\t", header=TRUE, stringsAsFactors = FALSE)

if(analysis == 'ccl_aggregate'){
  META <- file.path(PDIR, "ref", "onco_meta_df.rds")
  anno <- readRDS(META)
  anno$GDSC <- gsub(".cel$", "", anno$GDSC, ignore.case = TRUE)
  anno$aliquot_barcode <- apply(anno[,c('CCLE', 'GNE', 'GDSC')], 1, paste, collapse=",")
  colnames(anno)[grep("oncocode", colnames(anno))] <- 'cancer type'
} else if(analysis=='TCGA'){
  META <- file.path(PDIR, 'input/merged_sample_quality_annotations.trimmed.tsv')
  anno <- read.table(META, header=TRUE, check.names = FALSE, 
                     sep="\t", stringsAsFactors = FALSE)
}
segd <- formatSeg(segd, analysis, PDIR=dirname(PDIR))

#####################################
#### 2) Annotate Sample Segments ####
## Splits the master seg into individual seg files per sample
chr.size.dat <- getChrLength()
seqlevelsStyle(chr.size.dat) <- 'NCBI'
genes <- getGenes('hg19')
segl <- split(segd, f=segd$Sample)
chrdat <- sapply(segl, function(i) unique(i$chr))
if(class(chrdat) == 'list') segl <- segl[which(sapply(chrdat, length) >= 22)]

## Creates a dataframe of CN (TCN, HSCN_1/2) by gene per sample
genes_l <- lapply(segl, function(cnv) {
  print(paste0("Sample: ", unique(cnv$Sample), " (",
               grep(unique(cnv$Sample)[1], names(segl)), "/", length(segl), ")"))
  cnv <- makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
  seqlevels(cnv) <- as.character(unique(seqnames(cnv)))
  seqlevelsStyle(cnv) <- 'UCSC'
  annotateCNVs(cnv, genes$txdb, anno=genes$txdb.genes, 
               cols=c('Modal_HSCN_1', 'Modal_HSCN_2', 'Modal_Total_CN'))$genes
})
dir.create(file.path(PDIR, "output", "gene_matrix"), showWarnings = FALSE, recursive = TRUE)
save(genes_l, file=file.path(PDIR, "output", "gene_matrix", "genes_parallel.rda"))

##########################
#### 3) Create Matrix ####
load(file.path(PDIR, "output", "gene_matrix", "genes_parallel.rda"))

mapping <- getMapping()[,'SYMBOL', drop=FALSE]
colids <- c("Modal_HSCN_1", "Modal_HSCN_2", "Modal_Total_CN")
modal_genes <- lapply(setNames(colids, colids), function(g){
  gl_mat <- mapping
  
  for(sample in names(genes_l)){
    if((grep(sample, names(genes_l)) %% 100) == 0){
      print(paste0(sample, ": (", grep(sample, names(genes_l)), "/", length(genes_l), ")"))
    }
    gl <- genes_l[[sample]]
    midx <- match(mapping$SYMBOL, gl$SYMBOL)
    gl_mat <- cbind(gl_mat, gl[midx,g, drop=FALSE])
    colnames(gl_mat)[ncol(gl_mat)] <- sample
  }
  
  return(gl_mat)
})
saveRDS(modal_genes, file=file.path(PDIR, "output", "gene_matrix", "modal_genes.rds"))
#modal_genes <- readRDS(file.path(PDIR, "output", "gene_matrix", "modal_genes.rds"))

###############################
#### 4) Assign cancer type ####
allna <- lapply(modal_genes, function(i) apply(is.na(i[,-1]), 1, all))
narows <- unique(sort(c(sapply(allna, which), which(duplicated(modal_genes[[1]]$SYMBOL)))))
ctypes <- sapply(colnames(modal_genes[[1]]), function(sample){
  as.character(anno[grep(sample, anno$aliquot_barcode)[1],]$'cancer type')
})

if(analysis=='TCGA'){
  modal_genes2 <- lapply(modal_genes, function(mg){
    mg <- mg[-narows,]
    mg <- rbind(mg, t(ctypes))
    rownames(mg) <- c(mg$SYMBOL[-nrow(mg)], "cancer_type")
    mg <- t(mg[,-1])
    mg
  })
} else if(analysis=='ccl_aggregate'){
  ## Read in existing TCGA genes and extract the gene name and order
  modal_TCGA <-readRDS(file.path(gsub("ccl_aggregate", "TCGA", PDIR), 
                                 "output", "gene_matrix", "modal_genes2.rds"))
  genes_TCGA <- colnames(modal_TCGA[[1]])
  
  ## Subset the ccl_aggregate genes for just the TCCGA genes
  modal_genes2 <- lapply(modal_genes, function(mg){
    mg <- rbind(mg, t(ctypes))
    if(any(duplicated(mg$SYMBOL))) mg <- mg[-which(duplicated(mg$SYMBOL)),]
    rownames(mg) <- c(mg$SYMBOL[-nrow(mg)], "cancer_type")
    mg <- t(mg[,-1])
    
    mg[,match(genes_TCGA, colnames(mg))]
  })
}
saveRDS(modal_genes2, file=file.path(PDIR, "output", "gene_matrix", "modal_genes2.rds"))
#modal_genes2 <- readRDS(file.path(PDIR, "output", "gene_matrix", "modal_genes2.rds"))

lapply(names(modal_genes2), function(data_type){
  write.csv(modal_genes2[[data_type]], 
            file=file.path(PDIR, "output", "gene_matrix", paste0(data_type, "_matrix.csv")),
            quote = FALSE, row.names = TRUE, col.names = TRUE)
})

