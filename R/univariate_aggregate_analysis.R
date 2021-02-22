PDIR <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnvML/CCL/'
###################
#### FUNCTIONS ####
getPharmacoIDs <- function(meta_file, dataset, Xid, reverse=FALSE){
  meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
  if(reverse){
      search_id=dataset
  } else {
    search_id='PharmacoGX_ID'
  }
  pharmacoID <- meta[match(Xid, meta[,dataset]), search_id]
  
  return(pharmacoID)
}

aggregateSigUni <- function(analysis, threshold=0.01, padj_method='bonferroni'){
  all_drugs <- lapply(all_ds, function(dataset){
    OUTDIR <- file.path(PDIR, "output", "univariate", analysis)
    drug <- readRDS(file=file.path(OUTDIR, paste0(dataset, "_", analysis, ".rds")))
    return(drug)
  })
  
  drug_sig_list <- lapply(setNames(drugs_oi,drugs_oi), function(d){
    ## Metap method: fishers sumlog
    ds_aov <- lapply(all_drugs, function(ds){
      data.frame("gene"=rownames(ds[[d]]),
                 "AOV"=ds[[d]]$AOV)
    })
    ds_aov <- Reduce(function(x,y) merge(x,y,by='gene'), ds_aov)
    pval <- setNames(apply(ds_aov[,-1], 1, function(i) metap::sumlog(i)$p),
                     ds_aov[,1])
    qval <- p.adjust(pval, method=padj_method)
    return(qval[qval < 0.05])
  })
  return(drug_sig_list)
}

matchSigToRef <- function(drug_thresh_list, ref, col_id){
  sig_pos <- lapply(drug_thresh_list, function(sig){
    match_sig <- ref[match(names(sig), mcols(ref)[,col_id]),]
    match_sig$q <- sig
    return(match_sig)
  })
  return(sig_pos)
}

getMeanAov <- function(analysis){
  all_drugs <- lapply(all_ds, function(dataset){
    OUTDIR <- file.path(PDIR, "output", "univariate", analysis)
    drug <- readRDS(file=file.path(OUTDIR, paste0(dataset, "_", analysis, ".rds")))
    return(drug)
  })
  
  drugs_aov <- lapply(setNames(drugs_oi, drugs_oi), function(drug){
    allds_drug <- plyr::rbind.fill(lapply(all_drugs, function(drug_ds){
      as.data.frame(t(drug_ds[[drug]][,'AOV', drop=FALSE]))
    }))
    return(as.data.frame(t(allds_drug)))
  })
  drugs_aov_mu <- lapply(drugs_aov, rowMeans, na.rm=TRUE)
  return(drugs_aov_mu)
}

combineThresh <- function(dat){
  x <- dat[['0.15']]
  x <- lapply(x, function(i) {i$thresh <- '0.15'; return(i)})
  for(thr in rev(names(dat))){
    for(drug in drugs_oi){
      ov_idx <- findOverlaps(dat[[thr]][[drug]],
                             x[[drug]])
      x[[drug]]$thresh[subjectHits(ov_idx)] <- thr
    }
  }
  return(x)
}

###############
#### SETUP ####
## PATHS ##
gm_file <- file.path(PDIR, "input", "gene", "matrix", "Modal_Total_CN_matrix.csv")
bm_file <- file.path(PDIR, "input", "bin", "matrix", "Modal_Total_CN_matrix.csv")
pset_dir <- file.path(PDIR, "input", "drugs", "PSets")
meta_file <- file.path(PDIR, "input", "meta_df.csv")

## VARIABLES ##
analysis <- 'gene'
meta_df <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
all_ds <- setNames(c("CCLE", "GDSC", "GNE"), c("CCLE", "GDSC", "GNE"))
ds <- 'CCLE'
ds_map <- setNames(c("CTRPv2", "GDSC2", "gCSI"), names(all_ds))
patterns <- list("CCLE"="_GenomeWideSNP_", "GNE"="^Unk[0-9]+")
drugs_oi <- c("Docetaxel", "Entinostat", "Erlotinib",
              "Gemcitabine", "Lapatinib", "Paclitaxel", "Pictilisib", "Sirolimus",
              "Vorinostat")

## Read in Data ##
psets <- lapply(all_ds, function(ds){
  read.csv(file.path(pset_dir, paste0(ds_map[ds], "_aac.csv")), 
           check.names = FALSE, stringsAsFactors = FALSE)
})

m_file <- if(analysis == 'gene') gm_file else bm_file
gm <- read.csv(m_file, check.names = FALSE, stringsAsFactors = FALSE)
if(analysis == 'gene') {
  gm <- gm[-which(duplicated(gm[,1])),]
  rownames(gm) <- gm[,1]
  gm <- gm[,-1]
}

gm_split_ord <- rep('GDSC', nrow(gm))
for(id in names(patterns)){ gm_split_ord[grep(patterns[[id]], rownames(gm))] <- id }
gm_ds <- split(gm, gm_split_ord)

## Map IDs to PharmacoGX IDs
ds_pharma_ids <- lapply(all_ds, function(dataset){
  pharma_ids <- getPharmacoIDs(meta_file=meta_file, Xid=rownames(gm_ds[[dataset]]),
                               dataset=dataset, reverse=FALSE)
  return(pharma_ids)
})

##
OUTDIR <- file.path(PDIR, "output", "univariate", analysis)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

#############################
#### UNIVARIATE ANALYSIS ####
lapply(all_ds, function(dataset){
  print(dataset)
  if(!file.exists(file.path(OUTDIR, paste0(dataset, "_", analysis, ".rds")))){
    dataset_drugs <- lapply(drugs_oi, function(drug_x){
      cat(paste0("Drug: ", drug_x, "...\n"))
      # gm_ds[[dataset]]
      # ds_pharma_ids[[dataset]]
      
      pharma_x <- ds_pharma_ids[[dataset]]
      na_idx <- which(is.na(pharma_x))
      na_idx <- unique(c(na_idx, which(is.na(match(pharma_x, colnames(psets[[dataset]]))))))
      pharma_x <- pharma_x[-na_idx]
      gm_x <- gm_ds[[dataset]][-na_idx, -ncol(gm_ds[[dataset]])]
      
      aac_x <- psets[[dataset]][drug_x, match(pharma_x, colnames(psets[[dataset]]))]
      gm_xord <- gm_x[match(colnames(aac_x), pharma_x),]
      na_sum <- colSums(is.na(gm_xord))
      gm_xord <- gm_xord[,-which(na_sum > 100)]
      
      gain_idx <- gm_xord > 4
      loss_idx <- gm_xord < 2
      neut_idx <- gm_xord >= 2 & gm_xord <= 4
      gm_class <- gm_xord
      gm_class[loss_idx] <- 'LOSS'
      gm_class[neut_idx] <- 'NEUT'
      gm_class[gain_idx] <- 'GAIN'
      
      idx <- 1
      p_class <- apply(gm_class, 2, function(i){
        # cat(paste0("\t> ", colnames(gm_class)[idx], "\n"))
        idx <<- idx + 1
        
        #cor(as.numeric(i), as.numeric(aac_x), method = 'spearman', use = 'complete.obs')
        datax=data.frame("aac"=as.numeric(aac_x), "cn"=i)
        aov_result <- aov(aac ~ cn, datax)
        wilcox_result <- tryCatch({
          pairwise.wilcox.test(datax$aac, datax$cn)
        }, error=function(e){list("p.value"=rep(NA, 4))})
        
        if(is.matrix(wilcox_result$p.value)){
          wilcox_ids <- as.vector(sapply(colnames(wilcox_result$p.value), function(x){
            paste(gsub("^(.).*", "\\1",x), 
                  gsub("^(.).*", "\\1",rownames(wilcox_result$p.value)), sep="-")
          }))
        } else {
          wilcox_ids <- c('G-L', 'G-N', 'L-L', 'L-N')
        }
        
        c("AOV"=summary(aov_result)[[1]]$Pr[1], 
          setNames(round(as.vector(wilcox_result$p.value), 5), wilcox_ids))
      })
      if(is.list(p_class)){
        p_class <- plyr::rbind.fill(lapply(p_class, function(x) as.data.frame(t(x))))
        rownames(p_class) <- colnames(gm_class)[as.integer(rownames(p_class))]
      } else {
        p_class <- as.data.frame(t(p_class))
      }
      p_class <- p_class[order(p_class$AOV),]
      return(p_class)
    })
    names(dataset_drugs) <- drugs_oi
    saveRDS(dataset_drugs, file=file.path(OUTDIR, paste0(dataset, "_", analysis, ".rds")))
    return(NULL)
  }
})


##############################
#### SET UP REFERENCE POS ####
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(AnnotationDbi)
library(BSgenome.Hsapiens.UCSC.hg19)

## Gene Ref
genes0 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene )
genes0$symbol <- mapIds(org.Hs.eg.db,
                        keys=genes0$gene_id,
                        column="SYMBOL",
                        keytype="ENTREZID",
                        multiVals="first")
chrs <- paste0("chr", c(1:22, "X", "Y"))
chr_sizes <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[chrs]
cum_chr_sizes <- setNames(c(1, cumsum(as.numeric(chr_sizes))+1), c(chrs, 'chrZ'))
genes0_chr <- split(genes0, seqnames(genes0))[chrs]
genes0_chr <- lapply(names(genes0_chr), function(chr_id){
  cumstart <- start(genes0_chr[[chr_id]]) + (cum_chr_sizes[chr_id] - 1)
  cumend <- end(genes0_chr[[chr_id]]) + (cum_chr_sizes[chr_id] - 1)
  genes0_chr[[chr_id]]$cumstart <- cumstart
  genes0_chr[[chr_id]]$cumend <- cumend
  genes0_chr[[chr_id]]$cummid <- cumstart + ((cumend - cumstart)/2)
  genes0_chr[[chr_id]]
})
genes0 <- sort(unlist(as(genes0_chr, "GRangesList")))

## Bins Ref
bins_ref <- readRDS(file.path(PDIR, "input", 'bin', "bins_ref.rds"))
bins_ref$chrs <- sapply(strsplit(bins_ref$ID, split="_"), function(i){i[1]})
bins_ref$idx <- sapply(strsplit(bins_ref$ID, split="_"), function(i){i[2]})
bins_ref$cumend <- cumsum(as.numeric(width(bins_ref)))
bins_ref$cumstart <- c(1, (bins_ref$cumend[-length(bins_ref)]+1))
bins_ref$chg <- FALSE
bins_ref <- sapply(split(bins_ref, bins_ref$chrs), function(i){
  i[c(1, length(i)),]$chg <- TRUE
  return(i)
})
bins_ref <- sort(unlist(as(bins_ref, "GRangesList")))

############################
#### COMBINE UNIVARIATE ####
thresholds <- c(0.01, 0.05, 0.1, 0.15)
thresholds <- setNames(thresholds, thresholds)

drug_gene_list <- aggregateSigUni(analysis='gene', padj_method='bonferroni')
drug_bin_list <- aggregateSigUni(analysis='bin', padj_method='bonferroni')

tgenes_pos <- matchSigToRef(drug_gene_list, genes0, 'symbol')
tbins_pos <- matchSigToRef(drug_bin_list, bins_ref, 'ID')

# as.data.frame(tgenes_pos$Lapatinib[order(tgenes_pos$Lapatinib$thresh),]) # HER2 inhibitor (ERBB2)
# as.data.frame(tgenes_pos$Erlotinib[order(tgenes_pos$Erlotinib$thresh),]) # EGFR inhibitor (EGFR, EGFR-AS1)
# as.data.frame(tgenes_pos$Pictilisib[order(tgenes_pos$Pictilisib$thresh),]) # PI3K inhibitor (PTEN, PIK3CA, BRAF)
# as.data.frame(tgenes_pos$Sirolimus[order(tgenes_pos$Sirolimus$thresh),]) # rapamycin analog (PTEN, PIK3CA, BRAF) 

#####################
#### MAIN FIGURE ####
qval_breaks <- c(0, 0.001, 0.01, 0.05)
breaks_id <- levels(cut(seq(0.00001, 0.05, by=0.0001), qval_breaks))
thresh_cols <- list("genes"=setNames(rev(c("#fec44f", "#fe9929", "#d95f0e")),
                                     as.character(breaks_id)),
                    "bins"=setNames(rev(c("#e5f5e0", "#a1d99b", "#41ab5d")),
                                    as.character(breaks_id)))

FIG_OUTPUT <- file.path(PDIR, "output", "univariate", "figures")
dir.create(FIG_OUTPUT, recursive = TRUE, showWarnings = FALSE)
# png(file.path(FIG_OUTPUT, "univariate_genome.png"), units = 'in', 
#     width = 12, height = 5, res=200)
pdf(file.path(FIG_OUTPUT, "univariate_genome.pdf"), width = 12, height = 5)
par(mar=c(5.1, 7, 4.1, 2.1))
## Blank Drug x Chromosome plot
plot(0, type='n', xlim=c(1, max(bins_ref$cumend)), 
     ylim=c(0.5, (length(drugs_oi)+1)), yaxt='n', xaxt='n',
     ylab='', xlab='Chromosomes')
abline(v = cum_chr_sizes, lty=2, col="grey")

## Row/Col labels
axis(side=1, at = (cum_chr_sizes + c(diff(cum_chr_sizes)/2, 1e6))[c(TRUE, FALSE)], 
     labels=gsub("^chr", "", names(cum_chr_sizes))[c(TRUE, FALSE)], las=1, 
     tick = FALSE, cex.axis=0.8)
axis(side=1, at = (cum_chr_sizes + c(diff(cum_chr_sizes)/2, 1e6))[c(FALSE, TRUE)], 
     labels=gsub("^chr", "", names(cum_chr_sizes))[c(FALSE, TRUE)], las=1, 
     tick = FALSE, cex.axis=0.8, pos = -0.1)
axis(side=1, at=cum_chr_sizes, labels=rep("", length(cum_chr_sizes)))

axis(side = 2, at = seq_along(drugs_oi), labels = drugs_oi, las=1)

legend("topleft", pch=rep(16,4), col = thresh_cols[['genes']],
       legend = names(thresh_cols[['genes']]), horiz=TRUE, box.lwd = 0)
legend("topright", pch=rep(15,4), col = thresh_cols[['bins']],
       legend = names(thresh_cols[['bins']]), horiz=TRUE, box.lwd = 0)

## Fill in the individual sig genes (points) and bins
for(drug_idx in seq_along(drugs_oi)){
  print(paste0(drug_idx, " - ", drugs_oi[drug_idx]))
  sgp <- tgenes_pos[[drug_idx]]
  sbp <- tbins_pos[[drug_idx]]
  
  ## Map bins ID to get the cumulative position from bins_ref
  if(length(sbp)> 0){
    rect(xleft = sbp$cumstart, ybottom = drug_idx-0.3, 
         xright = sbp$cumend, ytop = drug_idx+0.3, 
         col=thresh_cols[['bins']][as.character(cut(sbp$q, qval_breaks))], 
         border = thresh_cols[['bins']][as.character(cut(sbp$q, qval_breaks))])
  }
  
  ## Find cumulative position for gene position
  if(length(sgp) > 0){
    points(y=rep(drug_idx, length(sgp)), x=sgp$cummid, cex=0.6,
           pch=16, col=thresh_cols[['genes']][as.character(cut(sgp$q, qval_breaks))])
  }
}
dev.off()

#####################
#### SUPP FIGURE ####
drugs_gene_aov_mu <- getMeanAov(analysis='gene')
drugs_bin_aov_mu <- getMeanAov(analysis='bin')

mu_genes_pos <- lapply(drugs_gene_aov_mu, function(gene) {
  geneX <- genes0[match(names(gene), genes0$symbol),]
  geneX$AOV <- gene
  return(geneX)
})
mu_bins_pos <- lapply(drugs_bin_aov_mu, function(bin) {
  binX <- bins_ref[match(names(bin), bins_ref$ID),]
  binX$AOV <- bin
  return(binX)
})


bin_pos <- thresh_bins_pos[['0.05']]
gene_pos <- thresh_genes_pos[['0.05']]
gene_pos <- tgenes_pos 
bin_pos <- tbins_pos 

for(borg in c('b_only', 'g_only')){
  gib_tbl <- sapply(setNames(drugs_oi, drugs_oi), function(drug){
    print(drug)
    bin <- bin_pos[[drug]]
    gene <- gene_pos[[drug]]
    
    # Identify the gene features [bin or gene] that overlaps/intersects with the bin
    ov_idx <- findOverlaps(gene, bin)
    bg_ov <- list("intersect"=gene[unique(queryHits(ov_idx)),])

    # Identify non-overlapping features for the gene
    if(length(gene) != length(unique(queryHits(ov_idx)))){
      # If there are gene features not overlapping the bin
      borg <- 'g_only'
      if(length(ov_idx) > 0){
        bg_ov[[borg]] <- gene[c(1:length(gene))[-queryHits(ov_idx)],]
      } else {
        bg_ov[[borg]] <- gene
      }
    } else {
      # If all gene features overlap the bin
      borg <- 'g_only'
      bg_ov[[borg]] <- gene[0,]
    }
    
    # Identify non-overlapping features for the bin
    if(length(bin) != length(unique(queryHits(ov_idx)))){
      # If there are bin features not overlapping the bin
      borg <- 'b_only'
      if(length(ov_idx) > 0){
        bg_ov[[borg]] <- bin[c(1:length(bin))[-subjectHits(ov_idx)],]
      } else {
        bg_ov[[borg]] <- bin
      }
    } else {
      # If all bin features overlap the bin
      borg <- 'b_only'
      bg_ov[[borg]] <- bin[0,]
    }
    
    ov_cnts <- sapply(bg_ov[c(2,1,3)], function(i){sum(width(i))})
    t(t(ov_cnts / sum(ov_cnts)))
  })
  rownames(gib_tbl) <- c('g_only', 'intersect', 'b_only')
  cols <- setNames(c("#d95f0e", "#grey", "#41ab5d"), rownames(gib_tbl))
  
  pdf(file.path(FIG_OUTPUT, "univariate-gib.pdf"), height = 5, width = 6)
  par(mar=c(7, 4.1, 4.1, 2.1))
  bp_x <- barplot(gib_tbl, col=cols,border = NA, las=2, xaxt='n')
  axis(side = 1, at = bp_x, labels = colnames(gib_tbl), las=2, tick = FALSE, line = 0)
  axis(side = 3, at = c(-0.2, bp_x), labels = c("n", sapply(bin_pos, length)), 
       las=1, tick = FALSE, line=-1)
  axis(side = 3, at = c(-0.2, bp_x), labels = c("n", sapply(gene_pos, length)), 
       las=1, tick = FALSE, line=0)
  dev.off()
  
  # summary metrics
  apply(gib_tbl, 1, function(x){c("mean"=mean(x), "sd"=sd(x))})
  #      g_only    intersect   b_only 
  # mean 0.1011525 0.2875470   0.6113005
  # sd   0.1183739 0.2056761   0.2720461
}

