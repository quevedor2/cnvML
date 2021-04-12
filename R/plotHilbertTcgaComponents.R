# Purpose: ?

# Description:  Uses the importance feature map created from the 
## GradCAM approach and maps it back to the individual samples
## to extract what the meaningful features are from each

library(HilbertCurve)
library(HilbertVis)
library(circlize)
library(DescTools)
library(GenomicRanges)
library(dplyr)
library(RcnvML)   # devtools::install_github("quevedor2/cnvML/RcnvML")


###################
#### Functions ####
overlapMatToRefHc <- function(alt_mat, ref_mat){
  nr2 <- nrow(alt_mat)
  nr1 <- nrow(ref_mat)
  lcm_nr <- pracma::Lcm(nr1, nr2)
  nr2_mult <- lcm_nr / nr2
  nr1_mult <- lcm_nr / nr1
  
  alt_matE <- t(apply(alt_mat, 1, rep, each=nr2_mult))
  alt_matE <- t(apply(alt_matE, 2, rep, each=nr2_mult))
  
  if(!exists("gpos_matE")){
    print("Generating gpos_matE...")
    gpos_matE <<- t(apply(ref_mat, 2, rep, each=nr1_mult))
    gpos_matE <<- t(apply(gpos_matE, 2, rep, each=nr1_mult))
    gpos_matE <<- data.frame("gpos"=as.vector(gpos_matE))
  }
  reduced_cam <- sapply(split(as.vector(alt_matE), gpos_matE$gpos), median)
  return(reduced_cam)
}

mapHcToGpos <- function(cam, ref, returngr=FALSE){
  ord_preds <- cam[match(ref$uid, names(cam))]
  ref$val <- ord_preds
  
  if(returngr){
    require(GenomicRanges)
    # Extract Chr, Start, End pos
    ref$st_chrom <- gsub(":.*", "", ref$start)
    ref$en_chrom <- gsub(":.*", "", ref$end)
    ref$start <- gsub("^.*:", "", ref$start)
    ref$end <- gsub("^.*:", "", ref$end)
    
    
    # Deal with cases where bins extend multiple Chrs
    chrs <- getChrLength()
    span_idx <- which(ref$st_chrom != ref$en_chrom)
    if(length(span_idx) > 1){
      multi_span <- ref[span_idx,]
      chr_spl <- split(multi_span, f=multi_span$st_chrom)[c(1:22, "X")]
      multi_span <- do.call(rbind, lapply(names(chr_spl), function(chr){
        end_pos <- end(chrs[seqnames(chrs) == gsub("^(chr)?", "chr", chr)])
        new_chr <- rbind(chr_spl[[as.character(chr)]], chr_spl[[as.character(chr)]])
        if(!is.null(new_chr)){
          new_chr[1,c('end', 'en_chrom')] <- c(end_pos, new_chr$st_chrom[1])
          new_chr[2,c('start', 'st_chrom')] <- c(1, new_chr$en_chrom[2])
        }
        return(new_chr)
      }))
      
      ref <- rbind(ref, multi_span)
    }
    
    refgr <- makeGRangesFromDataFrame(ref[-span_idx, -grep("en_chrom", colnames(ref))], 
                                      seqnames.field = "st_chrom", keep.extra.columns = TRUE)
    ref <- refgr
  }
  return(ref)
}

##################
#### VARIABLE ####
dataset <- 'TCGA'
PDIR <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/TCGA'

seg_files <- setNames(c('TCGA_mastercalls.abs_segtabs.fixed.txt'), c("TCGA"))
ctypes <- list.files(file.path(PDIR, "ctype"), pattern = "[^tmp|rds]")


##############
#### MAIN ####
segf <- file.path(PDIR, "input", seg_files[dataset])
segd <- read.table(segf, sep="\t", header=TRUE, stringsAsFactors = FALSE)
segd <- segd[-which(is.na(segd$Chromosome)),]
seg_gr <- makeGRangesFromDataFrame(segd, keep.extra.columns = TRUE,
                                   start.field = 'Start', end.field = 'End')
seqlevelsStyle(seg_gr) <- 'ENSEMBL'
seg_grl <- as(split(seg_gr, f=seg_gr$Sample), "GRangesList")


gbin_pos <- setupRefHcMatrix(order=8)
drug_preds_imp <- lapply(setNames(ctypes,ctypes), function(ctype){
  ## Read in predictions
  pred_file <- list.files(file.path(PDIR, "ctype", ctype), pattern='pred.csv$')
  preds <- read.csv(file.path(PDIR, "ctype", ctype, pred_file), 
                    header = FALSE, stringsAsFactors = FALSE)
  
  ## Read in all CAM files
  cam_files <- list.files(file.path(PDIR, "ctype", ctype), pattern='[^pred].csv$')
  pred_class <- unique(gsub("^.*?_", "", cam_files) %>% gsub("_.*", "", .))
  use_idx <- which(preds[,2] == pred_class)
  
  if(!file.exists(file.path(PDIR, "ctype", ctype, paste0("top_", h_cutoff, ".RDS")))){
    sample_ids <- gsub("^.*_", "", cam_files) %>% gsub(".csv", "", .)
    top_preds <- lapply(cam_files[use_idx], function(f){
      x <- as.matrix(read.csv(file.path(PDIR, "ctype", ctype, f), header = FALSE))
    })
    names(top_preds) <- sample_ids[use_idx]
    
    
    ## Threshold CAM binary files into important/not-important features
    cutoff <- 0.95; h_cutoff <- 100
    binary_preds <- lapply(top_preds, function(i) {
      binary_mat <- matrix(data=0, ncol=ncol(i), nrow=nrow(i))
      # binary_mat[i > quantile(as.vector(i), cutoff)] <- 254
      binary_mat[i > h_cutoff] <- 254
      return(binary_mat)
    })
    
    ## Create ordering via h-clust
    mat_dist <- dist(t(sapply(binary_preds, as.vector)))
    mat_dist <- dist(t(sapply(top_preds, as.vector)))
    hcl <- hclust(mat_dist)
    
    pred_idx <- 1
    top_preds_gr <- lapply(binary_preds, function(alt_mat){
      print(pred_idx)
      preds_linear <- overlapMatToRefHc(alt_mat = alt_mat, ref_mat = gbin_pos$mat)
      preds_gpos <- mapHcToGpos(cam=preds_linear, ref=gbin_pos$ord, returngr=TRUE)
      
      pg <- preds_gpos[preds_gpos$val > 0,]
      int_seg_gr <- intersectAndAnnotate(s_gr = seg_grl[[names(binary_preds)[pred_idx]]], 
                                         ref_gr=pg)
      pred_idx <<- pred_idx + 1
      return(int_seg_gr)
    })
    top_preds_gr <- setNames(top_preds_gr[hcl$order], names(binary_preds)[hcl$order])
    saveRDS(top_preds_gr, file = file.path(PDIR, "ctype", ctype, paste0("top_", h_cutoff, ".RDS")))
  } else {
    top_preds_gr <- readRDS(file.path(PDIR, "ctype", ctype, paste0("top_", h_cutoff, ".RDS")))
  }
  
  chrs <- getChrLength()
  pdf(file.path(PDIR, "output", "ctype", paste0(ctype, ".pdf")), 
      height=6, width=3)
  # top_preds_gr <<- top_preds_gr
  plot(x=0, type='n', xlim=c(0, (length(top_preds_gr)+1)), 
       ylim=c(0, max(chrs$cum.end)), 
       yaxt='n', xaxt='n', ylab='', xlab='')
  axis(side = 2, at = chrs$cum.mid, las=2,
       labels = gsub("^chr", "", as.character(seqnames(chrs))), cex.axis=0.7)
  
  cn_cols <- setNames(c("black", "thistle", "red", "darkred"), c(0:3))
  cn_cols <- setNames(c("blue", "deepskyblue", "grey", "salmon", "red"), 
                        c('deeploss', 'loss', 'neutral', 'gain', 'amp'))
  loh_cols <- setNames(c('white', 'black'), c('het', 'loh'))
  
  idx <<- 1
  lapply(top_preds_gr, function(gr0){
    if(length(gr0) > 0){
      gr0$classification <- 'neutral'
      gr0$classification[gr0$Modal_Total_CN < 2] <- 'loss'
      gr0$classification[gr0$Modal_Total_CN == 0] <- 'deeploss'
      gr0$classification[gr0$Modal_Total_CN > 2] <- 'gain'
      gr0$classification[gr0$Modal_Total_CN > 4] <- 'amp'
      gr0$loh <- 'het'
      gr0$loh[gr0$Modal_HSCN_1 == 0 | gr0$Modal_HSCN_2 == 0] <- 'loh'
      
      rect(xleft = idx-0.45, ybottom = gr0$cum.start - 250000, 
           xright = idx+0.45, ytop = gr0$cum.end + 250000, 
           col = loh_cols[gr0$loh], border = NA)
      rect(xleft = idx-0.3, ybottom = gr0$cum.start - 250000, 
           xright = idx+0.3, ytop = gr0$cum.end + 250000, 
           col = cn_cols[gr0$classification], border = NA)
    } else {
      idx <<- idx - 1
    }
    idx <<- idx + 1
  })
  abline(h = c(1, chrs$cum.end), lty=2, col="grey", cex=0.5)
  dev.off()
  
  # 
  # 
  # cols <- RColorBrewer::brewer.pal(if(length(top_preds) > 6) 6 else length(top_preds), "Dark2")
  # print("feat")
  # pdf(file.path(gsub("input", "output", PDIR), "drugs", dataset, "feat_mats", 
  #               paste0("feat_", drug, "_", actual, ".pdf")))
  # plot.new()
  # #for(i in c(1:length(cols))){
  # #rasterImage(alphaIt(top_preds[[i]], col=cols[i], 1), 0, 0, 1, 1)
  # rasterImage(alphaIt(abs(all_preds2), 'red', 1),   0, 0, 1, 1)
  # for(i in c(1:length(top_preds))){
  #   rasterImage(alphaIt(abs(top_preds[[i]] - mid_avg), col='black', scale/length(top_preds)), 0, 0, 1, 1)
  # }
  # dev.off()
})
