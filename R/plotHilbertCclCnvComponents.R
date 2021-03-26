library(HilbertCurve)
library(HilbertVis)
library(circlize)
library(DescTools)
library(GenomicRanges)

###################
#### Functions ####
## ## ##
## Visualization functions
alphaIt <- function(matX, col, alpha=0.5){
  #x <- abs(1-(matX/255))
  x <- (matX/255)
  
  col_ramp <- colorRampPalette(colors = c("white", col))(255)
  feat_cols <- col_ramp[matX+1]
  feat_rgb <- col2rgb(feat_cols)
  
  x2 <- abind::abind(matrix(feat_rgb[1,]/255, ncol=ncol(matX)), 
                     matrix(feat_rgb[2,]/255, ncol=ncol(matX)), 
                     matrix(feat_rgb[3,]/255, ncol=ncol(matX)), 
                     x/(1/alpha), along=3 )
  return(x2)
}

## ## ##
## Map Chr Pos to Hilbert
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

calcCumulativeGPos <- function(gr0, seqstyle='ENSEMBL'){
  chrs <- getChrLength()
  seqlevelsStyle(chrs) <- seqstyle
  
  gr0 <- unlist(as(lapply(split(gr0, seqnames(gr0)), function(chr_gr){
    # Matches chr position
    chr_id <- unique(as.character(seqnames(chr_gr)))
    ref_idx <- match(chr_id, as.character(seqnames(chrs)))
    # Calculates cumulative position
    chr_gr$cum.start <- start(chr_gr) + chrs[ref_idx,]$cum.start - 1
    chr_gr$cum.end<- end(chr_gr) + chrs[ref_idx,]$cum.start - 1
    return(chr_gr)
  }), "GRangesList"))
  return(gr0)
}

intersectAndAnnotate <- function(s_gr, ref_gr, cpos=TRUE){
  ## Intersects s_gr with ref_gr, then annotates all of s_gr component
  ## on the new intersect
  print(unique(s_gr$ID))
  end(s_gr) <- end(s_gr)-1 # Prevents merging of adjacent segments
  int_gr <- intersect(ref_gr, s_gr, ignore.strand=FALSE)
  ov_idx <- findOverlaps(int_gr, s_gr)
  mcols(int_gr) <- mcols(s_gr[subjectHits(ov_idx),])
  if(cpos) int_gr <- calcCumulativeGPos(int_gr)
  
  return(int_gr)
}

##################
#### VARIABLE ####
dataset <- 'CCLE'

segDIR <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines'
seg_files <- setNames(c('CCLE_cna_hg19.seg', 'GDSC_cna_hg19.seg', 'gCSI_cna_hg19.seg'), 
                      c("CCLE", "GDSC", "GNE"))
PDIR <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnvML/CCL/input'
meta_file <- file.path(PDIR, 'meta_df.csv')
drugs <- list.files(file.path(PDIR, "drugs", dataset), pattern = "[^tmp|rds]")
drug <- 'Docetaxel'


##############
#### MAIN ####
meta_df <- read.table(meta_file, sep=",", stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
segf <- file.path(segDIR, "ccl_aggregate", "input", seg_files[dataset])
segd <- read.table(segf, sep="\t", header=TRUE, stringsAsFactors = FALSE)
seg_gr <- makeGRangesFromDataFrame(segd, keep.extra.columns = TRUE)
seqlevelsStyle(seg_gr) <- 'ENSEMBL'
seg_grl <- as(split(seg_gr, f=seg_gr$ID), "GRangesList")


gbin_pos <- setupRefHcMatrix(order=8)

mid <- '5'
drug_preds_imp <- lapply(setNames(drugs,drugs), function(drug){
  out_rds <- file.path(gsub("input", "output", PDIR), "drugs", dataset, 
                       "feat_mats", paste0(drug, ".RDS"))
  
  dir <- file.path(PDIR, "drugs", dataset, drug, mid)
  campred <- read.table(file.path(dir, "CAM-pred.csv"), sep=",")
  pred_eq_idx <- which(campred[,2] == mid)
  mid_preds <- lapply(pred_eq_idx, function(idx){
    as.matrix(read.csv(file.path(dir, paste0("CAM-", drug, "_naive_", idx-1, ".csv")), header = FALSE))
  })
  mid_avg <- round(Reduce(function(x,y) x + y, mid_preds) / length(mid_preds),0)
  
  
  if(!file.exists(out_rds)){
    preds_imp <- lapply(setNames(c("0", "9"), c("0", "9")), function(actual){
      print(paste0(drug, " - ", actual))
      print(paste0(drug, "-", actual))
      dir <- file.path(PDIR, "drugs", dataset, drug, actual)
      
      campred <- read.table(file.path(dir, "CAM-pred.csv"), sep=",")
      pred_eq_idx <- which(campred[,2] == actual)
      samples <- as.character(campred[pred_eq_idx,1])
      ds_samples <- meta_df[match(samples, meta_df$PharmacoGX_ID), dataset]
      
      top_preds <- lapply(pred_eq_idx, function(idx){
        as.matrix(read.csv(file.path(dir, paste0("CAM-", drug, "_naive_", idx-1, ".csv")), header = FALSE))
      })
      
      scale <- 3
      top_preds2 <- lapply(top_preds, function(x) abs(x-mid_avg))
      top_preds2 <- lapply(top_preds2, function(i) {
        binary_mat <- matrix(data=0, ncol=ncol(i), nrow=nrow(i))
        binary_mat[i > quantile(as.vector(i), 0.99)] <- 254
        return(binary_mat)
      })
      pred_idx <- 1
      mat_dist <- dist(t(sapply(top_preds2, as.vector)))
      hcl <- hclust(mat_dist)
      
      top_preds_gr <- lapply(top_preds2, function(alt_mat){
        print(pred_idx)
        preds_linear <- overlapMatToRefHc(alt_mat = alt_mat, ref_mat = gbin_pos$mat)
        preds_gpos <- mapHcToGpos(cam=preds_linear, ref=gbin_pos$ord, returngr=TRUE)
        
        pg <- preds_gpos[preds_gpos$val > 0,]
        int_seg_gr <- intersectAndAnnotate(s_gr = seg_grl[ds_samples][[pred_idx]], 
                                           ref_gr=pg)
        pred_idx <<- pred_idx + 1
        return(int_seg_gr)
      })
      top_preds_gr <- setNames(top_preds_gr[hcl$order], ds_samples[hcl$order])
      
      
      
      
      all_preds2 <- all_avg <- round(Reduce(function(x,y) abs((x + y)-mid_avg), top_preds) / length(top_preds))
      all_preds2[,] <- 0
      all_preds2[all_avg > quantile(as.vector(all_avg), 0.99)] <- 254
      
      all_preds_linear <- overlapMatToRefHc(alt_mat = all_preds2, ref_mat = gbin_pos$mat)
      preds_gpos <- mapHcToGpos(cam=all_preds_linear, ref=gbin_pos$ord, returngr=TRUE)
      
      pg <- preds_gpos[preds_gpos$val > 0,]
      pg_red <- unlist(reduce(split(pg, ~val)))
      pg_red$val <- names(pg_red)
      names(pg_red) <- c(1:length(pg_red))
      # pg_red[seqnames(pg_red) == '7',]
      
      int_seg_gr <- lapply(seg_grl[ds_samples], intersectAndAnnotate, 
                           ref_gr=pg)
      chrs <- getChrLength()
      
      print("predbar")
      pdf(file.path(gsub("input", "output", PDIR), "drugs", dataset, "feat_mats", 
                    paste0("predbar_", drug, "_", actual, ".pdf")), 
          height=7, width=4)
      # top_preds_gr <<- top_preds_gr
      plot(x=0, type='n', xlim=c(0, (length(top_preds_gr)+1)), 
           ylim=c(0, max(chrs$cum.end)), 
           yaxt='n', xaxt='n', ylab='', xlab='')
      axis(side = 2, at = chrs$cum.mid, las=2,
           labels = gsub("^chr", "", as.character(seqnames(chrs))))
      idx <<- 1
      cn_cols <- setNames(c("black", "thistle", "red", "darkred"), c(0:3))
      
      #apply(int_seg_gr, function(gr0){
      lapply(top_preds_gr, function(gr0){
        gr0$nMajor[gr0$nMajor > 3] <- 3
        rect(xleft = idx-0.3, ybottom = gr0$cum.start, 
             xright = idx, ytop = gr0$cum.end, 
             col = cn_cols[as.character(gr0$nMajor)], border = NA)
        
        rect(xleft = idx, ybottom = gr0$cum.start, 
             xright = idx+0.3, ytop = gr0$cum.end, 
             col = cn_cols[as.character(gr0$nMinor)], border = NA)
        idx <<- idx + 1
      })
      abline(h = chrs$cum.end, lty=2, col="grey")
      dev.off()
      
      
      cols <- RColorBrewer::brewer.pal(if(length(top_preds) > 6) 6 else length(top_preds), "Dark2")
      print("feat")
      pdf(file.path(gsub("input", "output", PDIR), "drugs", dataset, "feat_mats", 
                    paste0("feat_", drug, "_", actual, ".pdf")))
      plot.new()
      #for(i in c(1:length(cols))){
      #rasterImage(alphaIt(top_preds[[i]], col=cols[i], 1), 0, 0, 1, 1)
      rasterImage(alphaIt(abs(all_preds2), 'red', 1),   0, 0, 1, 1)
      for(i in c(1:length(top_preds))){
        rasterImage(alphaIt(abs(top_preds[[i]] - mid_avg), col='black', scale/length(top_preds)), 0, 0, 1, 1)
      }
      dev.off()
      
      return(list("individual"=top_preds_gr,
                  "aggregate"=int_seg_gr))
    })
    saveRDS(preds_imp, file=out_rds)
  } else {
    preds_imp <- readRDS(out_rds)
  }
  
  return(preds_imp)
})




