library(dplyr)
library(gtools)
library(mclust)
library(dplyr)

###################
#### FUNCTIONS ####
plotClassification <- function(aac_dat, drug_id){
  col_breaks <- c(0:40)
  col_palette <- colorRampPalette(c("white", "#3182bd"))(length(col_breaks)-1)
  
  ## Plots confusion matrix using a 0-40 scale
  gplots::heatmap.2(t(table(aac_dat$Pred, aac_dat$Actual)), 
                    Rowv = FALSE, Colv = FALSE, trace = 'none', 
                    breaks=col_breaks, col = col_palette,
                    dendrogram = 'none', main=drug_id)
}

plotMixedAAC <- function(aac_dat, drug_id, mfit, cuts, quant_cutoff, Q){
  ## Plots AAC density and mixed gaussian representation + uncertainty
  # Initialize Blank Plot
  plot(0, type='n', xlim=c(0,1), ylim=c(0,2.25), main=drug_id, 
       las=1, ylab="Density", xlab="AAC", yaxt='n')
  
  # AAC Density
  dens_AAC <- density(aac_dat$AAC)
  dens_AAC$y <- dens_AAC$y/max(dens_AAC$y)
  lines(dens_AAC)
  segments(x0 = (dens_AAC$x), y0 = rep(0, length(dens_AAC$y)), 
           x1 = (dens_AAC$x), y1 = dens_AAC$y, col="#9ecae1")
  axis(side=2, at = c(0, 0.5, 1, 1.25, 1.75, 2.25), 
       labels = c(0, 0.5, 1, 0, 0.5, 1), las=1)
  
  # Add classification uncertainty track
  segments(x0 = mfit$data, y0 = rep(0, length(mfit$uncertainty)) + 1.25, 
           x1 = mfit$data, y1 = mfit$uncertainty + 1.25, col='grey')
  
  # Add mixed gaussian models (G=2)
  dX1 <- density(rnorm(n = 10000, mean = mfit$parameters$mean[1], sd = sqrt(mfit$parameters$variance$sigmasq)))
  dX1$y <- dX1$y / max(dX1$y)
  dX2 <- density(rnorm(n = 10000, mean = mfit$parameters$mean[2], sd = sqrt(mfit$parameters$variance$sigmasq)))
  dX2$y <- dX2$y / max(dX2$y)
  lines(dX1$x, dX1$y + 1.25, col="red")
  lines(dX2$x, dX2$y + 1.25, col="blue")
  
  abline(v = c(0, cuts), lty=2, col="darkgrey")
  rect(xleft = cuts[quant_cutoff], ybottom = -1, xright = 10, ytop = 5,
       col=scales::alpha("green", 0.1), border = NA)
  text(x = cuts - (diff(c(0, cuts))/2), y= rep(1.12, Q), c(0:(Q-1)))
}

plotF1Comp <- function(wF1, three_cols){
  par(mar=c(6.1, 4.1, 4.1, 2.1))
  
  max_xlim <- nrow(wF1)*(ncol(wF1)+1)
  bp_x <- barplot(t(wF1), beside=T, xaxt='n', las=1, 
                  xlim=c(0, (max_xlim + 5)), ylim=c(0,1),
                  col = three_cols[colnames(wF1)],
                  ylab="F1 score")
  axis(side = 1, at = bp_x[2,], las=2,
       labels = gsub("_AAC.*", "", rownames(wF1)))
  # 
  # summ_xidx <- sapply(colnames(wF1), function(colid){
  #   i <- wF1[,colid]
  #   xidx <- ((grep(colid, colnames(wF1)) - 1) * 3) + 1.5 + max_xlim
  #   xidx <- grep(colid, colnames(wF1)) + 0.5 + max_xlim
  #   w <- 0.4
  #   
  #   segments(x0 = xidx, y0 = min(i), x1 = xidx, y1 = max(i))
  #   rect(xleft = xidx - w, ybottom = quantile(i, 0.25), 
  #        xright = xidx + w, ytop = quantile(i, 0.75), col=three_cols[colid])
  #   segments(x0 = xidx-w, y0 = median(i), x1 = xidx+w, y1 = median(i))
  #   return(xidx)
  # })
  # axis(side = 1, at = summ_xidx, las=2, labels = colnames(wF1))
  
  legend("topleft", fill = three_cols[colnames(wF1)], 
         legend = colnames(wF1))
}

getPredF1Ratio <- function(df, three_cols){
  par(mfrow=c(1,3), mar=c(5.1, 2, 2.1, 1))
  summ_dat <- apply(df, 2, function(x){
    xy_df <- apply(df, 2, function(y){
      x/y
    })
    idx <- which(sapply(apply(xy_df, 2, unique), length)==1)
    boxplot(xy_df, ylim=c(0.5, 1.5), col=three_cols[colnames(xy_df)], las=1)
    abline(h = 1)
    
    data.frame("mean"=colMeans(xy_df), 
               "sd"=apply(xy_df, 2, sd),
               "p"=apply(xy_df, 2, function(i) {
                 tryCatch({t.test(xy_df[,idx], i, alternative = 'less')$p.val}, error=function(e){1})
               }))
  })
  return(summ_dat)
}

###################
#### MAIN ####
three_cols <- c('hilbert'='#a6cee3',
                'gene'='#1f78b4',
                'bin'='#b2df8a')
PDIR <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnvML/CCL'
Q <- 10 # number of quantiles to split into

type_ds_metrics <- lapply(c('hilbert', 'gene', 'bin'), function(data_type){
  predDIR=file.path(PDIR, "input", data_type)
  ds_metrics <- lapply(c("GDSC", "GNE"), function(ds){
    metrics <- sapply(list.files(file.path(predDIR, ds, "drugs"), pattern="csv$"), function(f){
      # data_type='hilbert'
      # predDIR=file.path(PDIR, "input", data_type)
      # ds <- 'GDSC'
      # f = list.files(file.path(predDIR, ds, "drugs"), pattern="csv$")[1]
      
      aac_dat = read.csv(file.path(predDIR, ds, "drugs", f), 
                         header = TRUE, stringsAsFactors = FALSE)
      rownames(aac_dat) <- aac_dat[,1]
      aac_dat <- aac_dat[,-1]
      
      
      i_grp <- gtools::quantcut(as.numeric(aac_dat$AAC), q=Q)
      cuts <- sapply(levels(i_grp), function(j) as.numeric(gsub("\\]", "", strsplit(j, split=",")[[1]][2])))
      #hist(aac_dat$AAC, breaks=100, xlim=c(0,1))
      #abline(v=cuts, col="red")
      
      mfit <- Mclust(as.numeric(aac_dat$AAC), G=2)
      aac_dat$uncertainty <- mfit$uncertainty
      aac_dat$n_uncertainty <- ((1-aac_dat$uncertainty) * 2)
      aac_dat$mixture <- mfit$classification
      
      cutoff <- aac_dat$AAC[which.max(aac_dat$uncertainty)]
      quant_cutoff <- sum(cuts <= cutoff) 
      # -1 to adjust for python 0-based numbering of classes
      aac_dat$pred_mixture <- as.integer(aac_dat$Pred > (quant_cutoff-1))
      aac_dat$actual_mixture <- as.integer(aac_dat$Actual > (quant_cutoff-1))
      
      TP_idx <- which(aac_dat$pred_mixture == 1  & aac_dat$actual_mixture == 1)
      FP_idx <- which(aac_dat$pred_mixture == 1  & aac_dat$actual_mixture == 0)
      FN_idx <- which(aac_dat$pred_mixture == 0  & aac_dat$actual_mixture == 1)
      w_precision <- sum(aac_dat$n_uncertainty[TP_idx]) / sum(aac_dat$n_uncertainty[c(TP_idx, FP_idx)]) 
      w_recall <- sum(aac_dat$n_uncertainty[TP_idx]) / sum(aac_dat$n_uncertainty[c(TP_idx, FN_idx)])
      w_F1 <- 2 * ((w_precision * w_recall) / (w_precision + w_recall))
      
      precision <- length(aac_dat$n_uncertainty[TP_idx]) / length(aac_dat$n_uncertainty[c(TP_idx, FP_idx)]) 
      recall <- length(aac_dat$n_uncertainty[TP_idx]) / length(aac_dat$n_uncertainty[c(TP_idx, FN_idx)])
      F1 <- 2 * ((precision * recall) / (precision + recall))
      
      
      
      drug_id <- gsub("_.*", "", f)
      nort <- gsub("^.*_", "", f) %>% gsub(".csv", "", .)
      pdf(file.path(PDIR, "output", 
                    paste0(data_type, "-", nort, '-', drug_id, "-", ds, "_confusion.pdf")))
      plotClassification(aac_dat, drug_id)
      dev.off()
      pdf(file.path(PDIR, "output", paste0(drug_id, "-", ds, "_param.pdf")))
      plotMixedAAC(aac_dat, drug_id, mfit, cuts, quant_cutoff, Q)
      dev.off()
      
      return(c("wF1"=w_F1, "wP"=w_precision, "wR"=w_recall, 
               "F1"=F1, "P"=precision, "R"=recall, 
               "cutoff"=cutoff, "Q"=quant_cutoff, "Qcutoff"=cuts[quant_cutoff]))
    })
    metrics <- as.data.frame(t(metrics))
  })
  names(ds_metrics) <- c("GDSC", "GNE")
  return(ds_metrics)
})
names(type_ds_metrics) <- c('hilbert', 'gene', 'bin')

feat_ratios <- lapply(setNames(c("GDSC", "GNE"), c("GDSC", "GNE")), function(dataset){
  wF1 <- sapply(type_ds_metrics, function(i) i[[dataset]]$wF1)
  rownames(wF1) <- rownames(type_ds_metrics$hilbert[[dataset]])
  naive_idx <- grep("naive", rownames(wF1))
  transfer_idx <- grep("transfer", rownames(wF1))
  
  pdf(file.path(PDIR, "output", paste0("wF1-", dataset, "_naive-drugs.pdf")), height = 5)
  plotF1Comp(wF1=wF1[naive_idx,,drop=FALSE], three_cols)
  naive_feat_ratio <- getPredF1Ratio(wF1[naive_idx,,drop=FALSE], three_cols)
  dev.off()
  
  pdf(file.path(PDIR, "output", paste0("wF1-", dataset, "_transfer-drugs.pdf")), height = 5)
  plotF1Comp(wF1=wF1[transfer_idx,,drop=FALSE], three_cols)
  transfer_feat_ratio <- getPredF1Ratio(wF1[transfer_idx,,drop=FALSE], three_cols)
  dev.off()
  
  wF1_comp <- sapply(colnames(wF1), function(id){
    data.frame("mean_t"=mean(wF1[transfer_idx,id]),
               "sd_t"=sd(wF1[transfer_idx,id]),
               "mean_n"=mean(wF1[naive_idx,id]),
               "sd_n"=sd(wF1[naive_idx,id]),
               "p"=t.test(wF1[transfer_idx,id], wF1[naive_idx,id], alternative='less')$p.val)
  })
  
  return(list("wF1"=wF1,
              "comp"=wF1_comp,
              "transfer"=transfer_feat_ratio,
              "naive"=naive_feat_ratio))
})

## Compare naive to transfer
lapply(feat_ratios, function(i) i$comp)
#        hilbert    gene      bin      
# mean_t 0.384067   0.4532638 0.4737231
# sd_t   0.1749989  0.1815586 0.1924812
# mean_n 0.532091   0.5012335 0.5012891
# sd_n   0.1599017  0.1590826 0.1464654
# p      0.05698452 0.2992856 0.3816191


## Summarize p-values
lapply(feat_ratios, function(x) sapply(x, function(i) i$bin$p))
# $GDSC
#                 transfer       naive
# [Hilbert-gene]  0.9988853 0.005837593
# [Hilbert-bin]   0.9995822 0.034213853
# [Gene-Hilbert]  0.003519639 0.9948244
# [Bin-Hilbert]   0.001713259 0.9669958
# [Bin-Gene]      0.079983154 0.3598906


# $GNE
#                 transfer     naive
# [Hilbert-gene]  0.2340196 0.1166115
# [Hilbert-bin]   0.8833488 0.1354264
# [bin-gene]      0.10344188 0.1331073
# [bin-Hilbert]   0.07283785 0.8390692
# [bin-Gene]      0.10344188 0.1331073