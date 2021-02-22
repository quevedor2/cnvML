## Creates a circle-line plot for analysis done on NETseq
#### Functions ####
blankPlot <- function(xlim, ylim){
  plot(0, type='n', xlab='', ylab='',
       xaxt='n', yaxt='n', axes=FALSE, 
       xlim=xlim, ylim=ylim)
}

reduceArmMatrix <- function(ctype_dat){
  chr_ord <- as.character(sapply(c(1:22,"X"), function(chr) paste0(chr, c("p", "q"))))
  
  ## Assemble sample by chr-arm matrix
  chr_arm_mat <- Reduce(function(x,y) merge(x,y, by='Arm', all=TRUE), ctype_dat)
  rownames(chr_arm_mat) <- chr_arm_mat[,1]
  colnames(chr_arm_mat) <- c("Arm", gsub("-TP", "", names(ctype_df)))
  chr_arm_mat <- chr_arm_mat[,-1,drop=FALSE]
  
  chrs_coverage <- (!chr_ord %in% rownames(chr_arm_mat))
  if(any(chrs_coverage)){
    filler_df <- data.frame(matrix(NA, nrow = sum(chrs_coverage), ncol=ncol(chr_arm_mat)))
    rownames(filler_df) <- chr_ord[which(chrs_coverage)]
    colnames(filler_df) <- colnames(chr_arm_mat)
    chr_arm_mat <- rbind(chr_arm_mat, filler_df)
  }
  
  ## Binarize the matrix
  na_idx <- is.na(chr_arm_mat)
  chr_arm_mat[] <- 1
  chr_arm_mat[na_idx] <- 0
  
  ## Reorder matrix
  ord <- order(factor(rownames(chr_arm_mat), levels=chr_ord))
  chr_arm_mat <- chr_arm_mat[ord,]
  return(chr_arm_mat)
}

#### Parameters ####
line.width <- 0.055 # Size of the line going across platforms
point.size <- 2 # Size of the dots for identity matrix
text.size <- 0.7
dataset.idx <- 4  # y loci for the horizontal line dividing dataset from platforms

#### Variables ####
PDIR <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/TCGA'
ctypes <- list.files(file.path(PDIR, "gistic_input"), pattern="TP.txt$")
ctypes <- setNames(ctypes, gsub(".txt", "", ctypes))



#### MAIN ####
ctype_df <- lapply(ctypes, function(ctype){
  df <- read.table(file.path(PDIR, "gistic_input", ctype), sep="\t", comment.char = "", 
                   header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  del_df <- df[df$`Del q-value` <= 0.05,]
  amp_df <- df[df$`Amp q-value` <= 0.05,]
  list("AMP"=amp_df[,c('Arm', 'Amp q-value')], "DEL"=del_df[,c('Arm', 'Del q-value')])
})

amp_mat <- reduceArmMatrix(lapply(ctype_df, function(x) x$AMP))
del_mat <- reduceArmMatrix(lapply(ctype_df, function(x) x$DEL))
ad_mats <- list("AMP"=amp_mat, "DEL"=del_mat)




pdf(file.path(PDIR, "output", "GISTIC2", "amp-del_gistic.pdf"), height = 11, width=9)
split.screen(matrix(c(0, 0.75, 0.8, 1.0,
                      0, 0.75, 0, 0.8,
                      0.75, 1, 0, 0.8), byrow = TRUE, ncol=4))

n.anal <- nrow(ad_mats[[1]])
n.plat <- ncol(ad_mats[[1]])

### Plot identity matrix
lapply(names(ad_mats), function(cna){
  chr_arm_mat <- ad_mats[[cna]]
  
  if(cna == names(ad_mats)[1]) {
    screen(2); par(mar=c(1, 2, 0.5, 0.2))
    blankPlot(xlim=c(0, n.plat), ylim=c(1, n.anal))
  }
  
  sapply(1:n.anal, function(i) {
    ## Background tracks
    if(cna == names(ad_mats)[1]){
      points(x=c(1:n.plat), y=rep(i, n.plat), 
             pch=16, col="grey", cex=point.size)
    }
    ## Foreground tracks based on identity matrix
    rect.track <- which(as.logical(chr_arm_mat[i,]))
    points(x=rect.track, y=rep(i, length(rect.track)),
           pch=16, col=if(cna=='AMP') 'red' else 'blue', 
           cex=if(cna == 'AMP') point.size else point.size * 0.66)
  })
})  

### Label Chr-Arms
screen(3); par(mar=c(1, 0, 0.5, 0));
blankPlot(xlim=c(0, 10), ylim=c(1, n.anal))
anal.idx <- seq_along(rownames(ad_mats[[1]]))
ad_cnts <- paste0(" [", apply(sapply(ad_mats, rowSums), 1, paste, collapse=","), "]")
text(y=anal.idx, x=0, pos = 4,
     labels = paste0(gsub("\\.", " ", rownames(ad_mats[[1]])), ad_cnts),
     cex=text.size)
abline(v = 0, lty=1, lwd=2, col="black")

### Label Cancer Types
screen(1); par(mar=c(0, 2, 0, 0.2));
blankPlot(xlim=c(0, n.plat), ylim=c(0, 10))
plat.idx <- seq_along(colnames(ad_mats[[1]]))
text(x=plat.idx, y=0, adj=0, srt=90,
     labels = colnames(ad_mats[[1]]),
     cex=text.size)


close.screen(all.screens=TRUE)
dev.off()
