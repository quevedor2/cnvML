library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(taRifx)

## Creates a circle-line plot for analysis done on NETseq
#### Functions ####
blankPlot <- function(xlim, ylim){
  plot(0, type='n', xlab='', ylab='',
       xaxt='n', yaxt='n', axes=FALSE, 
       xlim=xlim, ylim=ylim)
}

parseChrLoc <- function(x, chr_df){
  posx <- c("chr"=gsub(":.*", "", x),
                     "start"=gsub("^.*:([0-9]*)-.*", "\\1", x),
                     "end"=gsub("^.*-", "", x))
  posx['cumstart'] <- as.integer(posx['start']) + chr_df[match(posx['chr'], chr_df$chr),'cumstart']
  posx['cumend'] <- as.integer(posx['end']) + chr_df[match(posx['chr'], chr_df$chr),'cumstart']
  return(posx)
}


#### Parameters ####
line.width <- 0.055 # Size of the line going across platforms
point.size <- 2 # Size of the dots for identity matrix
text.size <- 0.7
dataset.idx <- 4  # y loci for the horizontal line dividing dataset from platforms

#### Variables ####
PDIR <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/TCGA'
ctypes <- list.files(file.path(PDIR, "gistic_input", "focal"), pattern="TP_.*.txt$")
ctypes <- split(ctypes, gsub("-TP.*", "", ctypes))


#### MAIN ####
chrs <- paste0("chr", c(1:22, "X", "Y"))
chr_sizes <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[chrs]
cum_chr_sizes <- setNames(c(1, cumsum(as.numeric(chr_sizes))+1), c(chrs, 'chrZ'))
chr_df <- data.frame("chr"=names(chr_sizes), 
                     "start"=rep(1, length(chr_sizes)), "end"=chr_sizes, 
                     "cumstart"=cum_chr_sizes[-length(cum_chr_sizes)], 
                     "cumend"=cum_chr_sizes[-1])
chr_df$cummid = with(chr_df, cumstart + (cumend - cumstart)/2)

ctype_df <- lapply(ctypes, function(ctype){
  print(ctype)
  q_thresh <- 0.05
  dfl <- lapply(ctype, function(ct){
    df <- read.table(file.path(PDIR, "gistic_input", "focal", ct), sep="\t", comment.char = "", 
                     header=FALSE, stringsAsFactors = FALSE, check.names = FALSE, nrows=4)
    rownames(df) <- df[,1]
    df <- as.data.frame(t(df[,-1]))
    if(nrow(df) > 0){
      df$cna <- gsub("^.*_([a-z]*).txt", "\\1", ct)
      df <- df[which(as.numeric(as.character(df$'q value')) <= q_thresh),]
      posdf<- do.call(rbind, lapply(df$'wide peak boundaries', parseChrLoc, chr_df=chr_df))
      df <- taRifx::remove.factors(cbind(df, posdf))
    } else { 
      df <- NA
    }
    return(df)
  })
  names(dfl) <- c('AMP', 'DEL')
  return(dfl)
})

lapply(ctype_df, function(ctype_x){
  
})

ov_fracs <- lapply(ctype_df, function(ctype_x){
  agg_dfx <- do.call(rbind, ctype_x)
  if(ncol(agg_dfx) != 1){
    if(any(is.na(agg_dfx$chr))) agg_dfx <- agg_dfx[-which(is.na(agg_dfx$chr)),]
    grx <- sort(makeGRangesFromDataFrame(agg_dfx))
    
    ov_frac <- sapply(ctype_df, function(ctype_y){
      agg_dfy <- do.call(rbind, ctype_y)
      if(ncol(agg_dfy) != 1){
        if(any(is.na(agg_dfy$chr))) agg_dfy <- agg_dfy[-which(is.na(agg_dfy$chr)),]
        gry <- sort(makeGRangesFromDataFrame(agg_dfy))
        
        int_sum <- sum(width(intersect(grx, gry)))
        x_sum <- sum(width(setdiff(grx, gry, ignore.strand=TRUE)))
        y_sum <- sum(width(setdiff(gry, grx, ignore.strand=TRUE)))
        int_sum / (x_sum + y_sum + int_sum)
      } else {
        0
      }
    })
  } else {
    ov_frac <- setNames(rep(0, length(ctype_df)), names(ctype_df))
  }
  return(ov_frac)
})
ov_dist <- abs(1-do.call(rbind, ov_fracs))
ord <- hclust(as.dist(ov_dist, diag=TRUE))$order

pdf(file.path(PDIR, "output", "GISTIC2", "amp-del_focal_gistic.pdf"), height = 11, width=9)
plot(0, type='n', ylim=c(0, length(ctype_df)+1), xlim=c(0, max(chr_df$cumend)),
     ylab='', xlab='Genomic position', axes=FALSE)
idx <<- 1
sapply(ctype_df[ord], function(ctype){
  print(names(ctype_df)[ord][idx])
  sapply(names(ctype), function(ct){
    col <- if(ct == 'AMP') 'red' else 'blue'
    if(class(ctype[[ct]]) == 'data.frame'){
      start_pos <- as.numeric(ctype[[ct]]$cumstart)
      end_pos<- as.numeric(ctype[[ct]]$cumend)
      if(length(start_pos) > 0){
        rect(xleft = start_pos, ybottom = idx - if(ct == 'DEL') 0.4 else 0, 
             xright = end_pos, ytop = idx + if(ct == 'AMP') 0.4 else 0, 
             col = col, border = col)
        
      }
    }
  })
  idx <<- idx + 1
})

abline(h = seq_along(ctype_df), col="grey")
abline(v = c(1, chr_df$cumend), col="grey", lty=2)
axis(side = 2, at=seq_along(ctype_df), labels = names(ctype_df[ord]), 
     las=1, tick = FALSE, cex.axis=0.7)
axis(side = 1, at=chr_df$cummid, labels = rownames(chr_df), 
     las=2, tick = FALSE, cex.axis=0.7)
axis(side=1, at=c(1, chr_df$cumend), labels=rep("", nrow(chr_df)+1), tick=TRUE)

dev.off()
