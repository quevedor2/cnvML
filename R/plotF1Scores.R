library(optparse)

option_list = list(
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="Directory containing the [GDSC/GNE/CCLE]/_F1_test.csv")
)
opt = parse_args(OptionParser(option_list=option_list))

DS <- c("GDSC", "CCLE", "GNE")
three_cols <- setNames(c('#d8b365','#e66101','#5e3c99'),
                       c('raw', 'naive', 'transfer'))
CATEGORIES = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD",
              "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
              "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
              "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
              "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
              "UCEC", "UCS", "UVM", "Normal")

if(is.na(opt$dir)){
  OUTDIR='.'
  OUTDIR='/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnvML/CCL/input/hilbert'
  OUTDIR='/cluster/home/quever/pughlab/projects/cancer_cell_lines/CCL/models/model4'
} else {
  OUTDIR <- opt$dir
}
setwd(OUTDIR)

## n-count of samples per cancer type
DATADIR <- file.path(gsub("/models.*", "", OUTDIR), "data")
ctype_cnts <- sapply(setNames(DS, DS), function(ds){
  ctype_cnt <- sapply(CATEGORIES, function(ctype) {
    length(list.files(file.path(DATADIR, ds, ctype)))
  })
  c("mean"=mean(ctype_cnt[ctype_cnt != 0]), "sd"=sd(ctype_cnt[ctype_cnt != 0]), 
    "sum"=sum(ctype_cnt), n=sum(ctype_cnt != 0))
})

## barplot of F1 scores per cancer type
pdf(file.path(OUTDIR, paste0("f1Plot_", basename(getwd()), ".pdf")), width = 10, height = 5.5)
par(mfrow=c(3,1), mar=c(4, 5, 1, 0), lwd = 0.2)
ds_pvals <- lapply(c("GDSC", "CCLE", "GNE"), function(ds){
  fls <- list.files(ds, pattern="F1_test.csv")
  f1 <- lapply(fls, function(f){
    read.csv(file.path(ds, f))
  })
  f1 <- t(sapply(f1, function(i) i$Frac))
  rownames(f1) <- gsub("_.*", "", fls)
  colnames(f1) <- CATEGORIES
  f1 <- f1[match(rownames(f1), c("raw", "naive", "transfer")),]
  sample_lim <- (nrow(f1)+1) * ncol(f1)
  bp <- barplot(f1, beside = T, col=three_cols, las=2, 
                ylim=c(0,1), xlim=c(0, (sample_lim + (4*nrow(f1)))),
                border = 'gray85',  ylab=ds, xaxt='n')
  axis(side = 1, at = bp[2,], labels = colnames(f1), las=2)
  
  ## Aggregate val and plot 
  vals <- apply(f1, 1, function(x){
    setNames(c(min(x, na.rm=TRUE), max(x, na.rm=TRUE),
               quantile(x, na.rm=TRUE,c(0.25, 0.5, 0.75)), 
               mean(x, na.rm=TRUE)), 
             c("min", "max", "Q1", "med", "Q3", "mean"))
  })
  idx <- 1
  apply(vals, 2, function(x){
    xpos <- (idx*nrow(f1)-1) + sample_lim
    segments(x0 = xpos, y0 = vals['min',idx], x1 = xpos, y1 = vals['max',idx])
    rect(xleft = xpos - 1, ybottom = vals['Q1',idx], 
         xright = xpos + 1, ytop = vals['Q3',idx], col=three_cols[idx])
    segments(x0 = xpos - 1, y0 = vals['med',idx], 
             x1 = xpos + 1, y1 = vals['med',idx])
    idx <<- idx + 1
  })
  axis(side = 1, at = sample_lim + (1:nrow(f1) * nrow(f1) - 1), 
       labels = names(three_cols), las=2)
  
  pvals <- apply(f1, 1, function(j){ 
      round(t.test(f1['transfer',],j, alternative='greater')$p.val, 4)
      #round(t.test(f1['transfer',],j, alternative='greater')$estimate,3)
  })
  return(pvals)
})
dev.off()