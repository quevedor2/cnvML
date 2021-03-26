#boxplot_tcgaF1
library(optparse)

option_list = list(
  make_option(c("-g", "--gene"), action="store", default=NA, type='character',
              help="Directory containing the GENES TCGA gene F1 csv"),
  make_option(c("-b", "--bin"), action="store", default=NA, type='character',
              help="Directory containing the BINS TCGA gene F1 csv"),
  make_option(c("-h", "--hilbert"), action="store", default=NA, type='character',
              help="Directory containing the HILBERT TCGA gene F1 csv")
)
opt = parse_args(OptionParser(option_list=option_list))

PDIR <-'/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnvML/TCGA'
cols <- c('gray34','gray85')
hilbert_col <- '#4393c3'
CATEGORIES = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD",
               "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
               "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
               "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
               "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
               "UCEC", "UCS", "UVM", "Normal")

if(any(sapply(opt, is.na))){
  paths <- list("bin"=file.path(PDIR, "input", 'f1_boxplot', 'tcga_bin'),
                "gene"=file.path(PDIR, "input", 'f1_boxplot', 'tcga_genes'),
                "hilbert"=file.path(PDIR, "input", 'f1_boxplot', 'tcga_hilbert'))
  #stop("Need to pass in directories")
} else {
  paths <- list("bin"=opt$bin,
                "gene"=opt$gene,
                "hilbert"=opt$hilbert)
}

all_f1s <- lapply(names(paths), function(p){
  if(p == 'hilbert'){
    f1 <- read.csv(file.path(paths[[p]], "F1_test.csv"))
    f1 <- cbind(f1$Frac, rep(NA, nrow(f1)))
    colnames(f1) <- c('CNN', '')
  } else {
    dirs <- c("ann", "lr")
    f1 <- lapply(dirs, function(f){
      read.csv(file.path(paths[[p]], f, "F1_test.csv"))
    })
    f1 <- sapply(f1, function(i) i$Frac)
    colnames(f1) <- toupper(dirs)
  }
  rownames(f1) <- CATEGORIES
  if(any(is.nan(f1))) f1[is.nan(f1)] <- NA
  return(f1)
})
names(all_f1s) <- names(paths)

##p-values for Fig 2.a
ann_lr_pval <- sapply(all_f1s[1:2], function(x) { t.test(x[,1], x[,2], alternative='greater')$p.val})
# bin       gene 
# 0.01235437 0.06726988 
dl_mat <- Reduce(function(x,y) cbind(x,y), lapply(all_f1s, function(x) x[,1, drop=FALSE]))
colnames(dl_mat) <- names(all_f1s)
cnn_ann_pval <- apply(dl_mat,2, function(x){apply(dl_mat, 2, function(y) t.test(x,y, alternative = 'greater')$p.val)})
# bin      gene   hilbert
# bin     0.5000000 0.5155623 0.6397957
# gene    0.4844377 0.5000000 0.6354212
# hilbert 0.3602043 0.3645788 0.5000000


pdf(file.path(PDIR, "output", "f1_boxplot", "tcga_F1.pdf"), width = 5, height = 3)
par(mfrow=c(4,1), mar=c(0, 5.1, 0, 2))
sapply(names(all_f1s), function(f1_id){
  f1_tmp <- all_f1s[[f1_id]]
  boxplot(f1_tmp, horizontal = TRUE, ylim=c(0,1), col=cols, axes=FALSE)
  axis(side = 2, at=c(1,2), labels=colnames(all_f1s[[f1_id]]), las=2)
  axis(side = c(1,3), at=c(-1, 2), lty = 3, labels=c('', ''), lwd=0.5)
  if(f1_id=='hilbert') axis(side = 1, at=seq(0, 1, by=0.2), labels=seq(0, 1, by=0.2))
})
dev.off()

##Fig 2.b
pdf(file.path(PDIR, "output", "f1_boxplot", "tcga_hilbertF1.pdf"), width = 9, height = 5)
barplot(all_f1s[['hilbert']][,1], ylim=c(0,1), las=2, col=hilbert_col, 
        border=NA, ylab="F1-score")
dev.off()


