#boxplot_tcgaF1
library(optparse)
library(reshape2)
library(ggplot2)

option_list = list(
  make_option(c("-g", "--gene"), action="store", default=NA, type='character',
              help="Directory containing the GENES TCGA gene F1 csv"),
  make_option(c("-b", "--bin"), action="store", default=NA, type='character',
              help="Directory containing the BINS TCGA gene F1 csv"),
  make_option(c("-h", "--hilbert"), action="store", default=NA, type='character',
              help="Directory containing the HILBERT TCGA gene F1 csv"),
  make_option(c("-c", "--cluster"), action="store", default='h4h', type='character',
              help="Custom parameter whether working on Mordor or H4H")
)
opt = parse_args(OptionParser(option_list=option_list))

ds <- ''
# ds <- 'CCLE'
PDIR <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnvML'
PDIR <- '/cluster/projects/pughlab/projects/cancer_cell_lines'
if(opt$cluster!='h4h'){
  PDIR <- if(nchar(ds) > 0) file.path(PDIR, 'CCL') else file.path(PDIR, 'TCGA')
  tcga <- if(nchar(ds) > 0) '' else 'tcga_'
}

cols <- c('gray34','gray85')
hilbert_col <- '#4393c3'
CATEGORIES = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD",
               "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
               "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
               "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
               "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
               "UCEC", "UCS", "UVM", "Normal")

if(any(sapply(opt, is.na))){
  if(opt$cluster == 'h4h'){
    paths <- list("bin"=file.path(PDIR, "TCGA_bin", "models"),
                  "gene"=file.path(PDIR, "TCGA_genes", "models"),
                  "hilbert_ascn"=file.path(PDIR, "TCGA", "models", "hilbert", "ASCN", "alexnet"),
                  "hilbert_tcn"=file.path(PDIR, "TCGA", "models", "hilbert", "TCN", "alexnet"),
                  "morton_ascn"=file.path(PDIR, "TCGA", "models", "morton", "ASCN", "alexnet"),
                  "diagonal_ascn"=file.path(PDIR, "TCGA", "models", "diagonal", "ASCN", "alexnet"),
                  "sweep_ascn"=file.path(PDIR, "TCGA", "models", "sweep", "ASCN", "alexnet"),
                  "scan_ascn"=file.path(PDIR, "TCGA", "models", "scan", "ASCN", "alexnet"),
                  "random_ascn"=file.path(PDIR, "TCGA", "models", "random", "ASCN", "alexnet"))
    outdir <- file.path(PDIR, "TCGA_plots", "f1_boxplot")
    #stop("Need to pass in directories")
  } else if(opt$cluster == 'mor'){
    paths <- list("bin"=file.path(PDIR, "input", 'f1_boxplot', ds, 'tcga_bin'),
                  "gene"=file.path(PDIR, "input", 'f1_boxplot', ds, 'tcga_genes'),
                  "hilbert_ascn"=file.path(PDIR, "input", 'f1_boxplot', ds, paste0(tcga, 'hilbert_ascn')),
                  "hilbert_tcn"=file.path(PDIR, "input", 'f1_boxplot', ds, paste0(tcga, 'hilbert_tcn')),
                  "morton_ascn"=file.path(PDIR, "input", 'f1_boxplot', ds, paste0(tcga, 'morton_ascn')),
                  "diagonal_ascn"=file.path(PDIR, "input", 'f1_boxplot', ds, paste0(tcga, 'diagonal_ascn')),
                  "sweep_ascn"=file.path(PDIR, "input", 'f1_boxplot', ds, paste0(tcga, 'sweep_ascn')),
                  "scan_ascn"=file.path(PDIR, "input", 'f1_boxplot', ds, paste0(tcga, 'scan_ascn')),
                  "random_ascn"=file.path(PDIR, "input", 'f1_boxplot', ds, paste0(tcga, 'random_ascn')))
    #stop("Need to pass in directories")
    outdir <- file.path(PDIR, "output", ds, "f1_boxplot")
  } else {
    stop("No cluster specified")
  }
} else {
  paths <- list("bin"=opt$bin,
                "gene"=opt$gene,
                "hilbert"=opt$hilbert)
}

all_f1s <- lapply(names(paths), function(p){
  print(p)
  if(grepl('tcn|ascn', p)){
    f1 <- tryCatch({
      f1 <- read.csv(file.path(paths[[p]], "F1_test.csv"))
      f1 <- cbind(f1$Frac, rep(NA, nrow(f1)))
      colnames(f1) <- c('CNN', '')
      return(f1)
    }, error=function(e){NULL})
  } else {
    dirs <- c("ann", "lr")
    f1 <- tryCatch({
      f1 <- lapply(dirs, function(f){
        read.csv(file.path(paths[[p]], f, "F1_test.csv"))
      })
      f1 <- sapply(f1, function(i) i$Frac)
      colnames(f1) <- toupper(dirs)
      return(f1)
    }, error=function(e){NULL})
  }
  if(!is.null(f1)){
    rownames(f1) <- CATEGORIES
    if(any(is.nan(f1))) f1[is.nan(f1)] <- NA
  }
  return(f1)
})
names(all_f1s) <- names(paths)
null_idx <- sapply(all_f1s, is.null)

##p-values for Fig 2.a
# Test whether ANN ([,1]) for gene/bins are more predictive than LR ([,2])
ann_lr_pval <- tryCatch({
  sapply(all_f1s[1:2], function(x) { t.test(x[,1], x[,2], alternative='greater')$p.val})
}, error=function(e){NULL})
# bin       gene 
# 0.01235437 0.06726988 

## Test whether the performance of one metric is better than another
# x <- matrix(c(11:20, 1:10, 1:10), ncol=3)
# apply(x, 2, function(i){
#   apply(x, 2, function(j) t.test(i,j, alternative='greater')$p.val)
# }) %>% round(., 3)
dl_mat <- Reduce(function(x,y) cbind(x,y), lapply(all_f1s, function(x) x[,1, drop=FALSE]))
colnames(dl_mat) <- names(all_f1s)[which(!null_idx)]
cnn_ann_pval <- apply(dl_mat,2, function(x){
  apply(dl_mat, 2, function(y) t.test(x,y, alternative = 'greater')$p.val)
})
#cnn_ann_pval
#               bin      gene   hilbert
# bin     0.5000000 0.5155623 0.6397957
# gene    0.4844377 0.5000000 0.6354212
# hilbert 0.3602043 0.3645788 0.5000000

dir.create(file.path(outdir), showWarnings = F, recursive = T)
#out_pdf <- file.path(PDIR, "output", "f1_boxplot", "tcga_F1-allSFC.pdf")
##Fig 2.a
# Boxplot comparing the performance metrics across all the different models
pdf(file.path(outdir, "tcga_F1.pdf"), width = 5, height = 7)
m_f1s <- as.data.frame(melt(all_f1s))
m_f1s <- m_f1s[-which(m_f1s$Var2==''),]
m_f1s$L1 <- factor(m_f1s$L1, levels=unique(m_f1s$L1))
ggplot(m_f1s, aes(x=value, y=L1, fill=Var2)) +
  facet_grid(rows=vars(L1), scales = 'free_y') +
  geom_boxplot() +
  theme_classic() +
  xlim(0,1) + ylab("") + xlab("F1") +
  scale_fill_discrete(name = "Algorithm") +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  )

par(mfrow=c(10,1), mar=c(0, 8, 0, 2))
sapply(names(all_f1s), function(f1_id){
  f1_tmp <- all_f1s[[f1_id]]
  boxplot(f1_tmp, horizontal = TRUE, ylim=c(0,1), col=cols, axes=FALSE)
  axis(side = 2, at=c(1,2), labels=colnames(all_f1s[[f1_id]]), las=2, pos=0, tick = FALSE)
  axis(side = 2, at=1.5, labels=f1_id, las=2, pos=-0.1, tick = FALSE)
  axis(side = c(1,3), at=c(-1, 2), lty = 3, labels=c('', ''), lwd=0.5)
  if(f1_id=='hilbert') axis(side = 1, at=seq(0, 1, by=0.2), labels=seq(0, 1, by=0.2))
})
dev.off()

##Fig 2.b
# Barplot of TCGA hilbert ASCN cancer type prediction
pdf(file.path(outdir, "tcga_hilbertF1.pdf"), width = 9, height = 5)
barplot(all_f1s[['hilbert_ascn']][,1], ylim=c(0,1), las=2, col=hilbert_col, 
        border=NA, ylab="F1-score")
dev.off()
pdf(file.path(outdir, "tcga_hilbertF1_all.pdf"), width = 9, height = 15)
m_f1s$L2 <- with(m_f1s, paste(L1, Var2, sep="_"))
m_f1s$L2 <- factor(m_f1s$L2, levels=unique(m_f1s$L2))
m_f1s$Var1 <- factor(m_f1s$Var1, levels=unique(m_f1s$Var1))

ggplot(m_f1s, aes(x=Var1, y=value, fill=Var2)) +
  facet_grid(rows=vars(L2), scales = 'free_y') +
  geom_bar(stat='identity') +
  theme_classic() +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  ylim(0,1) + ylab("F1") + xlab("TCGA cancer types") +
  scale_fill_discrete(name = "Algorithm")
dev.off()

## Fig 2.X
# Heatmap of t.test pvalues for model vs model performance comparison
b <- c(0, 0.05, 0.1, 0.5, 1)
pdf(file.path(outdir, "tcga_performanceF1.pdf"))
ggplot(melt(cnn_ann_pval), aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradientn(limits = c(min(b),max(b)),
                       colours=colorRampPalette(c("red", "yellow", "grey", "black"))(7),
                       breaks=b, labels=format(b),
                       values=c(0, seq(0.01, 0.5, length.out=5), 1)) +
  theme_classic() +
  xlab("Model_j") + ylab("Model_i") +
  ggtitle("One-sided greater than t-test p-value for Model-i compared to Model-j") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
