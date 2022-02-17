#boxplot_tcgaF1
library(optparse)
library(reshape2)
library(ggplot2)
library(cowplot)

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

ds <- c('CCLE', 'GDSC')
# ds <- 'CCLE'
PDIR <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/cnvML'
PDIR <- '/cluster/projects/pughlab/projects/cancer_cell_lines'
if(opt$cluster!='h4h'){
  PDIR <- if(nchar(ds) > 0) file.path(PDIR, 'CCL') else file.path(PDIR, 'TCGA')
  tcga <- if(nchar(ds) > 0) '' else 'tcga_'
}

cols <- c('gray34','gray85')
hilbert_col <- '#4393c3'
CATEGORIES = setNames(c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD",
                        "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
                        "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
                        "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
                        "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
                        "UCEC", "UCS", "UVM", "Normal"),
                      as.character(c(0:33)))

if(any(sapply(opt, is.na))){
  dataset <- 'CCL'
  
  if(opt$cluster == 'h4h'){
    paths <- list("hilbert_ascn"=file.path(PDIR, dataset, "models", "hilbert", "ASCN", "alexnet"),
                  "hilbert_tcn"=file.path(PDIR, dataset, "models", "hilbert", "TCN", "alexnet"),
                  "morton_ascn"=file.path(PDIR, dataset, "models", "morton", "ASCN", "alexnet"),
                  "diagonal_ascn"=file.path(PDIR, dataset, "models", "diagonal", "ASCN", "alexnet"),
                  "sweep_ascn"=file.path(PDIR, dataset, "models", "sweep", "ASCN", "alexnet"),
                  "scan_ascn"=file.path(PDIR, dataset, "models", "scan", "ASCN", "alexnet"),
                  "random_ascn"=file.path(PDIR, dataset, "models", "random", "ASCN", "alexnet"))
    outdir <- file.path(PDIR, paste0(dataset, "_plots"), "f1_boxplot")
    #stop("Need to pass in directories")
  } else if(opt$cluster == 'mor'){
    # Removed
    #stop("Need to pass in directories")
    outdir <- file.path(PDIR, "output", ds, "f1_boxplot")
  } else {
    stop("No cluster specified")
  }
} else {
  #Removed
}

all_f1s <- lapply(names(paths), function(p){
  print(p)
  if(grepl('tcn|ascn', p)){
    f1s <- lapply(ds, function(ds_i){
      path_i <- paths[[p]]
      path_i <- gsub("/([a-zA-Z]*$)", paste0("/", ds_i, "/\\1"), path_i)
      f1 <- tryCatch({
        f1 <- setNames(read.csv(file.path(path_i, "F1_test.csv"))[,c(1:2)],
                       c('Ctype', 'F1'))
        ids <- strsplit(path_i, split="/")[[1]]
        ids <- ids[(length(ids)-3):length(ids)]
        f1$Ctype <- CATEGORIES[as.character(f1$Ctype)]
        f1$model <- ids[1]
        f1$CN <- ids[2]
        f1$DS <- ds_i
        return(f1)
      }, error=function(e){NULL})
      return(f1)
    })
  } else {
    dirs <- c("ann", "lr")
    print("Bin or Gene...")
    f1 <- tryCatch({
      f1 <- lapply(dirs, function(f){
        read.csv(file.path(paths[[p]], f, "F1_test.csv"))
      })
      f1 <- sapply(f1, function(i) i$Frac)
      colnames(f1) <- toupper(dirs)
      return(f1)
    }, error=function(e){NULL})
  }
  return(f1s)
})
all_f1s <- unlist(all_f1s, recursive = F)
m_f1 <- do.call(rbind, all_f1s)

##p-values for Fig 2.a
# # Test whether ANN ([,1]) for gene/bins are more predictive than LR ([,2])
# ann_lr_pval <- tryCatch({
#   sapply(all_f1s[1:2], function(x) { t.test(x[,1], x[,2], alternative='greater')$p.val})
# }, error=function(e){NULL})

## Test whether the performance of one metric is better than another
# x <- matrix(c(11:20, 1:10, 1:10), ncol=3)
# apply(x, 2, function(i){
#   apply(x, 2, function(j) t.test(i,j, alternative='greater')$p.val)
# }) %>% round(., 3)

m_f1_ds <- split(m_f1, m_f1$DS)
cnn_ann_pvals <- lapply(m_f1_ds, function(m_f1_ds_i){
  f1_i <- split(m_f1_ds_i, with(m_f1_ds_i, paste(model, CN, sep="_")))
  dl_mat <- Reduce(function(x,y) merge(x,y,by='Ctype'), lapply(f1_i, function(i) i[,c('Ctype', 'F1')]))
  rownames(dl_mat) <- dl_mat[,1]
  dl_mat <- dl_mat[,-1]
  colnames(dl_mat) <- names(f1_i)
  cnn_ann_pval <- apply(dl_mat,2, function(x){
    apply(dl_mat, 2, function(y) t.test(x,y, alternative = 'greater')$p.val)
  })
  cnn_ann_pval
})

########################
#### Visualizations ####
dir.create(file.path(outdir), showWarnings = F, recursive = T)
## Fig 2.X
# Heatmap of t.test pvalues for model vs model performance comparison per transfer dataset
b <- c(0, 0.05, 0.1, 0.5, 1)
ggs <- lapply(names(cnn_ann_pvals), function(ds_i){
  cnn_pval_dsi <- cnn_ann_pvals[[ds_i]]
  ggplot(melt(cnn_pval_dsi), aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    scale_fill_gradientn(limits = c(min(b),max(b)),
                         colours=colorRampPalette(c("red", "yellow", "grey", "black"))(7),
                         breaks=b, labels=format(b),
                         values=c(0, seq(0.01, 0.5, length.out=5), 1)) +
    theme_classic() +
    xlab("Model_j") + ylab("Model_i") +
    ggtitle(paste0(ds_i, ": One-sided greater than t-test p-value \n 
                   for Model-i compared to Model-j")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
})
pdf(file.path(outdir, "ccl_performanceF1.pdf"), width = 10)
plot_grid(plotlist = ggs, ncol=2)
dev.off()


#out_pdf <- file.path(PDIR, "output", "f1_boxplot", "tcga_F1-allSFC.pdf")
##Fig 2.a
# Boxplot comparing the performance metrics across all the different models
pdf(file.path(outdir, "ccl_F1.pdf"), width = 5, height = 7)
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

