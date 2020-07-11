library(optparse)

option_list = list(
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="Directory containing the [GDSC/GNE/CCLE]/_F1_test.csv")
)
opt = parse_args(OptionParser(option_list=option_list))

three_cols <- c('#d8b365','#e66101','#5e3c99')
CATEGORIES = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD",
              "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
              "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
              "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
              "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
              "UCEC", "UCS", "UVM", "Normal")

if(is.na(opt$dir)){
  OUTDIR='.'
} else {
  OUTDIR <- opt$dir
}
setwd(OUTDIR)
pdf(file.path(OUTDIR, paste0("f1Plot_", basename(getwd()), ".pdf")), width = 10, height = 5.5)
par(mfrow=c(3,1), mar=c(4, 5, 1, 0), lwd = 0.2)
for(ds in c("GDSC", "CCLE", "GNE")){
  fls <- list.files(ds, pattern="F1_test.csv")
  f1 <- lapply(fls, function(f){
    read.csv(file.path(ds, f))
  })
  f1 <- t(sapply(f1, function(i) i$Frac))
  rownames(f1) <- gsub("_.*", "", fls)
  colnames(f1) <- CATEGORIES
  f1 <- f1[match(rownames(f1), c("raw", "naive", "transfer")),]
  bp <- barplot(f1, beside = T, col=three_cols, las=2, ylim=c(0,1), 
                border = 'gray85',  ylab=ds, xaxt='n')
  axis(side = 1, at = bp[2,], labels = colnames(f1), las=2)
}
dev.off()