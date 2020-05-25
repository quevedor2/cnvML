#### Code Updates ####
######################
# May-25-2020
# Looks at the TCGA/output/hilbert directory for each cancer_type (assumes)
# that hilbert_training_data.py and place_file_in_class.sh were run on the
# seg files. This script will calculate the "n" for each cancer_type and 
# map those files back to the metadata containign ploidy/purtiy values.
# Then it will calculate a summary statistic and plot n, purity, ploidy
# in a barplot for Figure 1 of the cnvML paper

#### Viz ####
#############
## Visualization function to add "boxplot" type bars
addBars <- function(x, y0, y1){
  segments(x0 = x, y0 = y0, x1 = x, y1 = y1)
  segments(x0 = x-0.2, y0 = y0, x1 = x+0.2, y1 = y0)
  segments(x0 = x-0.2, y0 = y1, x1 = x+0.2, y1 = y1)
}

#### Variables ####
###################
PDIR='/mnt/work1/users/pughlab/projects/cancer_cell_lines/TCGA'
purity_file = file.path('input', 'TCGA_mastercalls.abs_tables_JSedit.fixed.txt')
hilbert_files = file.path('output', 'hilbert')

#### MAIN ####
##############
## Read in dataframe containg purity and ploidy per TCGA file
purity_meta <- read.table(file.path(PDIR, purity_file), header = TRUE, sep="\t",
                          stringsAsFactors = FALSE, check.names = FALSE)

## Assign TCGA files to cancer type
cancer_types <- lapply(list.files(file.path(PDIR, hilbert_files)), 
                       function(cancer_class){
  files <- list.files(file.path(PDIR, hilbert_files, cancer_class))
  if(length(files) > 0){
    data.frame('array'=gsub(".png", "", files),
               "cancer_type"=rep(cancer_class, length(files)))
  }
})
names(cancer_types) <- list.files(file.path(PDIR, hilbert_files))
cancer_types <- do.call(rbind, cancer_types)

# Split by cancer type again
tcga_meta <- merge(cancer_types, purity_meta, by="array", all.x=TRUE)
tcga_metas <- split(tcga_meta, tcga_meta$cancer_type)

## Summarize the ploidy/purity for each cancer type and 
## stick Normal class at the end
tcga_summary <- as.data.frame(t(sapply(tcga_metas, function(tm){
  if(all(is.na(tm$purity))){
    summ <- c(nrow(tm), rep(2, 3), rep(0, 3))
  } else {
    summ <- c(nrow(tm),
            quantile(tm$ploidy, c(0.25, 0.50, 0.75), na.rm=TRUE),
            quantile(tm$purity, c(0.25, 0.50, 0.75), na.rm=TRUE))
  }
  
  setNames(summ, c("n", paste0("ploidy_", c(0.25, 0.50, 0.75)),
                   paste0("purity_", c(0.25, 0.50, 0.75))))
})))
normal_idx <- grep("Normal", rownames(tcga_summary))
tcga_summary <- tcga_summary[c(c(1:nrow(tcga_summary))[-normal_idx], 
                               normal_idx),]

## Viz that stuff!
pdf(file.path(PDIR, "output", "tcga_stats.pdf"), height = 6, width=5)
par(mfrow=c(3,1), mar=c(4, 4.1, 0.5, 2.1))
ctype_labels <- rownames(tcga_summary)
ctype_cols <- c(rep("#cccccc", nrow(tcga_summary))[-1], "#4dac26")
# Barplot number of samples and quantile cutoff
barplot(tcga_summary$n, ylab="n", las=1, col=ctype_cols)
abline(h = quantile(tcga_summary$n, 0.25), lty=2)
# Barplot median ploidy and 1st/3rd quantile
bp_x <- barplot(tcga_summary$'ploidy_0.5', ylab="Ploidy", 
                las=1, ylim=c(0,4), col=ctype_cols)
addBars(x=bp_x, y0=tcga_summary$'ploidy_0.25', y1=tcga_summary$'ploidy_0.75')
# Barplot median purity and 1st/3rd quantile
bp_x <- barplot(tcga_summary$'purity_0.5', ylab="Tumor purity", 
                las=1, ylim=c(0,1), col=ctype_cols)
addBars(x=bp_x, y0=tcga_summary$'purity_0.25', y1=tcga_summary$'purity_0.75')
axis(1, at=bp_x, labels = ctype_labels, las=2)
dev.off()