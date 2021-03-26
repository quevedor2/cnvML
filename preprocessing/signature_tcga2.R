###################
#### Functions ####
preprocessSegs <- function(segd, atype='pancan', 
                           balance=TRUE, quantile_cutoff=0.25,
                           min_n=FALSE, hard_cutoff=NULL){
  segd_by_sample <- split(segd, segd$Sample)
  ctype_cnts <- table(sapply(segd_by_sample, function(i) unique(i$ctype)))
  if(!is.null(hard_cutoff)){
    sample_to <- hard_cutoff
  } else {
    sample_to <- quantile(ctype_cnts, quantile_cutoff) ## Number of samples per cancer type
  }
  
  
  if(atype=='pancan'){
    if(balanced){
      print("Balanced cancer type representation...")
      
      segd_ctype <- split(segd, segd$ctype) ## Find each sample for each ctype
      segd_sample_ctype <- lapply(segd_ctype, function(s){
        names(split(s, s$Sample))
      })
      
      if(min_n){
        rm_idx <- sapply(segd_sample_ctype, length) < (sample_to / 4)
        print(paste0("Removing ", paste(names(segd_sample_ctype)[which(rm_idx)], collapse=",")))
        if(any(rm_idx)) segd_sample_ctype <- segd_sample_ctype[-which(rm_idx)]
      }
      
      balanced_samples <- lapply(segd_sample_ctype, function(s){
        if(length(s) >= sample_to){
          ## Undersample from cancer types with a lot of samples
          print(paste0("Undersampling ", unique(segd_by_sample[[s[1]]]$ctype)))
          segd_by_sample[sample(s, size = sample_to, replace= F)]
        } else {
          ## Sample to fill and create unique IDs
          print(paste0("Oversampling ", unique(segd_by_sample[[s[1]]]$ctype)))
          samples <- c(s, sample(s, size=(sample_to - length(s)), replace=T))
          while(any(duplicated(samples))){
            samples[duplicated(samples)] <- paste0(samples[duplicated(samples)], "X")
          }
          
          ## Change sample names in duplicates to match unique ID
          lapply(samples, function(sid){
            x <- segd_by_sample[[gsub("X*$", "", sid)]]
            x$Sample <- sid
            return(x)
          })
        }
      }) # Sample to balance
      segd2 <- do.call(rbind, lapply(balanced_samples, function(i) do.call(rbind,i)))
      
      segd2_by_sample <- split(segd2, segd2$Sample)
      ctype_cnts <- table(sapply(segd2_by_sample, function(i) unique(i$ctype)))
      
      print("Cancer type representation...")
      print(ctype_cnts)
      segd_ctype <- list("PANCAN"=segd2)
    } else {
      print("Cancer type representation...")
      print(ctype_cnts)
      segd_ctype <- list("PANCAN"=segd)
    }
  } else {
    print("No balancing, separating cancers by cancer type...")
    segd_ctype <- split(segd, segd$ctype)
    if(min_n){
      # Removing samples:
      rm_idx <- ctype_cnts < (sample_to / 4)
      rm_ctypes <- names(ctype_cnts)[which(rm_idx)]
      print(paste0("Removing ", paste(rm_ctypes, collapse=",")))
      if(any(rm_idx)) segd_ctype <- segd_ctype[-which(names(segd_ctype) %in% rm_ctypes)]
    }
  }
  return(segd_ctype)
}

formatSeg <- function(segd, analysis, PDIR=NULL){
  if(analysis=='TCGA'){
    colnames(segd)[2:4] <- c('chrom', 'loc.start', 'loc.end')
    
    normal_segf <- file.path(PDIR, '1000G', 'eacon', 'symlinks', 'seg', '1000g.seg')
    normal_segd <- read.table(normal_segf, sep="\t", header=TRUE, stringsAsFactors = FALSE)
    colnames(normal_segd) <- c('Sample', 'chrom', 'Chrom', 'loc.start', 'loc.end', 
                               'Length', 'Modal_Total_CN', 'Modal_HSCN_1', 'Modal_HSCN_2')
    segd2 <- plyr::rbind.fill(segd, normal_segd)
    segd <- segd2
  } else {
    colnames(segd)[c(1:3,6,11,7,8)] <- c('chrom', 'loc.start', 'loc.end',
                                         'Sample', 'Modal_Total_CN', 
                                         'Modal_HSCN_1', 'Modal_HSCN_2')
    if(any(segd$chrom %in% c('chrX', 'chrY'))){
      segd <- segd[-which(segd$chrom %in% c('chrX', 'chrY')),]
    }
  }
  
  if(sum(grepl("chrom", colnames(segd), ignore.case=TRUE)) > 1){
    segd <- segd[,-grep("Chrom", colnames(segd))]
  }
  if(any(is.na(segd$chrom))) segd <- segd[-which(is.na(segd$chrom)),]
  segd$chrom <- gsub("23", "X", segd$chrom) %>% gsub("24", "Y", .)
  return(segd)
}

##############
#### Main ####
#### 0) Parameters/Libraries  ####
library(sigminer)
library(NMF)
atype <- 'pancan'
atype <- 'ctype'
analysis <- 'TCGA'
balanced=TRUE

PDIR <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines'
if(analysis=='TCGA'){
  seg_files <- "TCGA_mastercalls.abs_segtabs.fixed.txt"
  META <- file.path(PDIR, analysis, 'input/merged_sample_quality_annotations.trimmed.tsv')
  anno <- read.table(META, header=TRUE, check.names = FALSE, 
                     sep="\t", stringsAsFactors = FALSE)
} else if(analysis=='ccl_aggregate'){
  seg_files <- c('CCLE_cna_hg19.seg', 'GDSC_cna_hg19.seg', 'gCSI_cna_hg19.seg')
  META <- file.path(PDIR, analysis, "ref", "onco_meta_df.rds")
  anno <- readRDS(META)
  anno$GDSC <- gsub(".cel$", "", anno$GDSC, ignore.case = TRUE)
  colnames(anno) <- gsub("oncocode", "cancer type", colnames(anno))
  colnames(anno) <- gsub("^GNE$", "gCSI", colnames(anno))
} else {
  stop("Analysis must be either 'TCGA' or 'ccl_aggregate'")
}



segf <- file.path(PDIR, analysis, "input", seg_files)
segd <- read.table(segf, sep="\t", header=TRUE, stringsAsFactors = FALSE)
sample_anno <- data.frame("sample"=unique(segd$Sample))
sample_anno$ctype <- sapply(sample_anno$sample, function(s){
  anno[grep(s, anno$aliquot_barcode)[1],]$'cancer type'
})
s_id <- setNames(sample_anno$ctype, sample_anno$sample)

#### 1) Address class imbalance - balance cancer types  ####
segd$ctype <- s_id[segd$Sample]
segd_ctype <- preprocessSegs(segd, atype, balanced, 0.50, min_n=TRUE)
# length(unique(segd_ctype[[1]]$Sample)) ## 4554


#### 2) Parse seg data into Sigminer:CN counts ####
cn_ctype <- lapply(segd_ctype, function(segd){
  print(unique(segd$ctype))
  cn <- sigminer::read_copynumber(segd, samp_col = "Sample", max_copynumber = 20L,
                                  seg_cols = c("Chromosome", "Start", "End", "Modal_Total_CN"),
                                  genome_build = "hg19", complement = FALSE, verbose = TRUE)
  return(cn)
})
save(cn_ctype, file=file.path(PDIR, 'TCGA', "output", "signatures",
                              paste0("cn_ctype_", atype, 
                                     "-", balanced, ".rda")))

#### 3) Convert CN counts into the feature-set ####
load(file.path(PDIR, 'TCGA', "output", "signatures", 
               paste0("cn_ctype_", atype, "-", balanced, ".rda"))) #cn_ctype
cn_tally_ctype <- lapply(cn_ctype, function(cn){
  cn_tally_W <- sig_tally(cn, method = "W")
  return(cn_tally_W)
})
save(cn_tally_ctype, file=file.path(PDIR, 'TCGA', "output", "signatures",
                                    paste0("cn_tally_ctype_",
                                    atype, "-", balanced, ".rda")))

#### 4) Deconvolute Feature-set matrix using NMF ####
library(NMF)
## ---- Bayesian NMF "Auto" version
load(file.path(PDIR, 'TCGA', "output", "signatures",
               paste0("cn_tally_ctype_",
                      atype, "-", balanced, ".rda"))) #cn_tally_ctype
sig_auto_ctype <- lapply(cn_tally_ctype, function(cn_tally_W){
  sig_auto_extract(cn_tally_W$nmf_matrix, cores=1, nrun = 5, K0=30,
                   method='L1W.L2H', strategy='stable', destdir = file.path(PDIR, "tmp"))
})
save(sig_auto_ctype, file=file.path(PDIR, 'TCGA', "output", "signatures",
                                    paste0("cn_autosig_ctype_", 
                                    atype, "-", balanced, ".rda")))

## ---- Manual method to select optimal K
load(file.path(PDIR, 'TCGA', "output", "signatures",
               paste0("cn_tally_ctype_",
                      atype, "-", balanced, ".rda"))) #cn_tally_ctype
cnt <- 1
sig_w_ctype <- lapply(cn_tally_ctype, function(cn_tally_W){
  print(paste0(">> NMF: ", cnt))
  cnt <<- cnt + 1
  sig_estimate(cn_tally_W$nmf_matrix, range=2:30, cores=2, keep_nmfObj=TRUE, use_random=TRUE,
               nrun = 30, verbose=TRUE, save_plots=TRUE, pConstant = 1e-13,
               plot_basename=file.path(PDIR, 'TCGA', "output", "signatures", "nmf"))
})
save(sig_w_ctype, file=file.path(PDIR, 'TCGA', "output", "signatures",
                                 paste0("cn_nmf_ctype_", 
                                        atype, "-", balanced, ".rda")))

for(nm in names(sig_w_ctype)){
  pdf(file.path(PDIR, 'TCGA', "output", "signatures", "nmf", paste0("ctypes_", nm, ".pdf")), width=12)
  nmf_dat <- sig_w_ctype[[nm]]$nmfEstimate
  plot(nmf_dat)
  dev.off()
}

#### 5) Visualize CN signatures ####
sig_method='auto' # 'auto' or 'ctype'
sig_type <- switch(sig_method,
                   auto='cn_autosig_ctype_',
                   ctype='cn_nmf_ctype_')
sig_obj_id <- load(file.path(PDIR, 'TCGA', "output", "signatures",
                             paste0(sig_type, atype, "-", balanced, ".rda")))
load(file.path(PDIR, 'TCGA', "output", "signatures",
               paste0(sig_type, atype, "-", balanced, ".rda")))
cn_sig <- switch(sig_method,
                 auto=sig_auto_ctype,
                 ctype=sig_w_ctype)

# Extract the Exposure (H) and Signature (W) matrices from NMF
sig_W <- sig_signature(cn_sig[[1]])
sig_H <- sig_exposure(cn_sig[[1]])
sig_Hdf <- get_sig_exposure(cn_sig[[1]])

## View the exposure and signature profiles
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

exposure_samples <- gsub("X*$", "", colnames(cn_sig[[1]]$Exposure))
exp_idx <- sapply(exposure_samples, function(i) grep(paste0("^", i), x=anno$aliquot_barcode)[1])
anno_grps <- anno[exp_idx,]$`cancer type`

pdf(file.path(PDIR, 'TCGA', "output", "signatures",
              paste0(sig_type, "TCGA-exposure.pdf")), width = 60, height = 8)
show_sig_exposure(cn_sig[[1]], palette=sample(col_vector, cn_sig[[1]]$K), 
                  sig_names=gsub("Sig", "", rownames(cn_sig[[1]]$Exposure)),
                  rm_panel_border=TRUE, rm_space=TRUE,
                  groups = anno_grps,
                  grp_order = sort(unique(anno_grps)))
dev.off()

pdf(file.path(PDIR, 'TCGA', "output", "signatures",
              paste0(sig_type, "TCGA-profile.pdf")), width = 10, height = 10)
show_sig_profile(sig_W,
                 mode = "copynumber",
                 normalize = "feature",
                 method = "W",
                 style = "cosmic")
dev.off()

#### 6) Fit 1000Genome to CN-signatures ####
sig_method='auto' # 'auto' or 'ctype'
sig_type <- switch(sig_method,
                   auto='cn_autosig_ctype_',
                   ctype='cn_nmf_ctype_')
sig_obj_id <- load(file.path(PDIR, 'TCGA', "output", "signatures",
                             paste0(sig_type, atype, "-", balanced, ".rda")))
load(file.path(PDIR, 'TCGA', "output", "signatures",
               paste0(sig_type, atype, "-", balanced, ".rda")))
cn_sig <- switch(sig_method,
                 auto=sig_auto_ctype,
                 ctype=sig_w_ctype)

segd2_wN <- formatSeg(segd = segd, analysis = 'TCGA', PDIR = PDIR)
max_sample <- 100
seg_wN_l <- split(segd2_wN, segd2_wN$Sample)
spl <- data.frame("start"=seq(1, length(seg_wN_l), by=max_sample),
                  "end"=c(seq(1, length(seg_wN_l), by=max_sample)[-1]-1, length(seg_wN_l)))
all_fits <- apply(spl, 1, function(i){
  print(paste0("Samples: ", i['start'], "-", i['end']))
  seg_wN_subset <- do.call(rbind, seg_wN_l[i['start']:i['end']])
  
  cn_wN <- sigminer::read_copynumber(seg_wN_subset, samp_col = "Sample", max_copynumber = 20L,
                                     seg_cols = c("chrom", "loc.start", "loc.end", "Modal_Total_CN"),
                                     genome_build = "hg19", complement = FALSE, verbose = TRUE)
  
  cn_tally_wN <- sig_tally(cn_wN, method = "W")
  
  sfit <- sig_fit(t(cn_tally_wN$nmf_matrix), sig = cn_sig[[1]], return_class = "data.table")
  return(sfit)
})

sig_mat <- do.call(rbind, lapply(all_fits, function(i) as.data.frame(i)[,-1]))
rownames(sig_mat) <- unlist(sapply(all_fits, function(i) as.character(i$sample)))
save(sig_mat, file=file.path(PDIR, 'TCGA', "output", "signatures", 'tcga-1000g_sigmat.rda'))

#### 7) Assign cancer types ####
load(file.path(PDIR, 'TCGA', "output", "signatures", 'tcga-1000g_sigmat.rda')) #sig_mat
ctypes <- s_id[rownames(sig_mat)] # (From Part 1))
ctypes[grep("^HG0", rownames(sig_mat))] <- 'Normal'
sig_mat$cancer_type <- ctypes

write.csv(sig_mat, file=file.path(PDIR, analysis, "output", "signatures", 
                                  paste0(atype, "_sig_matrix.csv")),
          quote = FALSE, row.names = TRUE, col.names = TRUE)

