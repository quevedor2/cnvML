###################
#### Functions ####
preprocessSegs <- function(segd, analysis='pancan', 
                           balance=TRUE, quantile_cutoff=0.25,
                           min_n=FALSE){
  if(analysis=='pancan'){
    segd_by_sample <- split(segd, segd$Sample)
    ctype_cnts <- table(sapply(segd_by_sample, function(i) unique(i$ctype)))
    
    if(balanced){
      print("Balanced cancer type representation...")
      sample_to <- quantile(ctype_cnts, quantile_cutoff) ## Number of samples per cancer type
      
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
    print("Separating cancers by cancer type...")
    segd_ctype <- split(segd, segd$ctype)
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
analysis <- 'pancan'
balanced=TRUE

PDIR <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/TCGA_signatures'
annof <- file.path(PDIR, "input", 'merged_sample_quality_annotations.tsv')
anno <- read.table(annof, sep="\t", stringsAsFactors = FALSE,  quote = '', 
                   header=TRUE, check.names = FALSE, fill = FALSE, 
                   comment.char = '', strip.white = TRUE)

segf <- file.path(PDIR, "input", "TCGA_mastercalls.abs_segtabs.fixed.txt")
segd <- read.table(segf, sep="\t", header=TRUE, stringsAsFactors = FALSE)
sample_anno <- data.frame("sample"=unique(segd$Sample))
sample_anno$ctype <- sapply(sample_anno$sample, function(s){
  anno[grep(s, anno$aliquot_barcode)[1],]$'cancer type'
})
s_id <- setNames(sample_anno$ctype, sample_anno$sample)

#### 1) Address class imbalance - balance cancer types  ####
segd$ctype <- s_id[segd$Sample]
segd_ctype <- preprocessSegs(segd, analysis, balanced, 0.50, min_n=TRUE)
# length(unique(segd_ctype[[1]]$Sample)) ## 4554


#### 2) Parse seg data into Sigminer:CN counts ####
cn_ctype <- lapply(segd_ctype, function(segd){
  print(unique(segd$ctype))
  cn <- sigminer::read_copynumber(segd, samp_col = "Sample", max_copynumber = 20L,
                                  seg_cols = c("Chromosome", "Start", "End", "Modal_Total_CN"),
                                  genome_build = "hg19", complement = FALSE, verbose = TRUE)
  return(cn)
})
save(cn_ctype, file=file.path(PDIR, "output", paste0("cn_ctype_", analysis, 
                                                     "-", balanced, ".rda")))

#### 3) Convert CN counts into the feature-set ####
load(file.path(PDIR, "output", paste0("cn_ctype_", analysis, "-", balanced, ".rda"))) #cn_ctype
cn_tally_ctype <- lapply(cn_ctype, function(cn){
  cn_tally_W <- sig_tally(cn, method = "W")
  return(cn_tally_W)
})
save(cn_tally_ctype, file=file.path(PDIR, "output", paste0("cn_tally_ctype_",
                                    analysis, "-", balanced, ".rda")))

#### 4) Deconvolute Feature-set matrix using NMF ####
library(NMF)
## ---- Bayesian NMF "Auto" version
load(file.path(PDIR, "output", paste0("cn_tally_ctype_",
                                      analysis, "-", balanced, ".rda"))) #cn_tally_ctype
sig_auto_ctype <- lapply(cn_tally_ctype, function(cn_tally_W){
  sig_auto_extract(cn_tally_W$nmf_matrix, cores=1, nrun = 5, K0=30,
                   method='L1W.L2H', strategy='stable', destdir = file.path(PDIR, "tmp"))
})
save(sig_auto_ctype, file=file.path(PDIR, "output", paste0("cn_autosig_ctype_", 
                                    analysis, "-", balanced, ".rda")))

## ---- Manual method to select optimal K
load(file.path(PDIR, "output", paste0("cn_tally_ctype_",
                                      analysis, "-", balanced, ".rda"))) #cn_tally_ctype
sig_w_ctype <- lapply(cn_tally_ctype, function(cn_tally_W){
  sig_estimate(cn_tally_W$nmf_matrix, range=2:30, cores=2, keep_nmfObj=TRUE, use_random=TRUE,
               nrun = 30, verbose=TRUE, save_plots=TRUE, pConstant = 1e-13,
               plot_basename=file.path(PDIR, "output", "nmf"))
})
save(sig_w_ctype, file=file.path(PDIR, "output", paste0("cn_nmf_ctype_", 
                                                        analysis, "-", balanced, ".rda")))

#### 5) Visualize CN signatures ####
sig_method='auto' # 'auto' or 'ctype'
sig_type <- switch(sig_method,
                   auto='cn_autosig_ctype_',
                   ctype='cn_nmf_ctype_')
sig_obj_id <- load(file.path(PDIR, "output", paste0(sig_type, analysis, "-", balanced, ".rda")))
load(file.path(PDIR, "output", paste0(sig_type, analysis, "-", balanced, ".rda")))
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

pdf(file.path(PDIR, "output", paste0(sig_type, "TCGA-exposure.pdf")), width = 60, height = 8)
show_sig_exposure(cn_sig[[1]], palette=sample(col_vector, cn_sig[[1]]$K), 
                  sig_names=gsub("Sig", "", rownames(cn_sig[[1]]$Exposure)),
                  rm_panel_border=TRUE, rm_space=TRUE,
                  groups = anno_grps,
                  grp_order = sort(unique(anno_grps)))
dev.off()

pdf(file.path(PDIR, "output", paste0(sig_type, "TCGA-profile.pdf")), width = 10, height = 10)
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
sig_obj_id <- load(file.path(PDIR, "output", paste0(sig_type, analysis, "-", balanced, ".rda")))
load(file.path(PDIR, "output", paste0(sig_type, analysis, "-", balanced, ".rda")))
cn_sig <- switch(sig_method,
                 auto=sig_auto_ctype,
                 ctype=sig_w_ctype)

segd2_wN <- formatSeg(segd = segd, analysis = 'TCGA', PDIR = dirname(PDIR))
cn_wN <- sigminer::read_copynumber(segd2_wN, samp_col = "Sample", max_copynumber = 20L,
                                seg_cols = c("chrom", "loc.start", "loc.end", "Modal_Total_CN"),
                                genome_build = "hg19", complement = FALSE, verbose = TRUE)
save(cn_wN, file=file.path(PDIR, "output", 'tcga-1000g_CN.rda'))
cn_tally_wN <- sig_tally(cn_wN, method = "W")



#### 7) Fit CCL samples to CN-signatures ####




load(file.path(PDIR, "output", paste0("cn_ctype_", analysis, "-", balanced, ".rda")))
load(file.path(PDIR, "output", paste0("cn_tally_ctype_",
                                      analysis, "-", balanced, ".rda")))
load(file.path(PDIR, "output", paste0("cn_autosig_ctype_", 
                                      analysis, "-", balanced, ".rda")))
load(file.path(PDIR, "output", paste0("cn_nmf_ctype_", 
                                      analysis, "-", balanced, ".rda")))




glm_model <- MASS::rlm(cophenetic ~ poly(rank, degree = 3), 
                       data=sig_w_ctype$PANCAN$survey)
nmf_loess <- loess(sig_w_ctype$PANCAN$survey$cophenetic ~ sig_w_ctype$PANCAN$survey$rank, span=1.5)
hweights <- data.frame(rank = sig_w_ctype$PANCAN$survey$rank, resid = glm_model$resid, weight = glm_model$w)
hweights[order(glm_model$w), ]


plot(cophenetic ~ rank, data=sig_w_ctype$PANCAN$survey)
lines(cophenetic ~ rank, data=sig_w_ctype$PANCAN$survey)
lines(y=predict(glm_model, list(rank=seq(2, 15, by=0.1))), 
      x=seq(2, 15, by=0.1), col="green")
lines(y=predict(nmf_loess, seq(2, 15, by=0.1)), x=seq(2, 15, by=0.1), col="blue")


barplot(abs(setNames(diff(predict(nmf_loess, 2:15)), 2:14)))
sig_w <- sig_extract(cn_tally_ctype[[1]]$nmf_matrix, n_sig = 2, pConstant = 1e-13)

sig_w_ctype <- lapply(cn_tally_ctype, function(cn_tally_W){
  sig_w <- sig_extract(cn_tally_W$nmf_matrix, n_sig = 2, pConstant = 1e-13)
  sig_w
})


options(sigminer.copynumber.max = 10)




cn_tally_W <- sig_tally(cn, method = "W")



segl <- split(segd[,c(2:6,9)], segd$Sample)
segl <- segl[-which(sapply(segl, nrow)==1)]

library(dplyr)
library(flexmix)
library(NMF)
library(YAPSA)
this_path <- "~/git/cnsignatures"
source(paste(this_path,"helper_functions.R",sep="/"))
source(paste(this_path,"main_functions.R",sep="/"))

# setwd(file.path(this_path, "manuscript_Rmarkdown"))
# eset.dir <- '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/seg_files/esets'
# handle <- 'gCSI'
# eset <- file.path(eset.dir, paste0(handle, "_bin_eset.Rdata"))
# load(eset) # cl.eset


PDIR <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/TCGA_signatures'
segf <- file.path(PDIR, "input", "TCGA_mastercalls.abs_segtabs.fixed.txt")
segd <- read.table(segf, sep="\t", header=TRUE, stringsAsFactors = FALSE)
colnames(segd)[2:4] <- c('chrom', 'loc.start', 'loc.end')

chr.size.dat <- getChrLength()
seqlevelsStyle(chr.size.dat) <- 'NCBI'
scale <- 1e5




order=7



#read.table(file.path(tcga_dir, "input", "TCGA_mastercalls.abs_segtabs.fixed.txt"))

tcga_dir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/TCGA_signatures'
segl <- split(segd[,c(2:6,9)], segd$Sample)
segl <- segl[-which(sapply(segl, nrow)==1)]
CN_features <- extractCopynumberFeatures(segl, pset.assay='Modal_Total_CN', cores = 8)

############################################
#### Decompose into Mixture of Gaussian ####
seed=77777
min_prior=0.001
model_selection="BIC"
nrep=1
niter=1000

dat<-as.numeric(CN_features[["segsize"]][,2])
segsize_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                         min_prior=min_prior,niter=niter,nrep=nrep,min_comp=10,max_comp=10)
save(segsize_mm, file=file.path(PDIR, "output", "newfit", "mm_segsize.rda"))

dat<-as.numeric(CN_features[["bp10MB"]][,2])
bp10MB_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                        min_prior=min_prior,niter=niter,nrep=nrep,min_comp=3,max_comp=3)
save(bp10MB_mm, file=file.path(PDIR, "output", "newfit", "mm_bp10MB.rda"))

dat<-as.numeric(CN_features[["osCN"]][,2])
osCN_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                      min_prior=min_prior,niter=niter,nrep=nrep,min_comp=3,max_comp=3)
save(osCN_mm, file=file.path(PDIR, "output", "newfit", "mm_osCN.rda"))

dat<-as.numeric(CN_features[["bpchrarm"]][,2])
bpchrarm_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                          min_prior=min_prior,niter=niter,nrep=3,min_comp=5,max_comp=5)
save(bpchrarm_mm, file=file.path(PDIR, "output", "newfit", "mm_bpchrarm.rda"))

dat<-as.numeric(CN_features[["changepoint"]][,2])
changepoint_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                             min_prior=min_prior,niter=niter,nrep=nrep,min_comp=7,max_comp=7)
save(changepoint_mm, file=file.path(PDIR, "output", "newfit", "mm_changepoint.rda"))

dat<-as.numeric(CN_features[["copynumber"]][,2])
copynumber_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                            nrep=nrep,min_comp=8,max_comp=8,min_prior=0.005,niter=2000)
save(copynumber_mm, file=file.path(PDIR, "output", "newfit", "mm_copynumber.rda"))

CN_components<-list(segsize=segsize_mm,
                    bp10MB=bp10MB_mm,
                    osCN=osCN_mm,
                    changepoint=changepoint_mm,
                    copynumber=copynumber_mm,
                    bpchrarm=bpchrarm_mm)
save(CN_components, file=file.path(PDIR, "output", "newfit", "components.rda"))



comp.mat <- generateSampleByComponentMatrix(CN_features, CN_components,
                                            cores=1, subcores=num_cores)
save(comp.mat, file=file.path(PDIR, "output", "newfit", "sample_comp_mat.rda"))

pdf(file.path(PDIR, "output", "newfit", "heatmap_tcga.pdf"))
NMF::aheatmap(comp.mat,
              fontsize = 7, Rowv=FALSE, Colv=FALSE,
              legend = T, breaks=c(seq(0,199,2),500), 
              main="Component x Sample matrix")
dev.off()


nsig<-7
nmfalg<-"brunet"
num_cores <- 4
seed=77777
sigs<-NMF::nmf(t(comp.mat),
               nsig,seed=seed,nrun=1000,
               method=nmfalg,.opt = paste0("p",num_cores))

coefmap(sigs,Colv="consensus",tracks=c("basis:"), main=paste0(handle, ": Patient x Signature matrix"))
basismap(sigs,Rowv=NA,main=paste0(handle, ": Signature x Component matrix"))



save(comp.mat, CN_components, CN_features, sigs, 
     file=file.path("~", paste0("cnsignatures_", handle, ".RData")))
britroc_sample_component_matrix<-generateSampleByComponentMatrix(CN_features, CN_components,
                                                                 cores=1, subcores=num_cores)



#### Aggregate and compare multiple samples ####
require(dplyr)
fls <- list.files("~", pattern=("cnsignatures"))
handles <- gsub("cnsignatures_", "", fls) %>% gsub(".RData", "", .)
handles <- handles[handles != 'gcsi']

load(paste0("~/cnsignatures_", 'gCSI', ".RData"))
ref.components <- CN_components

all.sigs <- lapply(handles, function(h){
  # "comp.mat": sample by component matrix
  # "CN_features": Raw CN feature data
  # "CN_components": Mixed Gaussian applied to CN features to get components
  # "sigs": NMF decomposed Sample x Component matrix into Signatures
  nsig<-7
  nmfalg<-"brunet"
  num_cores <- 4
  seed=77777
  
  load(paste0("~/cnsignatures_", h, ".RData"))
  h.component_matrix <- generateSampleByComponentMatrix(CN_features,ref.components)
  sigs <- NMF::nmf(t(h.component_matrix), nsig, seed=seed, 
                   nrun=1000,method=nmfalg,.opt = "p16")
  comp_sig_mat<-basis(sigs)     # Signature x Component
  patient_sig_mat<-coef(sigs)   # Patient x Signature
  
  list("comp.mat"=h.component_matrix,
       "sigs"=sigs,
       "basis"=comp_sig_mat,
       "coef"=patient_sig_mat)
})

all.sigs <- lapply(handles, function(h){
  # "comp.mat": sample by component matrix
  # "CN_features": Raw CN feature data
  # "CN_components": Mixed Gaussian applied to CN features to get components
  # "sigs": NMF decomposed Sample x Component matrix into Signatures
  load(paste0("~/cnsignatures_", h, ".RData"))
  
  #britroc feat_sig matrix
  comp_sig_mat<-basis(sigs)     # Signature x Component
  patient_sig_mat<-coef(sigs)   # Patient x Signature
  
  list("basis"=comp_sig_mat,
       "coef"=patient_sig_mat)
})
names(all.sigs) <- handles
all.sigs[['gCSI']][['basis']]

sig1 <- sapply(all.sigs, function(i) i[['coef']][,'SW948'])
cor(sig1, method = 'spearman')


#### Compare signatures across gCSI, CCLE, and GDSC ####
reord_components<-c(11:13,24:31,17:23,32:36,14:16,1:10)

reord_britroc<-as.integer(c(2,6,5,4,7,3,1))
names(reord_britroc)<-paste0("s",1:7)
feat_sig_mat<-feat_sig_mat[,reord_britroc]
colnames(feat_sig_mat)<-paste0("s",1:nsig)
sig_feat_mat<-t(feat_sig_mat)

#pcawg feat_sig matrix
feat_sig_mat_pcawg<-basis(pcawg_sigs)[,]
sig_feat_mat_pcawg<-t(feat_sig_mat_pcawg)
colnames(feat_sig_mat_pcawg)<-paste0("s",1:nsig)

#tcga feat_sig matrix
feat_sig_mat_tcga<-basis(tcga_sigs)[,]
sig_feat_mat_tcga<-t(feat_sig_mat_tcga)
colnames(feat_sig_mat_tcga)<-paste0("s",1:nsig)


#determine matching signatures and their correlation
reord_pcawg<-apply(feat_sig_mat,2,function(x){which.max(apply(feat_sig_mat_pcawg,2,cor,x,method="s"))})
sigcor_pcawg<-apply(feat_sig_mat,2,function(x){max(apply(feat_sig_mat_pcawg,2,cor,x,method="s"))})

reord_tcga<-apply(feat_sig_mat,2,function(x){which.max(apply(feat_sig_mat_tcga,2,cor,x,method="s"))})
sigcor_tcga<-apply(feat_sig_mat,2,function(x){max(apply(feat_sig_mat_tcga,2,cor,x,method="s"))})

#plot the feat_sig matrices side by side
par(mfrow=c(1,3))
basismap(sigs,Rowv=NA,Colv=reord_britroc,main="BritROC",tracks=NA)
basismap(pcawg_sigs,Rowv=NA,Colv=reord_pcawg,main="PCAWG",tracks=NA)
basismap(tcga_sigs,Rowv=NA,Colv=reord_tcga,main="TCGA",tracks=NA)
