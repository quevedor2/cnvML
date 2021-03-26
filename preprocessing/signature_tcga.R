library(dplyr)
library(flexmix)
library(NMF)
library(YAPSA)
this_path <- "~/git/cnsignatures"
source(paste(this_path,"helper_functions.R",sep="/"))
source(paste(this_path,"main_functions.R",sep="/"))

###################
#### Functions ####
getChrLength <- function(){
  require(BSgenome.Hsapiens.UCSC.hg19)
  chr.lengths = seqlengths(Hsapiens)[1:24]
  chr.len.gr <- makeGRangesFromDataFrame(data.frame("chrom"=names(chr.lengths),
                                                    "loc.start"=rep(1, length(chr.lengths)),
                                                    "loc.end"=chr.lengths))
  chr.len.gr$cum.end <- cumsum(as.numeric(end(chr.len.gr)))
  chr.len.gr$cum.start <- chr.len.gr$cum.end - (end(chr.len.gr) -1)
  chr.len.gr$cum.mid <- (chr.len.gr$cum.start + ((chr.len.gr$cum.end - chr.len.gr$cum.start)/2))
  return(chr.len.gr)
}

addCumPos <- function(dat, ref, dat.type){
  m.row.idx <- match(as.character(dat$chrom), as.character(seqnames(ref)))
  if(dat.type=='data'){
    dat$cpos <- ref[m.row.idx,]$cum.start +  dat$pos - 1
    dat$chr.stat
  } else if(dat.type == 'seg'){
    dat$cloc.start <- ref[m.row.idx,]$cum.start +  dat$loc.start - 1
    dat$cloc.end <- ref[m.row.idx,]$cum.start +  dat$loc.end - 1
  }
  dat$chr.stat <- (m.row.idx %% 2) + 1
  return(dat)
}

##############
#### Main ####

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
scale <- 1e4
order=8



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
