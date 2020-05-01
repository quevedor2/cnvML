library(dplyr)
library(flexmix)
library(NMF)
library(YAPSA)
this_path <- "~/git/cnsignatures"
source(paste(this_path,"helper_functions.R",sep="/"))
source(paste(this_path,"main_functions.R",sep="/"))

setwd(file.path(this_path, "manuscript_Rmarkdown"))


eset.dir <- '/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/seg_files/esets'
handle <- 'gCSI'
eset <- file.path(eset.dir, paste0(handle, "_bin_eset.Rdata"))
load(eset) # cl.eset

CN_features <- extractCopynumberFeatures(cl.eset, cores = 8)


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

dat<-as.numeric(CN_features[["bp10MB"]][,2])
bp10MB_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                        min_prior=min_prior,niter=niter,nrep=nrep,min_comp=3,max_comp=3)

dat<-as.numeric(CN_features[["osCN"]][,2])
osCN_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                      min_prior=min_prior,niter=niter,nrep=nrep,min_comp=3,max_comp=3)

dat<-as.numeric(CN_features[["bpchrarm"]][,2])
bpchrarm_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                          min_prior=min_prior,niter=niter,nrep=3,min_comp=5,max_comp=5)

dat<-as.numeric(CN_features[["changepoint"]][,2])
changepoint_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                             min_prior=min_prior,niter=niter,nrep=nrep,min_comp=7,max_comp=7)

dat<-as.numeric(CN_features[["copynumber"]][,2])
copynumber_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                            nrep=nrep,min_comp=8,max_comp=8,min_prior=0.005,niter=2000)

CN_components<-list(segsize=segsize_mm,
                    bp10MB=bp10MB_mm,
                    osCN=osCN_mm,
                    changepoint=changepoint_mm,
                    copynumber=copynumber_mm,
                    bpchrarm=bpchrarm_mm)




comp.mat<-generateSampleByComponentMatrix(CN_features, CN_components,
                                          cores=1, subcores=num_cores)

NMF::aheatmap(comp.mat,
              fontsize = 7, Rowv=FALSE, Colv=FALSE,
              legend = T, breaks=c(seq(0,199,2),500), 
              main="Component x Sample matrix")


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
