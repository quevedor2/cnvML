library(HilbertCurve)
library(HilbertVis)
library(circlize)
library(DescTools)
library(GenomicRanges)
library(dplyr)
library(optparse)
library(RcnvML)   # devtools::install_github("quevedor2/cnvML/RcnvML", ref='dev')

####################
#### Parameters ####
option_list <- list( 
  make_option(c("-c", "--cntype"), type="character", default='ASCN',
              help="Copy number to plot: TCN or ASCN [default]"),
  make_option(c("-d", "--dir"), type="character", 
              default='/mnt/work1/users/pughlab/projects/cancer_cell_lines',
              help="Path to seg files [dir]/[dataset]/input/[seg_file] [default]"),
  make_option(c("-o", "--order"), type="integer", default=8,
              help="Order for space filling curve (SFC) [default]"),
  make_option(c("-m", "--maxcn"), type="integer", default=5,
              help="Max CN to plot to [default]"),
  make_option(c("-s", "--sfc"), type="character", default='hilbert',
              help="Space-filling curve to use, 'sweep' or 'hilbert' [default]"),
  make_option(c("-d", "--dataset"), type="character", default='ccl_aggregate',
              help="Dataset to use, 'ccl_aggregate' or 'TCGA' [default]"),
  make_option(c("-v", "--verbose"), action="store_false", default=FALSE,
              help="Print extra output [default]")
)
opt <- parse_args(OptionParser(option_list=option_list))


# Set up parameters
PDIR <- opt$dir                   # Path to seg files [file.path(PDIR, analysis, input, seg_i)]
analysis <- opt$dataset           # dataset to analyze (TCGA or ccl_aggregate)
chr.size.dat <- getChrLength()    # hg19 chromosome sizes
seqlevelsStyle(chr.size.dat) <- 'NCBI'  # Chr labelling style
scale <- 1e4                      # Divides genomic loci by this value
order=opt$order                   # Order/level for space-filling curve (SFC)
a_col=c("black", "red")           # Range of color for the A-allele (or TCN)
b_col=c("black", "blue")          # Range of color for the B-allele
breaks <- seq(0, opt$maxcn, by=0.1) # Where to add the breaks for Chr Integer CN (Increase for TCN)
add_1000g <- TRUE                 # Adds 1000G data to the TCGA seg file
sfc <- opt$sfc                    # whether to use a hilbert SFC or sweep
verbose <- FALSE                  # Verbose, TRUE or FALSE
cntype <- opt$cntype                # Plot ASCN or TCN as colours

# Asserts
assert_that(length(cntype)==1, any(c('TCN', 'ASCN') %in% cntype), 
            msg='--cntype must be either "TCN" or"ASCN"')
assert_that(length(sfc)==1, any(c('sweep', 'hilbert') %in% sfc), 
            msg='--sfc must be either "sweep" or "hilbert"')
assert_that(length(analysis)==1, any(c('ccl_aggregate', 'TCGA') %in% analysis), 
            msg='--dataset must be either "ccl_aggregate" or "TCGA"')
assert_that(is.integer(order), length(order)==1, order>3,
            msg='--order must be a single integer greater than 3')
assert_that(is.integer(order), length(order)==1, order>3,
            msg='--order must be a single integer greater than 3')
assert_that(is.integer(opt$maxcn), length(opt$maxcn)==1, opt$maxcn>2,
            msg='--maxcn must be a single integer greater than 2')

###########################
#### Read in Seg files #### 
annotate <- FALSE
seg_1000G <- file.path(PDIR, '1000G', 'eacon', 'symlinks', 'seg', '1000g.seg')
seg_files <- switch(analysis,
                    'TCGA'="TCGA_mastercalls.abs_segtabs.fixed.txt",
                    'ccl_aggregate'=c('CCLE_cna_hg19.seg', 
                                      'GDSC_cna_hg19.seg', 
                                      'gCSI_cna_hg19.seg'),
                    stop("Analysis must be either 'TCGA' or 'ccl_aggregate'"))

####################
#### Seg to SFC #### 
for(seg_i in seg_files){
  segf <- file.path(PDIR, analysis, "input", seg_i)
  segd <- read.table(segf, sep="\t", header=TRUE, stringsAsFactors = FALSE)
  
  ##########################
  #### Format seg files ####
  if(analysis=='TCGA'){
    colnames(segd)[2:4] <- c('chrom', 'loc.start', 'loc.end')
    
    if(add_1000g){
      # add 1000Genome to TCGA seg
      normal_segd <- read.table(seg_1000G, sep="\t", header=TRUE, stringsAsFactors = FALSE)
      colnames(normal_segd) <- c('Sample', 'chrom', 'Chrom', 'loc.start', 'loc.end', 
                                 'Length', 'Modal_Total_CN', 'Modal_HSCN_1', 'Modal_HSCN_2')
      segd2 <- plyr::rbind.fill(segd, normal_segd)
      #segd <- segd2[grep("^HG", segd2$Sample),]
      segd <- segd2
    }
  } else {
    colnames(segd)[c(1:3,6,11,7,8)] <- c('chrom', 'loc.start', 'loc.end',
                                  'Sample', 'Modal_Total_CN', 
                                  'Modal_HSCN_1', 'Modal_HSCN_2')
    if(any(segd$chrom %in% c('chrX', 'chrY'))){
      segd <- segd[-which(segd$chrom %in% c('chrX', 'chrY')),]
    }
  }
  
  #####################################
  #### Create Space Filling Curves ####
  segl <- split(segd, segd$Sample)
  hilberts <- lapply(segl, function(seg){
    # Print Sample ID
    sample_id <- unique(seg$Sample) 
    print(paste0(seg_i, "-Sample: ", sample_id, " (", 
                 grep(sample_id, names(segl)), "/", length(segl), ")"))
    
    if(!all(is.na(seg$chrom))){
      ## Set up the chromosome IDs
      seg$chrom <- gsub("23", "X", seg$chrom) %>%
        gsub("24", "Y", .) %>%
        gsub("^chr", "", ., ignore.case = TRUE)
      
      ## Setup the chromosome coloring
      a1 = circlize::colorRamp2(breaks = breaks,
                                colors = colorRampPalette(a_col)(length(breaks)))
      a2 = circlize::colorRamp2(breaks = breaks,
                                colors = colorRampPalette(b_col)(length(breaks)))
      
      ## Plot Allele specific copy-number
      png(file.path(PDIR, analysis, "output", "hilbert", paste0(sample_id, ".png")), 
          width=300, height=300)
      hc_obj <- setupRefHcMatrix(order = order, scale=scale)
      if(sfc == 'sweep'){
        if(verbose) print("Using a Sweep SFC")
        bins <- mapSPC(spc='sweep', hc_ord=hc_obj$ord)
        hc_obj$hc@POS <- bins[-nrow(bins),c('x1', 'y1', 'x2', 'y2')]
      }
      
      
   
      # Adds hte colour track for the A- and B-allele CN
      # Ensure that Red/Blue channels range from 0-255. MixColor halves both 
      # values when mixing
      chr.data <- addCumPos(seg, chr.size.dat, dat.type='seg')
      ir <- IRanges(chr.data$cloc.start/scale, chr.data$cloc.end/scale)
      if(cntype='TCN'){
        # Expecting column: Modal_Total_CN
        assert_that(all(c('Modal_Total_CN') %in% colnames(chr.data)),
                    msg='Could not find both Modal_Total_CN in chr.data')
        cols <- a1(chr.data$Modal_Total_CN)
      } else (cntype='ASCN'){
        # Expecting column: Modal_HSCN_1, Modal_HSCN_2
        assert_that(all(c('Modal_HSCN_1', 'Modal_HSCN_2') %in% colnames(chr.data)),
                    msg='Could not find both Modal_HSCN_1/2 in chr.data')
        cols <- DescTools::MixColor(col1 = a1(chr.data$Modal_HSCN_1), 
                                    col2 = a2(chr.data$Modal_HSCN_2))
        cols <- apply(col2rgb(cols),2,function(x) rgb(x[1]*2, x[2]*2, x[3]*2, 
                                                      maxColorValue = 255))
      }
      hc_layer(hc_obj$hc, ir, mean_mode = "absolute", col = cols)
      
      if(annotate){
        # Adds lines separating chromosomes and text for chromosomes
        hc_polygon(hc_obj$hc, gp=gpar(col='white'),
                   x1=chr.size.dat$cum.start/scale, x2=chr.size.dat$cum.end/scale)
        hc_text(hc_obj$hc, x1=chr.size.dat$cum.start/scale, x2=chr.size.dat$cum.end/scale, 
                labels = c(1:22, "X", "Y"), gp = gpar(fontsize = 10, col='white'), 
                centered_by = "polygon")
      }
      dev.off()
      
    }
  })
}






