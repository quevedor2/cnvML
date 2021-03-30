# cnvML
This repo is meant to document the code used for: `Boosting aneuploidy based predictions through models that utilize an allele-specific encoding of copy number`. This paper focuses on cancer type and drug predictions in TCGA and GDSC/CCLE/gCSI using ML models

# Set-up
The `RcnvML` package needs to be installed as it contains many functions that supports the analysis:
```
devtools::install_github("quevedor2/cnvML/RcnvML")
```

# Analysis
## Generating CN Features
### Input files
Seg files are fed as the input for all gene/bin/SFC mappings. The code is expecting the following seg files for:
* `ccl_aggregate`
  * CCLE_cna_hg19.seg
  * gCSI_cna_hg19.seg
  * GDSC_cna_hg19.seg
* `TCGA`
  * TCGA_mastercalls.abs_segtabs.fixed.txt
  * 1000g.seg
  
### Mapping Seg files to Bins

### Mapping Seg files to Genes

### Mapping Seg files to Space-Filling Curves
The R script [plotSFC.R](https://github.com/quevedor2/cnvML/blob/master/R/plotSFC.R) can be used to create the SFC png images, but expects an input structure as follows:
```
[dataset]/
└── input
    ├── [seg_file1.seg]
    ├── [seg_file2.seg]
    └── [seg_file3.seg]
```

The following commands were used to create the SFC images
```sh
# ASCN plots using a Sweep SFC
Rscript plotSFC.R --cntype 'ASCN' --order 8 --maxcn 5 --sfc sweep --dataset ccl_aggregate
Rscript plotSFC.R --cntype 'ASCN' --order 8 --maxcn 5 --sfc sweep --dataset TCGA

# ASCN plots using a Hilbert SFC
Rscript plotSFC.R --cntype 'ASCN' --order 8 --maxcn 5 --sfc hilbert --dataset ccl_aggregate
Rscript plotSFC.R --cntype 'ASCN' --order 8 --maxcn 5 --sfc hilbert --dataset TCGA

# TCN plots using a Hilbert SFC
Rscript plotSFC.R --cntype 'TCN' --order 8 --maxcn 8 --sfc hilbert --dataset ccl_aggregate
Rscript plotSFC.R --cntype 'TCN' --order 8 --maxcn 8 --sfc hilbert --dataset TCGA
```

## Mapping features to cancer type
The following commands were used to create the SFC images