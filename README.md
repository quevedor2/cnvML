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
The R script [assignOncocode.R](https://github.com/quevedor2/cnvML/blob/master/R/assignOncocode.R) is used to redirect the otput of the plotSFC png images into a directory structure of oncocodes. This is meant to be fed into pickle for downstream neural network building. The mapping is obtained from the data() files from the RcnvML package:
```
# ASCN plots using a Sweep SFC
Rscript assignOncocode.R --cntype 'ASCN' --sfc sweep --dataset ccl_aggregate
Rscript assignOncocode.R --cntype 'ASCN' --sfc sweep --dataset TCGA

# ASCN plots using a Hilbert SFC
Rscript assignOncocode.R --cntype 'ASCN' --sfc hilbert --dataset ccl_aggregate
Rscript assignOncocode.R --cntype 'ASCN' --sfc hilbert --dataset TCGA

# TCN plots using a Hilbert SFC
Rscript assignOncocode.R --cntype 'TCN' --sfc hilbert --dataset ccl_aggregate
Rscript assignOncocode.R --cntype 'TCN' --sfc hilbert --dataset TCGA
```

## Preprocessing input for tensorflow
The python script [pickle_hilbert.py](https://github.com/quevedor2/cnvML/blob/master/common/pickle_hilbert.py) is used to create pickle structures of the SFC png images that have been sectioned off into their oncocodes. It will also reshape the images to have a consistent img size (300), and randomly shuffle the data (while maintaining the labels) to allow for training/test groups at a later stage.
```
# Pickle the TCGA data
python pickle_hilbert.py -a TCGA -s sweep -c ASCN -d /cluster/projects/pughlab/projects/cancer_cell_lines/
python	pickle_hilbert.py -a TCGA -s hilbert -c ASCN -d /cluster/projects/pughlab/projects/cancer_cell_lines/
python	pickle_hilbert.py -a TCGA -s hilbert -c TCN -d /cluster/projects/pughlab/projects/cancer_cell_lines/

# Pickle the CCL data
python	pickle_hilbert.py -a CCL -s sweep -c ASCN -d /cluster/projects/pughlab/projects/cancer_cell_lines/
python	pickle_hilbert.py -a CCL -s hilbert -c ASCN -d /cluster/projects/pughlab/projects/cancer_cell_lines/
python	pickle_hilbert.py -a CCL -s hilbert -c TCN -d /cluster/projects/pughlab/projects/cancer_cell_lines/
```
