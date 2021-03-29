# R code for cnvML
## Overview

| Code | Figure | Description |
|------|------|------|
| compareSFC.R | Fig S1 | Sweep vs Hilbert SFC |

## Pre-processing
* plotSFC.R is the main function to plot the Hilbert spacefilling curves for TCGA and the cancer cell line seg files
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

## Main Figures

## Supplementary Figures
### compareSFC.R 
This piece of code is used to generate **Supplementary Figure 1: Comparison between Hilbert and Sweep space filling curves**.

Using several different orders/levels (4, 6, 8), it will generate a Hilbert and Sweep space-filling curve using the standard hg19 genome, colored by chromosomes. It will then run a Euclidean distance between adjacent regions in the 2D space for both SFC.

