# R code for cnvML
## Overview

| Code | Section | Figure | Description |
|------|------|------|------|
| <a href="#plotsfcr">plotSFC.R</a> | Preprocessing |  | Plots Hilbert/Sweep ASCN and TCN SFC |
|  <a href="#arrowssfcr">arrowsSFC.R</a> | Supplementary Figure | Fig S1 | Mapping of SFC |
| <a href="#comparesfcr">compareSFC.R</a> | Supplementary Figure | Fig S2 | Sweep vs Hilbert SFC |

## Pre-processing
* plotSFC.R is the main function to plot the Hilbert spacefilling curves for TCGA and the cancer cell line seg files
```
Options:
	-c CNTYPE, --cntype=CNTYPE
		Copy number to plot: TCN or ASCN [Default=ASCN]
	-p PDIR, --pdir=PDIR
		Path to seg files [dir]/[dataset]/input/[seg_file] [Default=/mnt/work1/users/pughlab/projects/cancer_cell_lines]
	-o ORDER, --order=ORDER
		Order for space filling curve (SFC) [Default=8]
	-m MAXCN, --maxcn=MAXCN
		Max CN to plot to [Default=5]
	-s SFC, --sfc=SFC
		Space-filling curve to use, 'sweep' or 'hilbert' [Default=hilbert]
	-d DATASET, --dataset=DATASET
		Dataset to use, 'ccl_aggregate' or 'TCGA' [Default=ccl_aggregate]
	-v, --verbose
		Print extra output [Default=FALSE]
	-h, --help
		Show this help message and exit
```

## Main Figures

## Supplementary Figures
### arrowsSFC.R 
This piece of code is used to generate **Supplementary Figure 2: Mapping of Hilbert and Sweep space filling curves**.

Using a 4th and 6th orders/level, it repurposes the `HilbertCurve::HilbertCurve()` function to generate a directional graph of either a HilbertCurve or a SweepCurve. This adapted function allows for the input of a different positional space (@POS), which allows for a Sweep Curve to be used. If no `sfc_pos` positional space is given, it defaults to the innate Hilbert Curve space.

### compareSFC.R 
This piece of code is used to generate **Supplementary Figure 1: Comparison between Hilbert and Sweep space filling curves**.

Using several different orders/levels (4, 6, 8), it will generate a Hilbert and Sweep space-filling curve using the standard hg19 genome, colored by chromosomes. It will then run a Euclidean distance between adjacent regions in the 2D space for both SFC.

