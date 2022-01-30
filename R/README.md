# R code for cnvML
## Overview

| Code | Section | Figure | Description |
|------|------|------|------|
| <a href="#segtobins">segToBins.R</a> | Preprocessing |  | Converts the ASCN segs to genomic bins that match SFC segmentation |
| <a href="#segtogenes">segToGenes.R</a> | Preprocessing |  | Annotates all the genes in the genome with their corresponding ASCN value |
| <a href="#plotsfcr">plotSFC.R</a> | Preprocessing |  | Plots Hilbert/Sweep ASCN and TCN SFC |
| <a href="#assignoncocoder">assignOncocode.R</a> | Preprocessing |  | Moves SFC pngs to oncocode directory structure |
| <a href="#arrowssfcr">arrowsSFC.R</a> | Supplementary Figure | Fig S1 | Mapping of SFC |
| <a href="#comparesfcr">fillsSFC.R</a> | Supplementary Figure | Fig S2 | 2D to 1D SFCs |
| <a href="#comparesfcr">compareSFC.R</a> | Supplementary Figure | Fig S3 | Sweep vs Hilbert SFC |

## Pre-processing
### segToBins.R
* segsToBins.R is the main function to convert the ASCN seg data and bin them into genomic segments that are equivalent to the SFC segmentation. Three Sample-by-Bin matrices are generated from this script, one for Total_CN, HSCN_1 and HSCN_2, as well as a reference bins file that creates the mapping between genomic bin ID and coordinates.
```
Options:
	-p PDIR, --pdir=PDIR
		Path to seg files [dir]/[dataset]/input/[seg_file] [/mnt/work1/users/pughlab/projects/cancer_cell_lines]

	-d DATASET, --dataset=DATASET
		Dataset to use, 'ccl_aggregate' or 'TCGA' [ccl_aggregate]

	-s SEGFILE, --segfile=SEGFILE
		If --dataset is set to custom, indicate the name of the seg_file [NULL]

	-m METAFILE, --metafile=METAFILE
		If --dataset is set to custom, indicate the meta RDS file path [NULL]

	-h, --help
		Show this help message and exit
```

### segToGenes.R
* segsToGenes.R is the main function to convert the ASCN seg data into gene-level ASCN representation. Three Sample-by-Gene matrices are generated from this script, one for Total_CN, HSCN_1 and HSCN_2.
```
Options:
	-p PDIR, --pdir=PDIR
		Path to seg files [dir]/[dataset]/input/[seg_file] [/mnt/work1/users/pughlab/projects/cancer_cell_lines]

	-d DATASET, --dataset=DATASET
		Dataset to use, 'ccl_aggregate' or 'TCGA' [ccl_aggregate]

	-s SEGFILE, --segfile=SEGFILE
		If --dataset is set to custom, indicate the name of the seg_file [NULL]

	-m METAFILE, --metafile=METAFILE
		If --dataset is set to custom, indicate the meta RDS file path [NULL]

	-h, --help
		Show this help message and exit
```

### plotSFC.R
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

### assignOncocode.R
* assignOncocode.R uses the datasets built into RcnvML to categorize each outputed SFC png file into a directory for their oncocode. This is meant to facilitate the directory structure needed to `pickle` the data, done with `common/pickle_hilbert.py` script.
```
Options:
	-c CNTYPE, --cntype=CNTYPE
		Copy number to plot: TCN or ASCN [ASCN]
	-p PDIR, --pdir=PDIR
		Path to seg files [dir]/[dataset]/input/[seg_file] [/mnt/work1/users/pughlab/projects/cancer_cell_lines]
	-s SFC, --sfc=SFC
		Space-filling curve to use, 'sweep' or 'hilbert' [hilbert]
	-d DATASET, --dataset=DATASET
		Dataset to use, 'ccl_aggregate' or 'TCGA' [ccl_aggregate]
	-v, --verbose
		Print extra output [FALSE]
	-h, --help
		Show this help message and exit
```

## Main Figures

## Supplementary Figures
### arrowsSFC.R 
This piece of code is used to generate **Supplementary Figure 1: Mapping of space filling curves**.

Using a 4th and 6th orders/level, it repurposes the `HilbertCurve::HilbertCurve()` function to generate a directional graph of either a HilbertCurve, a sweep, scan, random, diagonal or Morton SFC vurve. This adapted function allows for the input of a different positional space (@POS), which allows for the various SFC curves to be used. If no `sfc_pos` positional space is given, it defaults to the innate Hilbert Curve space.

### fillsSFC.R 
This piece of code is used to generate **Supplementary Figure 2: Mapping the 2D to 1D relation of space filling curves**.

Using a 6th orders/level, it repurposes the `HilbertCurve::HilbertCurve()` function to generate the various SFC curves. Then, by simulating a blocked out region in the 2D space, it will represent the corresponding 1D regions that is reflected by that blocked out region for three conditions: A row block, a column black, and a square block.

### compareSFC.R 
This piece of code is used to generate **Supplementary Figure 3: Comparison between Hilbert and Sweep space filling curves**.

Using several different orders/levels (4, 6, 8), it will generate a Hilbert and Sweep space-filling curve using the standard hg19 genome, colored by chromosomes. It will then run a Euclidean distance between adjacent regions in the 2D space for both SFC.

