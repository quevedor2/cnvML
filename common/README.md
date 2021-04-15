# R code for cnvML
## Overview

| Code | Section | Figure | Description |
|------|------|------|------|
| <a href="#pickle_hilbertpy">pickle_hilbert.py</a> | Preprocessing |  | Shuffles, resizes, and pickles the input images |


## Pre-processing
### pickle_hilbert.py
* pickle_hilbert.py is designed to work the outputted directory structure from `R/assignOncocode.R`. It has 4 main functions: 1) to search each directory within the path for pngs and classify them according to their directory, 2) to reshape the image to a standard size (300x300), 3) to shuffle the order of the images (`X`) while maintaining their associated label (`y`), and 4) to output a pickle file that will be used for downstream tensorflow applications. As such, it will output 4 files, `X.pickle`, `Xids.pickle`, `y.pickle`, and `yids.pickle` to the `[dir]/[analysis]/data/[sfc]/[cntype]` directory.

```
pickle_hilbert.py -a <analysis> -s <sfc> -c <cntype> -d <dir>

Options:
	-a CNTYPE, --cntype=CNTYPE
		Copy number to plot: TCN or ASCN
	-d DIR, --dir=DIR
		Path to png oncocode separated files [dir]/[analysis]/data/[sfc]/[cntype]
	-s SFC, --sfc=SFC
		Space-filling curve to use, 'sweep' or 'hilbert'
	-a ANALYSIS, --analysis=ANALYSIS
		Dataset to use, 'CCL' or 'TCGA'

```

## Main Figures

## Supplementary Figures
