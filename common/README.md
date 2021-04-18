# R code for cnvML
## Overview

| Code | Section | Figure | Description |
|------|------|------|------|
| <a href="#pickle_hilbertpy">pickle_hilbert.py</a> | Preprocessing |  | Shuffles, resizes, and pickles the input images |
| <a href="#build_cnn_modelpy">build_cnn_model.py</a> | Preprocessing |  | Builds CNN models and evaluates performance |


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

### build_cnn_model.py
 * build_cnn_model.py is designed to work with the outputted pickle files from `common/pickle_hilbert.py`. It runs the four main wrapper functions: 
   * `load_data.readPickle`: to read in the pickle files for X, Xids, and y
   * `load_data.balanceAndFormatData`: to balance upsample and downsample some samples in the groups to ensure a base level of samples per group (target is median number of samples across all groups, minimum group inclusion is median/4), convert the y values (oncocode dir names) to one-hot encoding, run a train test split (0.2 fraction, seed=1234), and to normalize all X values to a 0-1 scale (from 0-255)
   * `model.buildModel`: Builds the CNN model using keras tensorflow and ADAM optimizer. It will fit and evaluate the data before saving the model to an `.h5` file. It will also outplot the loss_accuracy (cnn_performance.png)
   * `anal.spotcheckModel`: Will take int he model built from buildModel (or read it in if it already exists), and it will plot the confusion matrix for the test class (20% sectioned off from balanceAndFormatData, plot the F1 scores per class, and output those scores as a csv.
```
build_cnn_model.py -m <model> -l <lr> -s <sfc> -c <cntype> -d <dir>

Options:
	-c CNTYPE, --cntype=CNTYPE
		Copy number to plot: TCN or ASCN
	-d DATASET, --dataset=DATASET
		TCGA or CCL
	-s SFC, --sfc=SFC
		Space-filling curve to use, 'sweep' or 'hilbert'
	-l LR, --lr=LR
		Learning rate to pass in for model training (0.001)
	-m MODEL, --model=MODEL
		Which prebuilt CNN model to use (model4)

```
## Main Figures

## Supplementary Figures
