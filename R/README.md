# R code for cnvML
## Overview

| Code | Figure | Description |
|------|------|------|
| compareSFC.R | Fig S1 | Sweep vs Hilbert SFC |

## Pre-processing

## Main Figures

## Supplementary Figures
### compareSFC.R 
This piece of code is used to generate **Supplementary Figure 1: Comparison between Hilbert and Sweep space filling curves**.

Using several different orders/levels (4, 6, 8), it will generate a Hilbert and Sweep space-filling curve using the standard hg19 genome, colored by chromosomes. It will then run a Euclidean distance between adjacent regions in the 2D space for both SFC.

