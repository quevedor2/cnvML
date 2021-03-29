


<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
    </li>
    <li>
      <a href="#usage">Usage</a>
      <ul>
        <li><a href="#importing-biomedical-data">Importing biomedical data</a></li>
        <li><a href="#data-visualization">Data visualization</a></li>
        <li><a href="#group-comparisons">Group comparisons</a></li>
      </ul>
    </li>
  </ol>
</details>




## About The Project

`RocheTest` package is part of the Computational Biology take home test from Roche. This project has 3 main requirements:
1. Write a function to import biomedical data (TCGA-CDR)
2. Write a function to visualize this data
3. Write a function to conduct statistical test for difference between 2 groups 

## Getting Started

`RocheTest` requires that several packages are installed. However, all 
dependencies are available from CRAN.

```{r install-pkg, eval=FALSE, results='hide'}
library(devtools)
devtools::install_github('quevedor2/RocheTest')
```

## Usage
### Importing biomedical data
The Cancer Genome Atlas - Clinical Data Resource (TCGA-CDR) can be downloaded from [https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81](https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81). The first table in the excel spreadsheet can be saved as csv file and imported using the `loadCDR()` function. Additionally, this data is already preloaded and saved as a data object.

As part of this project, the saved data object and default parameters for this function have a preset filter for patients with a known tumor status.

```{r read-in-dataset, eval=FALSE, results='hide'}
library(RocheTest)
cdr_path <- '/path/to/tcga-cdr.csv'
loadCDR(cdr_path = cdr_path, filter_col = 'tumor_status', filter_cat = 'WITH TUMOR')
```

To access the prebuilt data:
```{r load-data, eval=TRUE, results='hide'}
library(RocheTest)
data('tcgacdr')
```

### Data visualization
The second part of this project is to create a function to visualize and compare specific metadata column for a set of samples. To approach this function, `RocheTest` implements ggplot boxplot from a melted data structure. While more information can be found in the manual, it should be noted that either `patients` or `cancer_type` parameter must be given. Additionally, only one of the two must be given, the other one is set to NULL by default.

To get a list of available cancer_types or patients, you can explore the `tcgacdr` object:
```{r explore-cdr, eval=TRUE}
# cancer_type
head(table(tcgacdr$type),5)

# patients
head(tcgacdr$bcr_patient_barcode, 5)
```

To visualize the difference between Progression Free Interval (PFI) between several cancer types for Stage I and Stage III tumor stages:
```{r viz-pfi-ctype, eval=TRUE}
vizCDR(tcgacdr, cancer_type=c('ACC', 'BLCA', 'BRCA'),
       metadata_col = 'PFI.time', meta_col = 'ajcc_pathologic_tumor_stage',
       metadata_levels =  c('Stage III', 'Stage IV'))
```


To visualize the difference between Age for a subset of `patients` with a new tumor event::
```{r viz-age-tumorevents, eval=TRUE}
vizCDR(tcgacdr, patients=sample(tcgacdr$bcr_patient_barcode, size=1000),
       metadata_col = 'age_at_initial_pathologic_diagnosis', 
       meta_col = 'new_tumor_event_type',
       metadata_levels =  c('Distant Metastasis', 'Locoregional Recurrence'))
```

### Group comparisons
The final part of this project is to create a function to compare 2 groups based on a specific metadata column. Given that the standard approach for a 2-group comparison is a parametric or a nonparametric test, the t-test and wilcoxon rank-sum tests are implemented in the `compareGroups()` function by default. Additionally, as a non-parametric is the safer test for datasets of an unknown distribution, it is the default option for `stat` parameter.

For this example, I am running a parametric two-tailed t-test to look for differences between age for glioblastoma (GBM) compared to lung adenocarcinoma (LUAD):
```{r test-age-param, eval=TRUE}
age_test <- compareGroups(tcgacdr, 
       metadata_col = 'age_at_initial_pathologic_diagnosis', 
       meta_col = 'type',
       metadata_levels =  c('GBM', 'LUAD'), 
       stat='ttest')

print(paste0("LUAD mean age: ", mean(age_test$data$LUAD, na.rm=T)))
print(paste0("LUAD mean age: ", mean(age_test$data$GBM, na.rm=T)))
print(paste0("t-statistic: ", age_test$stat$statistic))
print(paste0("p-val: ", age_test$p))
```

To compare values using a non-parametric test, I am using the `new_tumor_event_dx_days_to` column as it is more count-based data (of days) and represents more of a Poisson distribution over Normal (evidenced by difference between median and mean). Using the wilcoxon rank-sum test, I can compare the differences between Asians and African Americans:
```{r test-dx-nonparam, eval=TRUE}
lvl <- c('ASIAN', 'BLACK OR AFRICAN AMERICAN')
dx_test <- compareGroups(tcgacdr, 
       metadata_col = 'new_tumor_event_dx_days_to', 
       meta_col = 'race',
       metadata_levels =  lvl, 
       stat='wilcox')

qstat <- matrix(c(median(dx_test$data$ASIAN, na.rm=T),
                  median(dx_test$data$`BLACK OR AFRICAN AMERICAN`, na.rm=T),
                  mean(dx_test$data$ASIAN, na.rm=T),
                  mean(dx_test$data$`BLACK OR AFRICAN AMERICAN`, na.rm=T)),
                byrow=T, ncol=2)
colnames(qstat) <- lvl
rownames(qstat) <- c('median', 'mean')

qstat
print(paste0("W-statistic: ", dx_test$stat$statistic))
print(paste0("p-val: ", dx_test$p))
```

To complement the default ttest and u-test approach, `compareGroups()` allows for users to specific their own function using the `stat='custom'` method.  The requirements for this approach is that the `statfun` must be a function with an `x` and `y` parameter that returns a `list` object with a `p.value` element. 

This customization method is highlighted here where I implement a one-sided t-test to assess whether MALE's have a longer overal survival time than FEMALEs:

```{r test-os-custom, eval=TRUE}
ttestOne <- function(x,y){
  t.test(x,y, alternative = 'greater')
}

os_test <- compareGroups(cdr=tcgacdr, 
       metadata_col = 'OS.time', 
       meta_col = 'gender',
       metadata_levels =  c('MALE', 'FEMALE'), 
       stat='custom', statfun=ttestOne)

print(paste0("MALE mean OS.time: ", mean(os_test$data$MALE, na.rm=T)))
print(paste0("FEMALE mean OS.time: ", mean(os_test$data$FEMALE, na.rm=T)))
print(paste0("t-statistic: ", os_test$stat$statistic))
print(paste0("p-val: ", os_test$p))
```
