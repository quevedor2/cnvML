#install.packages("devtools")
library(devtools)

pkg=file.path("~/git", "cnvML", "RcnvML")
setwd(pkg)


#use_description()

#### Assembling data ####
# usethis::use_testthat()               # unit testing
# usethis::use_vignette("vignette")     # vignette documentation
# usethis::use_data_raw()               # raw data

usethis::use_package("assertthat")
usethis::use_package("utils")
usethis::use_package("stats")
usethis::use_package("methods")
usethis::use_package("HilbertCurve")
usethis::use_package("GenomicRanges")
usethis::use_package("intervals")
usethis::use_package("GenomeInfoDb")
usethis::use_package("BSgenome.Hsapiens.UCSC.hg19")
usethis::use_package("ComplexHeatmap")
usethis::use_package("grid")
usethis::use_package("abind")
usethis::use_package("grDevices")

#### Building ####
devtools::load_all()
devtools::document(pkg)
devtools::check(pkg)
#devtools::build_vignettes(pkg)

devtools::build(pkg)
devtools::install(pkg)
# devtools::install_github("quevedor2/aneuploidy_score")
devtools::install_github("quevedor2/RocheTest", ref = "dev")
