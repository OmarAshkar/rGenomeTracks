[![](https://img.shields.io/badge/devel%20version-0.99.3-blue.svg)](https://github.com/OAshkar/rGenomeTracks)
[![](https://img.shields.io/badge/download-15/total-blue.svg)](https://bioconductor.org/packages/stats/bioc/rGenomeTracks)

rGenomeTracks is an R package the leverage the power of pyGenomeTracks to create quick and tidy epigenomic visualizations in R.


## Installation
Installation from bioconductor

```r
BiocManager::install("rGenomeTracks")
```
or from github 

```r
## install.packages("devtools")
devtools::install_github("OAshkar/rGenomeTracks")
```


## Basic Usage

```r
library(rGenomeTracks)

# Get track directory
 tads_dir <- system.file("extdata", "tad_classification.bed", package = "rGenomeTracks")
 
 # Create TADS track
 tads <- track_domains(
   file = tads_dir, border_color = "black",
   color = "none", height = 5,
   line_width = 5,
   show_data_range = FALSE,
 )

#Plot the track
par(mar = c(1,1,1,1))
plot_gtracks(tads, chr = "X", start = 25 * 10^5, end = 31 * 10^5)
```

Please visit the package's vignette for more examples.

## Cite 

``` r
citation(rGenomeTracks)
```

## Suggestions & Bug Reports
Feel free to open a new issue here: https://github.com/OAshkar/rGenomeTracks/issues
