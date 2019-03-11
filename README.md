# refflat
This package is designed to plot the locations of genes within a small genomic
interval. It is useful, for example, to draw genes below a GWAS signal
for context. See also this [example](https://github.com/anthony-aylward/islet-cytokines-outline/blob/master/nfatc2-ENSG00000101096.pdf) of a plot created using `refflat`.

## Installation
```r
library(devtools)
install_github("anthony-aylward/refflat")
```

## Usage

For a very basic example of what `refflat` does, try the following:

```r
library(refflat)
plot_refflat("chr16", 535e5, 542e5)
```

For more details, see the included [vignette](vignettes/example_plot.md)
