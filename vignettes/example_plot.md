---
title: "Example `refflat` Plot"
author: "Anthony Aylward"
date: "2019-03-11"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Example `refflat` plot}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




```r
library(RColorBrewer)
library(refflat)
palette <- brewer.pal(3, "Paired") 
par(mar = c(2,4,0,0), mfcol = c(2, 1))
plot(
  t2d_fto[["POS"]],
  -1 * log10(t2d_fto[["PVALUE"]]),
  pch = 21,
  col = palette[2],
  bg = palette[1],
  xlab = "",
  xaxt = "n",
  ylab = "-log10 p-value"
)
plot_refflat("chr16", t2d_fto[1, "POS"], t2d_fto[nrow(t2d_fto), "POS"])
```

![plot of chunk without_buffer](figure/without_buffer-1.png)


```r
library(RColorBrewer)
library(refflat)
palette <- brewer.pal(3, "Paired") 
par(mar = c(2,4,0,0), mfcol = c(2, 1))
plot(
  t2d_fto[["POS"]],
  -1 * log10(t2d_fto[["PVALUE"]]),
  pch = 21,
  col = palette[2],
  bg = palette[1],
  xlab = "",
  xaxt = "n",
  ylab = "-log10 p-value"
)
plot_refflat(
  "chr16",
  t2d_fto[1, "POS"],
  t2d_fto[nrow(t2d_fto), "POS"],
  buffer = 1e3
)
```

![plot of chunk with_buffer](figure/with_buffer-1.png)

