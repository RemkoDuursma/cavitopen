# CAVITOPEN

An R implementation of the CAVITOPEN model of Herve Cochard et al.


## Installation

On Windows, you must have [Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed. Also, install the `devtools` R package. Then,

```
devtools::install_github("remkoduursma/cavitopen")
```

## Use

This following code reproduces the example in the original spreadsheet.

```
library(cavitopen)
run1 <- cavitopen()
plot(run1)
```
