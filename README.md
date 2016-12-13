[![Travis-CI Build Status](https://travis-ci.org/sdray/adespatial.svg?branch=master)](https://travis-ci.org/sdray/adespatial)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/fh9p9cxd80961036/branch/master?svg=true)](https://ci.appveyor.com/project/sdray/adespatial/branch/master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/adespatial)](http://cran.r-project.org/package=adespatial)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/adespatial)](https://cran.r-project.org/package=adespatial)


# adespatial
Multivariate Multiscale Spatial Analysis

This package contains some new functions and many others that were included in development packages hosted in the sedaR project on [R-Forge](https://r-forge.r-project.org/R/?group_id=195).

For instance, `adespatial` includes the `forward.sel` function (formerly in `packfor`) and all functions of `spacemakeR`.

Installing *adespatial*
-------------
To install the development version from github:

1. Install the release version of `devtools` from CRAN with `install.packages("devtools")`.

2. Make sure you have a working development environment.
    * **Windows**: Install [Rtools](http://cran.r-project.org/bin/windows/Rtools/).
    * **Mac**: Install Xcode from the Mac App Store.
    * **Linux**: Install a compiler and various development libraries (details vary across different flavors of Linux).
    
Then:

```r
library(devtools)
install_github("sdray/adespatial")
```

The stable version can be installed from CRAN using:

```r
install.packages("adespatial")
```

Once installed, the package can be loaded using:

```r
library("adespatial")
```
