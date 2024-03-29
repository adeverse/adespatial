# adespatial <img src='man/figures/logo.svg' align="right" height="139" />

[![R-CMD-check](https://github.com/adeverse/adespatial/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/adeverse/adespatial/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version-ago/adespatial)](http://cran.r-project.org/package=adespatial)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/adespatial)](https://cran.r-project.org/package=adespatial)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Multivariate Multiscale Spatial Analysis

---------------------------

**Please note! Since January 2024, this repository has belonged to the *adeverse* organization.**
To avoid confusion, we strongly recommend updating any existing local clones to point to the new 
repository URL. You can do this by using `git remote` on the command line:

`git remote set-url origin git@github.com:adeverse/adespatial.git`

or 

`git remote set-url origin https://github.com/adeverse/adespatial.git`

---------------------------


This package contains some new functions and many others that were included in development packages hosted in the sedaR project on [R-Forge](https://r-forge.r-project.org/R/?group_id=195).

For instance, `adespatial` includes the `forward.sel` function (formerly in `packfor`) and all functions of `spacemakeR`. To have an overview of the package, read the vignette after installing the package by:

```r
vignette("tutorial", package = "adespatial")
```

Installing the development version of `adespatial`
-------------

- Install the release version of `remotes` from CRAN with `install.packages("remotes")`.

- Make sure you have a working development environment.
    * **Windows**: Install [Rtools](http://cran.r-project.org/bin/windows/Rtools/).
    * **Mac**: Install Xcode from the Mac App Store.
    * **Linux**: Install a compiler and various development libraries (details vary across different flavors of Linux).
    
Then:

```r
remotes::install_github("adeverse/adespatial")
```

If you do not wish to install the development environments Rtools (Windows) / XCode (Mac), we can supply binary packages of the development version of `adespatial` on request. 


Installing the stable version of `adespatial`
-------------

```r
install.packages("adespatial")
```


Loading `adespatial`
-------------

```r
library("adespatial")
```
