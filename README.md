[![CRAN_Release_Badge](http://www.r-pkg.org/badges/version-ago/rbioacc)](http://cran.r-project.org/package=rbioacc)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/rbioacc)](https://cran.r-project.org/package=rbioacc)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI:10.1016/j.ecoenv.2022.113875](http://img.shields.io/badge/DOI-10.1016/j.ecoenv.2022.113875-orange.svg)](https://doi.org/10.1016/j.ecoenv.2022.113875)

# rbioacc

Inference and Prediction of ToxicoKinetic (TK) Models


### Little hack

To load all internal function of a package during dev: `devtools::load_all()`

### A lighter package build

To make the package lighter, we have to remove the vignettes: see file `.Rbuildignore`

### Error to recompile during package development

Sometimes, there is an Error to recompile during development after change of .stan files.
A solution is to remove the `rbioacc` folder in R repository of the win-library (see the path written in the error message).

An other solution is to build the package from the terminal using `R CMD -preclean INSTALL rbioacc` from parent directory of `rbioacc`.



### Note 

- S3 Object System: http://adv-r.had.co.nz/S3.html
- Google's R Style Guide: https://google.github.io/styleguide/Rguide.html
- testthat: https://github.com/r-lib/testthat
- covr: https://github.com/r-lib/covr
- to add a package `xxr`: `usethis::use_package("xxr")`
- to add data set `datar`: `usethis::use_data(datar)`
