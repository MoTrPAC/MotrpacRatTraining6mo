# Notes for developers

## Load the package in development mode
Set the working directory to this repository:  
```r
setwd("~/src/MOTRPAC/MotrpacRatTraining6mo")
library(devtools)
```
Load all of the functions defined in `./R`:  
```r
load_all()
```

## Create a test
Use `usethis::use_r()` and `usethis::use_test()`. 
See details [here](https://r-pkgs.org/testing-basics.html#create-a-test). 

## Create a vignette 
Use `usethis::use_vignette("vignette-title")`, edit the .Rmd in `vignettes/`, 
and build it with `devtools::build_rmd()`. See details [here](https://r-pkgs.org/vignettes.html). 
