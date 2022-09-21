# MotrpacRatTraining6mo

## Table of Contents
* [Overview](#overview)
  * [About this package](#about-this-package)
  * [About MoTrPAC](#about-motrpac)
* [Installation](#installation)
* [Getting help](#getting-help)
* [Acknowledgements](#acknowledgements)
* [Data use agreement](#data-use-agreement)
* [Citing MoTrPAC data](#citing-motrpac-data)

## Overview

### About this package 
This package provides functions to fetch, explore, and reproduce the processed data and downstream
analysis results presented in the main paper for the first 
large-scale multi-omic multi-tissue endurance exercise training study conducted 
in young adult rats by the Molecular Transducers of Physical Activity Consortium 
(MoTrPAC). 
**See the [vignette](https://motrpac.github.io/MotrpacRatTraining6mo/articles/MotrpacRatTraining6mo.html) for examples of how to use this package.**
A [bioRxiv](https://www.biorxiv.org/) link to the corresponding 
preprint will be added shortly. 

While some of the functions in this package can be used by themselves, they
were primarily written to analyze data in the 
[MotrpacRatTraining6moData](https://motrpac.github.io/MotrpacRatTraining6moData)
R package. See examples of how these data can be analyzed *without* this package in the 
[MotrpacRatTraining6moData vignette](https://motrpac.github.io/MotrpacRatTraining6moData/articles/MotrpacRatTraining6moData.html).

### About MoTrPAC
MoTrPAC is a national research consortium designed to discover and perform 
preliminary characterization of the range of molecular transducers (the 
"molecular map") that underlie the effects of physical activity in humans. 
The program's goal is to study the molecular changes that occur during and after 
exercise and ultimately to advance the understanding of how physical activity 
improves and preserves health. The six-year program is the largest targeted NIH 
investment of funds into the mechanisms of how physical activity improves health 
and prevents disease. See [motrpac.org](https://www.motrpac.org/) and 
[motrpac-data.org](https://motrpac-data.org/) for more details. 

## Installation
Install this package with `devtools`:
```r
if (!require("devtools", quietly = TRUE)){
  install.packages("devtools")
}
devtools::install_github("MoTrPAC/MotrpacRatTraining6mo")
```

## Getting help 
**See the [vignette](https://motrpac.github.io/MotrpacRatTraining6mo/articles/MotrpacRatTraining6mo.html) 
for examples of how to use this package.**
Still have questions? For questions, bug reporting, and feature requests for this package, please 
[submit a new issue](https://github.com/MoTrPAC/MotrpacRatTraining6mo/issues) 
and include as many details as possible. 

If the concern is related to data provided in the 
[MotrpacRatTraining6moData](https://github.com/MoTrPAC/MotrpacRatTraining6moData)
package, please submit an issue 
[here](https://github.com/MoTrPAC/MotrpacRatTraining6moData/issues) instead. 

## Acknowledgements 
MoTrPAC is supported by the National Institutes of Health (NIH) Common
Fund through cooperative agreements managed by the National Institute of Diabetes and
Digestive and Kidney Diseases (NIDDK), National Institute of Arthritis and Musculoskeletal
Diseases (NIAMS), and National Institute on Aging (NIA). 
Specifically, the MoTrPAC Study is supported by NIH grants U24OD026629 (Bioinformatics Center), 
U24DK112349, U24DK112342, U24DK112340, U24DK112341, U24DK112326, U24DK112331, U24DK112348 (Chemical Analysis Sites), 
U01AR071133, U01AR071130, U01AR071124, U01AR071128, U01AR071150, U01AR071160, U01AR071158 (Clinical Centers), 
U24AR071113 (Consortium Coordinating Center), U01AG055133, U01AG055137 and U01AG055135 (PASS/Animal Sites).

## Data use agreement 
Recipients and their Agents agree that in publications using **any** data from MoTrPAC public-use data sets 
they will acknowledge MoTrPAC as the source of data, including the version number of the data sets used, e.g.:

* Data used in the preparation of this article were obtained from the Molecular Transducers of Physical Activity 
Consortium (MoTrPAC) database, which is available for public access at motrpac-data.org. 
Specific datasets used are [version numbers].

* Data used in the preparation of this article were obtained from the Molecular Transducers of Physical Activity 
Consortium (MoTrPAC) MotrpacRatTraining6moData R package [version number]. 

## Citing MoTrPAC data 
A citation will be added when the preprint is available on [bioRxiv](https://www.biorxiv.org/). 
