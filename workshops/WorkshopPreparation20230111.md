# Preparation for MoTrPAC PASS1B R package workshop

> TLDR; complete the steps in this document *before* the workshop if you want to use 
the MoTrPAC R packages to explore PASS1B data. 

Thank you for your interest in the MoTrPAC Workshop: PASS1B R packages & JupyterHub (1/11/23 at 12:00 PM PT)!

The **first part** of this workshop will focus on the **BIC's interactive JupyterHub Dashboard**,
which allows users to explore PASS1B data and results without any code. **This part of 
the workshop is suitable for people without any programming experience.** However, to follow along
on your own computer, you need an authorized [MoTrPAC Data Portal](https://motrpac-data.org/) account. 

During the **second part** of the workshop, we will work through practical examples of how to
use the **[MotrpacRatTraining6mo](https://motrpac.github.io/MotrpacRatTraining6mo/) and 
[MotrpacRatTraining6moData](https://motrpac.github.io/MotrpacRatTraining6moData/) R packages** to explore
and analyze PASS1B data in R. This part of the workshop is geared towards people who are familiar with R, 
but I will do my best to accommodate new programmers/R users as well. 
**If you are interested in this part of the workshop, please complete the steps in this document beforehand**. 
Time **will not** be allotted at the beginning of the workshop to perform these installation steps. 

### Contents
* [Install R](#install-r)  
* [Install RStudio Desktop](#install-rstudio-desktop)
* [Install the MoTrPAC R packages](#install-the-motrpac-r-packages)
  * [Install MotrpacRatTraining6moData (v1.7.0)](#install-motrpacrattraining6modata-v170)
  * [Install MotrpacRatTraining6mo (v1.4.0)](#install-motrpacrattraining6mo-v140)
* [Install other R packages for the workshop](#install-other-r-packages-for-the-workshop)

### Install R  
Optional if you already have R installed. You must be running at least 3.5.0, but
I recommend a newer version since it's up to 4.2.2. Be warned that you must re-install all
packages every time you update R. 

1. Click on the link corresponding to your operating system to download R 4.2.2.
If you aren't sure which kind of Mac you are using, go to `[Apple Icon] > About This Mac`. 
If it says `Intel` under `Processor`, use the second macOS link. 
    * [macOS 11 (Big Sur) and higher, Apple silicon Macs (M1 and higher)](https://cran.r-project.org/bin/macosx/big-sur-arm64/base/R-4.2.2-arm64.pkg)
    * [macOS 10.13 (High Sierra) and higher, Intel 64-bit (older Macs)](https://cran.r-project.org/bin/macosx/base/R-4.2.2.pkg)
    * [Windows](https://cran.r-project.org/bin/windows/base/R-4.2.2-win.exe)
    * [Linux - follow links here](https://cran.r-project.org/bin/linux/)

2. Click on `.pkg` file in Downloads and follow the Installer prompts to install R 4.2.2. 

### Install RStudio Desktop
Optional if you already have RStudio Desktop installed. 

1. Click on the link corresponding to your operating system to download RStudio Desktop (Version: 2022.12.0+353)
    * [Windows 10/11](https://download1.rstudio.org/electron/windows/RStudio-2022.12.0-353.exe)
    * [macOS 10.15+](https://download1.rstudio.org/electron/macos/RStudio-2022.12.0-353.dmg)
    * [Linux - select installer here](https://posit.co/download/rstudio-desktop/)
  
2. Click on the file in Downloads to install RStudio. 

### Install the MoTrPAC R packages 
**Even if you have previously installed these R packages, I recommend re-installing them to 
ensure you are using the newest versions.** We will be using new functions added on 1/5/23. 

#### Install MotrpacRatTraining6moData (v1.7.0)
Note that this package takes 5-10 minutes to install.
If you encounter any errors, follow the troubleshooting instructions [here](https://github.com/MoTrPAC/MotrpacRatTraining6moData#troubleshooting). 

1. Open RStudio  
2. Copy/paste the following command in the `Console`:  

    ```r
    if (!require("devtools", quietly = TRUE)){
      install.packages("devtools")
    }
    options(timeout=1e5)
    devtools::install_github("MoTrPAC/MotrpacRatTraining6moData")
    ```
  
**If you use this data R package in any capacity, I highly recommend "watching" the corresponding
GitHub repository.** This will send you email notifications about new releases, which
can include changes to the data. In the top-right corner of [this page](https://github.com/MoTrPAC/MotrpacRatTraining6moData), 
choose `Watch` > `Custom` > `Releases`. This requires a GitHub account. You can sign up
for a free account [here](https://github.com/join). 

#### Install MotrpacRatTraining6mo (v1.4.0)
This package should install very quickly if `MotrpacRatTraining6moData` is already installed. 

1. Open RStudio  
2. Copy/paste the following command in the `Console`:  

    ```r
    devtools::install_github("MoTrPAC/MotrpacRatTraining6mo")
    ```

### Install other R packages for the workshop 
To avoid having to install additional dependencies *during* the workshop, install
them ahead of time by copy/paste-ing the following commands in your RStudio `Console`. 
```r
install.packages("foreach")
install.packages("gprofiler2")
install.packages("doParallel")
install.packages("reshape2")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("IHW")
```

**Technical aside:** You may be wondering why you need to install these packages
manually if `MotrpacRatTraining6mo` depends on them, especially when many other packages
were installed when you installed `MotrpacRatTraining6mo`. This is because there
are (at least) two kinds of R package dependencies: `Imports` and `Suggests`. Dependencies
listed under `Imports` must be installed when the R package is installed. Dependencies listed
under `Suggests` are required by lesser-used functions, so to reduce overhead during R package
installation, these dependencies are not automatically installed. Instead, when you call one of 
these lesser-used functions and a dependency is missing, you get an error message instructing
you to install an additional package. If you are curious 
about which packages are required or suggested for `MotrpacRatTraining6mo`, take a look at the
[DESCRIPTION](https://github.com/MoTrPAC/MotrpacRatTraining6mo/blob/main/DESCRIPTION) file. 

That's it! Thank you for getting this far. I look forward to the workshop on 1/11/23! 

---

Author: Nicole Gay  
Updated: 1/5/23
