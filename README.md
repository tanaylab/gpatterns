The gpatterns Package - Analysis of Methylation Patterns
========================================================

This package provides a set of tools for the analysis of methylation patterns. it includes basic pipeline to import data from raw fastq, as well as functions to analyse average methylation and epipolymorphism of multiple samples.

### Code

Source code can be found at: <https://bitbucket.org/tanaylab/gpatterns>

### Installation

#### Requirements

-   *Perl*
-   *Python* (&gt;=2.7)
-   *SAMtools*
-   R packages:
    -   *devtools*.
    -   *misha*(&gt;=3.5.6).

#### Installing misha package:

``` r
# Download and install misha package
install.packages("http://www.wisdom.weizmann.ac.il/~aviezerl/gpatterns/misha_3.5.6.tar.gz", repos=NULL) 
```

#### Installing gpatterns package:

Download and install *gpatterns*:

``` r
devtools::install_bitbucket("tanaylab/gpatterns", ref='default')
library(gpatterns)
```
Important: 
Until the package would be open for everyone, use: 
``` r
devtools::install_bitbucket("tanaylab/gpatterns", ref='default', auth_user='username', password='password')
library(gpatterns)
```
Where 'username' is your bitbucket username and 'password' is your bitbucket password.

#### Using the package

Please refer to the package vignettes for usage and workflow.

``` r
browseVignettes('import') 
browseVignettes('analysis')
```