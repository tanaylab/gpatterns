# The gpatterns pacakge - Analysis of methylation patterns#

This package provides a set of tools for the analysis of methylation patterns. 

### Installation 
#### Requirements 
- _Perl_
- R packages:
    * remotes.
    * _misha_.


#### Installing misha package:
```
#!r
install.packages("http://www.wisdom.weizmann.ac.il/~aviezerl/gpattrerns/misha_3.5.4.tar.gz", repos=NULL) # Download and install misha package
```

#### Installing gpatterns package:
Download and install *gpatterns*: 
```
#!r
remotes::install_bitbucket("tanaylab/gpatterns", ref='default', vignette = TRUE)
library(gpatterns)
```

#### Using the package
Please refer to the package vignette for usage and workflow.
```
#!r
browseVignettes('gpatterns') 
```