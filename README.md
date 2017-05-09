# The gpatterns Package - Analysis of Methylation Patterns 

This package provides a set of tools for the analysis of methylation patterns.
it includes basic pipeline to import data from raw fastq, as well as functions 
to analyse average methylation and epipolymorphism of multiple samples. 

### Code
Source code can be found at: https://bitbucket.org/tanaylab/gpatterns


### Installation 
#### Requirements 
- _Perl_
- _Python_ (>=2.7)
- _SAMtools_ 
- R packages:
    * _devtools_.
    * _misha_(>=3.5.6).


#### Installing misha package:
```
#!r
install.packages("http://www.wisdom.weizmann.ac.il/~aviezerl/gpatterns/misha_3.5.6.tar.gz", repos=NULL) # Download and install misha package
```

#### Installing gpatterns package:
Download and install *gpatterns*: 
```
#!r
devtools::install_bitbucket("tanaylab/gpatterns", ref='default')
library(gpatterns)
```

#### Using the package
Please refer to the package vignettes for usage and workflow. 
```
#!r
browseVignettes('import') 
browseVignettes('analysis')
```