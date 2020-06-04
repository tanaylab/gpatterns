#!/usr/bin/env Rscript
library(misha)
suppressPackageStartupMessages(suppressWarnings(library(gpatterns)))
library(optparse)

########################################################################
option_list = list(
  make_option(c('-i', "--input"), action="store", default='-', type='character',
              help="tidy cpgs input file. a value of '-' indicates that the file would be read from stdin"),
  make_option(c('-o', "--output"), action="store", default='-', type='character',
              help="output file. a value of '-' indicates that the file would be read from stdout"), 
  make_option(c('-f', "--frag_intervs"), action="store", type='character', default=NA,
              help="intervals set of the fragment"),   
  make_option(c("--groot"), action="store", default=NA, type='character',
              help="root of misha db"),  
  make_option(c("--maxdist"), action="store", type='numeric', default=NA,
              help="maximal distance from fragments"),
  make_option(c("--rm_off_target"), action="store_true", type='logical', default=FALSE,
              help="remove reads with distance > maxdist from frag_intervs (withut this option those reads would be left unchanged)"))  

parser <- OptionParser(usage = "%prog -i <input> -o <output> -f <fragment intervals> --groot <misha root>", 
        description = "Change tidy_cpgs coordinates to fragments coordinates (e.g. MSP1 fragments)",
        option_list=option_list)

opt <- parse_args(parser, args=commandArgs(TRUE))

 if(is.na(opt$frag_intervs) || is.na(opt$groot)) {   
    message(parser@usage)
    quit()
 }

 if (opt$input == '-'){
    opt$input <- 'cat /dev/stdin'
}

if (opt$output == '-'){
    opt$output <- ''
}
 
misha::gsetroot(opt$groot)
calls <- dplyr::as_tibble(gpatterns::fread(opt$input, header=TRUE, showProgress = FALSE, verbose=FALSE))
if (nrow(calls) == 0){
    res <- calls
} else {
    res <- gpatterns::gpatterns.adjust_read_pos(calls, frag_intervs=opt$frag_intervs, maxdist=opt$maxdist, rm_off_target=opt$rm_off_target)
}
data.table::fwrite(res, opt$output, sep=',', quote=FALSE)
