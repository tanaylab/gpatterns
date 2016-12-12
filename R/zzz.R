.onLoad <- function(libname, pkgname) {
    gpatterns.set_parallel(parallel::detectCores())
    ggplot2::theme_set(ggplot2::theme_bw(base_size = 15))
    ggplot2::theme_update(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))    
}