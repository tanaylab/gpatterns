############################################################################################
#' Effiecient implementation of kmeans++ algorithm
#'
#' @param tab table
#' @param K K
#' @param id_column id_column
#' @param allow_nas allow_nas
#' @param verbose verbose
#' @param fn fn
#' @param tab_file tab_file
#' @param TGL.kmeans.bin TGL.kmeans.bin
#'
#' @return
#' @export
#'
#' @examples
TGL_kmeans <- function(tab,
                        K,
                        id_column = FALSE,
                        allow_nas = TRUE,
                        verbose = FALSE,
                        fn = NULL,
                        tab_file = FALSE,
                        TGL.kmeans.bin = sprintf("%s/%s",
                                                 system.file("bin", package="gpatterns"),
                                                 'TGLKMeans_static')){
    if(!tab_file){
        if (is.null(fn)){
            fn <- tempfile()
        }
        if (id_column){
            id.col.name <- colnames(tab)[1]
            colnames(tab)[1] <- 'id'
        } else {
            tab <- data.frame(id = 1:nrow(tab), tab)
        }

        message("writing table...")
        write.table(tab, fn, sep="\t", quote=F, row.names=FALSE)
    } else {
        fn <- tab
    }

    message("clustering...")
    system(
        qq('@{TGL.kmeans.bin} @{fn} @{K} euclid -allow_nas=@{as.numeric(allow_nas)}'),
        ignore.stdout = !verbose,
        ignore.stderr = !verbose)

    tab_k <- fread(paste(fn, "kclust", sep="."), header=T, colClasses=c('character', 'numeric'))
    if (id_column){
        colnames(tab_k)[1] <- id.col.name
    }
    return(tab_k)
}


############################################################################################
#' Alternative interface for TGL_kmeans
#' Returned value is more similar to native R kmeans
#' Allows keeping the output of the kmeans exec within the R return value
#'
#'
#' @param data
#' @param keep.log
#' @param tmpdir
#' @param tmpfn
#' @param TGL.kmeans.bin
#'
#' @return
#' @export
#'
#' @examples
TGL_kmeans2 <- function(data,
                        k,
                        keep.log=FALSE,
                        tmpdir=NULL, tmpfn=NULL,
                        TGL.kmeans.bin = "inst/bin/TGLKMeans_static")
{
    library(data.table)

    tmpfile <- ''
    if (is.null(tmpdir)) {
        tmpdir <- tempdir()
    }
    if (is.null(tmpfn)) {
        tmpfile <- tempfile(pattern='tgl_kmeans.', tmpdir=tmpdir)
    } else {
        tmpfile <- file.path(tmpdir, tmpfn)
    }

    write.table(data, tmpfile, sep='\t', quote=FALSE, col.names=FALSE)

    km_log <- system2(TGL.kmeans.bin, c(tmpfile, k, 'euclid', '-allow_nas=1'), stdout=keep.log, stderr=keep.log)
    kclust <- as.data.frame(fread(paste0(tmpfile, '.kclust'), sep='\t', header=TRUE))
    center <- as.data.frame(fread(paste0(tmpfile, '.center'), sep='\t', header=FALSE))

    kclust$clust <- kclust$clust + 1

    km <- list()

    km$cluster <- kclust$clust
    names(km$cluster) <- kclust$id
    km$centers <- center[,-1]
    km$size <- tapply(km$clust, km$clust, length)

    if (keep.log) {
        km$log <- km_log
    }

    return(km)
}
