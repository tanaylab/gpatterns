
#' Effiecient implementation of kmeans++ algorithm
#'
#' @param tab table
#' @param K K
#' @param id_column id_column
#' @param allow_nas allow_nas
#' @param verbose verbose
#' @param fn fn
#' @param tab_file tab_file
#' @param order_rows order_rows
#' @param clust_order_func clust_order_func
#' @param row_order_func row_order_func
#' @param method method
#' @param row_order_column row_order_column 
#' @param TGL_kmeans_bin TGL_kmeans_bin
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
                        order_rows = TRUE,
                        clust_order_func = partial(mean, na.rm=T),
                        row_order_func = NULL,
                        method='ward.D2',
                        row_order_column = FALSE,
                        TGL_kmeans_bin = sprintf("%s/%s",
                                                 system.file("bin", package="gpatterns"),
                                                 'TGLKMeans_static')){
    if(!tab_file){
        if (is.null(fn)){
            fn <- tempfile()
            on.exit(system(qq('rm -f @{fn}')))
        }
        if (id_column){
            id.col.name <- colnames(tab)[1]
            colnames(tab)[1] <- 'id'
        } else {
            tab <- data.frame(id = 1:nrow(tab), tab)
        }

        message("writing table...")        
        fwrite(tab, fn, row.names=FALSE, sep="\t", quote=FALSE, na='NA')
    } else {
        fn <- tab
    }

    message("clustering...")
    system(
        qq('@{TGL_kmeans_bin} @{fn} @{K} euclid -allow_nas=@{as.numeric(allow_nas)}'),
        ignore.stdout = !verbose,
        ignore.stderr = !verbose)

    tab_k <- fread(paste(fn, "kclust", sep='.'), header=T, colClasses=c('character', 'numeric'))    
    tab$clust <- tab_k$clust    

    summarise_clust <- function(x, clust_order_func) {
        tibble(clust = x$clust[1], m = x %>% select(-clust) %>% as.matrix() %>% clust_order_func())
    }
    clust_ord <- tab %>% select(-id) %>% group_by(clust) %>% do({summarise_clust(., clust_order_func)}) %>% ungroup %>% arrange(m) %>% mutate(new_clust = 1:n()) %>% select(-m)
    tab <- tab %>% left_join(clust_ord, by='clust') %>% select(-clust) %>% rename(clust=new_clust) %>% arrange(clust)
    
    if (order_rows){
        if (!is.null(row_order_func)) {        
            tab <- tab %>% plyr::ddply(.(clust), function(x) {
                    x$m <- apply(x %>% select(-id, -clust), 1, row_order_func)
                    x %>% arrange(m) %>% mutate(ord=1:nrow(x)) %>% select(-m)} ) %>%
                select(id, clust, ord, everything())
            tab <- tab %>% group_by(clust) %>% mutate(ord = 1:n()) %>% ungroup %>% select(id, clust, ord, everything()) %>% arrange(clust, ord)
        } else {
            tab <- tab %>% group_by(clust) %>% do(.hclust_order(., c('id', 'clust'), NULL, NULL, tidy=FALSE, method=method)) %>% ungroup %>% arrange(clust, ord)
        }
    } else {
        tab <- tab %>% group_by(clust) %>% mutate(ord = 1:n()) %>% ungroup %>% select(id, clust, ord, everything()) %>% arrange(clust, ord)
    }
   
    if (id_column){
        colnames(tab)[1] <- id.col.name
    }

    if (!row_order_column){
        tab <- tab %>% select(-ord)
    }

    return(tab)
}



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
                        TGL.kmeans.bin = sprintf("%s/%s",
                                                 system.file("bin", package="gpatterns"),
                                                 'TGLKMeans_static'))
{

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
