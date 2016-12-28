# Utility functions

#' @export
qq <- GetoptLong::qq

#' @export
qqv <- function(text) {GetoptLong::qq(text, envir=parent.frame(), collapse=FALSE)}

#' @export
fread <- partial(data.table::fread, data.table=FALSE)

#' @export
fwrite <- data.table::fwrite

#############################################################################
.is_tidy_cpgs <- function(obj){
    return('data.frame' %in% class(obj) &&
               all(
                   c(
                       "read_id",
                       "chrom",
                       "start",
                       "end",
                       "strand",
                       "umi1",
                       "umi2",
                       "insert_len",
                       "num",
                       "cg_pos",
                       "meth",
                       "qual"
                   ) %in% colnames(obj)
               ))
}

#############################################################################
#' @export
comify <- scales::comma

#############################################################################
#' @inheritParams gpatterns::gcluster.run2
#' @export
gcluster.run3 <- partial(gcluster.run2, script = .gpatterns.sg_script)

#############################################################################
.gpatterns.base_dir <- function(track)
{
    map(track, ~  c(gdir.cwd(), strsplit(.x, '.', fixed=TRUE)[[1]])) %>%
        map(~ do.call(file.path, as.list(.x))) %>%
        as_vector %>%
        return()
}

########################################################################
#' @export
gpatterns.track_exists <- function(track){
    dir.exists(paste0(.gpatterns.base_dir(track), '/tidy_cpgs'))
}

########################################################################
.gpatterns.intervals2files <- function(intervals, files, dir){
        intervals %>% unite('coord', chrom:end, sep='_') %>%
            mutate(coord = paste0(dir, coord, '.tcpgs.gz')) %>%
            filter(coord %in% files) %>% .$coord
}

#############################################################################
.gpatterns.force_chromosomes <- function(df){
    f <- df$chrom %in%  as.character(gintervals.all()$chrom )
    if (!all(f)){
        chroms <- df$chrom[!f] %>% unique %>% paste(collapse=', ')
        warning(qq("removed the following chromosomes that did not exist in genome db: @{chroms}"))
    }
    return(df %>% filter(f))
}

#############################################################################
#' Set parallel threads
#'
#' @param thread_num number of threads. use '1' for non parallel behaviour
#'
#' @return None
#'
#' @examples
#' gpatterns.set_parallel(8)
#'
#' @export
gpatterns.set_parallel <- function(thread_num) {
    if (1 == thread_num) {
        options(gpatterns.parallel = FALSE)
    } else {
        doMC::registerDoMC(thread_num)
        options(gpatterns.parallel = TRUE)
        options(gpatterns.parallel.thread_num = thread_num)
    }
}
