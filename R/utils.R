# Utility functions

#' @export
qq <- GetoptLong::qq

#' @export
qqv <- partial(GetoptLong::qq, collapse = TRUE)

#' @export
fread <- partial(data.table::fread, data.table=FALSE)


#############################################################################
comify <- scales::comma

#############################################################################
.gpatterns.base_dir <- function(track)
{
    base_path <- c(gdir.cwd(), strsplit(track, '.', fixed=TRUE)[[1]])
    return(do.call(file.path, as.list(base_path)))
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