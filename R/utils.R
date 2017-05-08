# Utility functions ------------------------------------------------

#' @export
qq <- GetoptLong::qq

#' @export
qqv <- function(text) {GetoptLong::qq(text, envir=parent.frame(), collapse=FALSE)}

#' @export
fread <- partial(data.table::fread, data.table=FALSE)

#' @export
fwrite <- data.table::fwrite

#' @export
comify <- scales::comma


#' Smoothscatter that works with the pipe operator
#' @export
p_smoothScatter <- function(.data, .x, .y, xlab = NULL, ylab = NULL, ...){
    x <- lazyeval::lazy_eval(lazyeval::lazy(.x), .data)
    y <- lazyeval::lazy_eval(lazyeval::lazy(.y), .data)
    xlab <- xlab %||% as.character(lazyeval::lazy(.x)$expr)
    ylab <- ylab %||% as.character(lazyeval::lazy(.y)$expr)
    smoothScatter(x, y, xlab=xlab, ylab=ylab, ...)
}

#' Returns a list of track names
#'
#' @description Returns a list of gpatterns tracks in the genomic database
#'
#' @param regex regular expression of track names
#' @param ignore.case see 'grep'
#' @param perl see 'grep'
#' @param fixed see 'grep'
#' @param useBytes see 'grep'
#'
#' @return An array that contains the names of tracks that match the supplied patterns
#' @export
#'
#' @examples
gpatterns.ls <- function(regex='', ignore.case = FALSE, perl=FALSE, fixed=FALSE, useBytes=FALSE){
    cov_tracks <- gtrack.ls(paste0('.*', regex, '.*\\.cov$'),
                            ignore.case = ignore.case,
                            perl = perl,
                            fixed = fixed,
                            useBytes = useBytes) %>% gsub('\\.cov', '', .)
    meth_tracks <- gtrack.ls(paste0('.*', regex, '.*\\.meth$'),
                            ignore.case = ignore.case,
                            perl = perl,
                            fixed = fixed,
                            useBytes = useBytes) %>% gsub('\\.meth', '', .)
    unmeth_tracks <- gtrack.ls(paste0('.*', regex, '.*\\.unmeth$'),
                             ignore.case = ignore.case,
                             perl = perl,
                             fixed = fixed,
                             useBytes = useBytes) %>% gsub('\\.unmeth', '', .)
    tracks <- tibble(track = cov_tracks) %>% inner_join(tibble(track = meth_tracks), by='track') %>% inner_join(tibble(track = unmeth_tracks), by='track') %>% .$track
    return(tracks)
}

# Patterns utils -------------------------------------------



#
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

.gpatterns.tcpgs_colClasses <- function(uniq=TRUE){
    if (uniq){
        colClasses <- c(
            read_id = 'character',
            chrom = 'character',
            start = 'numeric',
            end = 'character',
            strand = 'character',
            umi1 = 'character',
            umi2 = 'character',
            insert_len = 'numeric',
            num = 'numeric',
            cg_pos = 'numeric',
            meth = 'character',
            qual = 'numeric')
    } else {
        colClasses <- c(
            read_id = 'character',
            chrom = 'character',
            start = 'numeric',
            end = 'character',
            strand = 'character',
            umi1 = 'character',
            umi2 = 'character',
            insert_len = 'numeric',
            cg_pos = 'numeric',
            meth = 'character',
            qual = 'numeric')
    }
    return(colClasses)
}

.gpatterns.get_tidy_cpgs_from_dir <- function(dir,
                                              intervals = NULL,
                                              uniq = TRUE){
    if (stringr::str_sub(dir, -1) != '/') {
        dir <- paste0(dir, '/') 
    }

    files <- list.files(gsub('/$', '', dir), full.names=TRUE, pattern='tcpgs.gz')

    .get_tcpgs <- function(files){
        res <- files %>%
            map_df(function(f) fread(
                qq('gzip -d -c @{f}'),
                colClasses = .gpatterns.tcpgs_colClasses(uniq)) %>%
                    tbl_df %>%
                    mutate(cg_pos = cg_pos, meth = ifelse(meth == 'Z', 1, 0)))
        if (0 == nrow(res)){
            return(NULL)
        }
        return(res)
    }

    .intervals2files <- function(intervals, files) {
        .gpatterns.intervals2files(intervals=intervals, files=files, dir = dir)
    }


    if (!is.null(intervals)){
        tidy_intervals <- .gpatterns.get_tidy_cpgs_intervals(path=dir)        
        if (!is.character(intervals)){            
            if (all(
                unite(intervals, 'coord', chrom:end)$coord %in%
                unite(tidy_intervals, 'coord', chrom:end)$coord)
            ){
                return(.intervals2files(intervals, files) %>% .get_tcpgs)
            }
        }

        tcpgs <- tidy_intervals %>%
            gintervals.filter(intervals) %>%
            .intervals2files(files) %>%
            .get_tcpgs()
        
        f <- tcpgs %>%
            select(chrom, start=cg_pos) %>%
            mutate(end = start+1) %>%
            gintervals.neighbors1(intervals) %>%
            mutate(f = dist == 0) %>% .$f
        return(tcpgs %>% filter(f))
    }

    return(.get_tcpgs(files))
}



#' @inheritParams gpatterns::gcluster.run2
#' @export
gcluster.run3 <- partial(gcluster.run2, script = .gpatterns.sg_script)

#
.gpatterns.base_dir <- function(track)
{
    map(track, ~  c(gdir.cwd(), strsplit(.x, '.', fixed=TRUE)[[1]])) %>%
        map(~ do.call(file.path, as.list(.x))) %>%
        as_vector %>%
        return()
}


#' check if tidy_cpgs exist for track
#' @export
gpatterns.track_exists <- function(track){
    dir.exists(paste0(.gpatterns.base_dir(track), '/tidy_cpgs'))
}


.gpatterns.intervals2files <- function(intervals, files, dir){
        intervals %>% unite('coord', chrom:end, sep='_') %>%
            mutate(coord = paste0(dir, coord, '.tcpgs.gz')) %>%
            filter(coord %in% files) %>% .$coord
}

#
.gpatterns.force_chromosomes <- function(df){
    f <- df$chrom %in%  as.character(gintervals.all()$chrom )
    if (!all(f)){
        chroms <- df$chrom[!f] %>% unique %>% paste(collapse=', ')
        warning(qq("removed the following chromosomes that did not exist in genome db: @{chroms}"))
    }
    return(df %>% filter(f))
}


# Plotting functions



#' wrapper around pheatmap
#' 
#' @param pmat data frame to plot as matrix. can have an additional id column (if id_column == TRUE)
#' @param annotation data frame with 'samp' field, and additional annotations
#' @param annotation_colors tidy data frame with 3 fields: 'type' - name of the annotation, e.g. 'tissue', and pairs of 'variable' and 'color', e.g. variable='tumor', color='green'
#' @param annotation_col character vector with column names from `annotation` to use to annotate columns
#' @param annotation_row character vector with column names from `annotation` to use to annotate rows
#' @param id_column TRUE if pmat's first column contains row ids
#' @inheritParams pheatmap::pheatmap
#' 
#' @export
pheatmap1 <- function(pmat, annotation=NULL, annotation_colors=NULL, annotation_col=NULL, annotation_row = NULL, id_column = TRUE, border_color=NA, clustering_callback = function(hc, ...){dendsort::dendsort(hc)}, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", min_vals_row=1, min_vals_col=1, ...){
    pmat <- as.data.frame(pmat)
    id_colname <- 'samp'
    if (id_column){
        id_colname <- colnames(pmat)[1]
        ids <- pmat[, 1]
        pmat <- pmat[, -1]
        rownames(pmat) <- ids
    }
    
    if (!is.null(annotation)){                
        annots <- annotation %>%
            filter(samp %in% colnames(pmat) | samp %in% rownames(pmat)) %>%
            mutate_if(is.character, factor) %>%
            as.data.frame
        rownames(annots) <- annots$samp
        annots <- annots %>% select(-samp)

        if (!is.null(annotation_colors)){
            annot_cols <- plyr::dlply(annotation_colors,
                                      .(type),
                                      function(x) x %>%
                                        select(-type) %>%
                                        spread(variable, color) %>%
                                         unlist)  
            annot_cols <- annot_cols[colnames(annots)]         
        } else {
            annot_cols <- NA
        }
    } else {
        annots <- NA
    }

    if (!is.null(annotation_col)){
        annotation_col <- annots[, annotation_col]
    } else {
        annotation_col <- NA
    }

    if (!is.null(annotation_row)){
        annotation_row <- annots[, annotation_row]
    } else {
        annotation_row <- NA
    }

    if (min_vals_row > 1){
        pmat <- pmat[apply(pmat, 1, function(x) sum(!is.na(x)) >= 100), ]
    }

    if (min_vals_col > 1){
        pmat <- pmat[, apply(pmat, 2, function(x) sum(!is.na(x)) >= 100)]   
    }

    pheatmap::pheatmap(pmat,
                            border_color=border_color,
                            annotation_col=annotation_col,                            
                            annotation_colors=annot_cols,
                            annotation_row=annotation_row,   
                            clustering_callback=clustering_callback, 
                            clustering_distance_rows=clustering_distance_rows,
                            clustering_distance_cols=clustering_distance_cols,                
                             ...)    
}



#' Build a gradient pallete based on colors at specific values.
#' 
#' @param shades dataframe with columns:
#'     - point - the values at which the pure colors should be set
#'     - color - the pure colors 
#' @param length number of colors
#' 
#' @export
build_pallette <- function(shades, length)

{

    from <- min(shades$point)

    to <- max(shades$point)

    values <- seq(from, to, length.out=length)

    

    colors <- shades$color[order(shades$point)]

    points <- sapply(sort(shades$point), function(p) which.min(abs(values-p)))

    lens <- points[-1] - points[-length(points)] + 1



    palette <- lapply(1:(length(points)-1), function(i) colorRampPalette(c(colors[i], colors[i+1]), space="Lab")(lens[i]))

    palette[-1] <- lapply(palette[-1], function(p) p[-1])

    return (do.call(c, palette))

}





.check_tracks_exist <- function(tracks, suffixes = c("")) {
    for (suffix in suffixes) {
        trs <- paste0(tracks, ".", suffix)
        if (any(!gtrack.exists(trs))) {
            message("The following tracks do not exist:")
            print(trs[!gtrack.exists(trs)])
            stop()
        }
    }
}

#' @export
do.call_ellipsis <- function(f, additional_params=list(), ...){
    f_args <- names(as.list(args(f)))
    elipsis <- list(...)        
    if (!is.null(names(elipsis))){
        new_elipsis <- list()        
        for (x in names(elipsis)[names(elipsis) %in% f_args]) { 
            new_elipsis[[x]] <- elipsis[[x]] 
        }                
        do.call(f, c(additional_params, new_elipsis))
    } else {
        do.call(f, additional_params)
    }
}

.hclust_order <- function(d, keys, variable, value, tidy=TRUE, ...){
    if (tidy){
        d_mat <- d %>%
        select(one_of(c(keys, variable, value))) %>%
        spread_(variable, value)
    } else {
        d_mat <- d
    }

    # remove empty columns
    d_mat <- d_mat[, apply(d_mat, 2, function(x) any(!is.na(x)))]

    hc <- d_mat  %>%
        select(-one_of(keys)) %>%
        as.matrix %>%
        dist %>%
        hclust(...)
    d <- d_mat %>%
        select(one_of(keys)) %>%
        mutate(ord=hc$ord) %>%
        right_join(d, by=keys) %>%
        arrange(ord)
    return(d)
}


# # Misha overrides ------------------------------------------------

# #' overrides gtrack.create_sparse without rescaning database
# #' @export
# gtrack.create_sparse <- function (track = NULL, description = NULL, intervals = NULL,
#            values = NULL, rescan = TRUE)
# {
#     if (is.null(substitute(track)) || is.null(description) ||
#         is.null(intervals) || is.null(values))
#         stop("Usage: gtrack.create_sparse(track, description, intervals, values)",
#              call. = F)
#     .gcheckroot()
#     trackstr <- do.call(.gexpr2str, list(substitute(track)),
#                         envir = parent.frame())
#     intervalsstr <- deparse(substitute(intervals), width.cutoff = 500)[1]
#     valuesstr <- deparse(substitute(values), width.cutoff = 500)[1]
#     trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.",
#                                                            "/", trackstr), sep = "/"))
#     direxisted <- file.exists(trackdir)
#     if (!is.na(match(trackstr, get("GTRACKS"))))
#         stop(sprintf("Track %s already exists", trackstr), call. = F)
#     .gconfirmtrackcreate(trackstr)
#     success <- FALSE
#     tryCatch({
#         .gcall("gtrack_create_sparse", trackstr, intervals, values,
#                new.env(parent = parent.frame()), silent = TRUE)
#         gdb.reload(rescan = rescan)
#         .gtrack.attr.set(trackstr, "created.by", sprintf("gtrack.create_sparse(%s, description, %s, %s)",
#                                                          trackstr, intervalsstr, valuesstr), T)
#         .gtrack.attr.set(trackstr, "created.date", date(), T)
#         .gtrack.attr.set(trackstr, "description", description,
#                          T)
#         success <- TRUE
#     }, finally = {
#         if (!success && !direxisted) {
#             unlink(trackdir, recursive = TRUE)
#             gdb.reload()
#         }
#     })
#     retv <- 0
# }

# #' @export
# gtrack.rm <- function (track = NULL, force = FALSE, rescan = TRUE)
# {
#     if (is.null(substitute(track)))
#         stop("Usage: gtrack.rm(track, force = FALSE)", call. = F)
#     .gcheckroot()
#     trackname <- do.call(.gexpr2str, list(substitute(track)),
#                          envir = parent.frame())
#     if (is.na(match(trackname, get("GTRACKS")))) {
#         if (force)
#             return(invisible())
#         stop(sprintf("Track %s does not exist", trackname), call. = F)
#     }
#     answer <- "N"
#     if (force)
#         answer <- "Y"
#     else {
#         str <- sprintf("Are you sure you want to delete track %s (Y/N)? ",
#                        trackname)
#         cat(str)
#         answer <- toupper(readLines(n = 1))
#     }
#     if (answer == "Y" || answer == "YES") {
#         dirname <- sprintf("%s.track", paste(get("GWD"), gsub("\\.",
#                                                               "/", trackname), sep = "/"))
#         unlink(dirname, recursive = TRUE)
#         if (file.exists(dirname))
#             cat(sprintf("Failed to delete track %s\n", trackname))
#         else gdb.reload(rescan = rescan)
#     }
# }

# gtrack.array.import <- function (track = NULL, description = NULL, ...)
# {
#     args <- as.list(substitute(list(...)))[-1L]
#     if (is.null(substitute(track)) || is.null(description) ||
#         !length(args))
#         stop("Usage: gtrack.array.import(track, description, [src]+)",
#              call. = F)
#     .gcheckroot()
#     trackstr <- do.call(.gexpr2str, list(substitute(track)),
#                         envir = parent.frame())
#     srcs <- c()
#     colnames <- list()
#     for (src in args) {
#         src <- do.call(.gexpr2str, list(src), envir = parent.frame())
#         srcs <- c(srcs, src)
#         if (is.na(match(src, get("GTRACKS"))))
#             colnames[[length(colnames) + 1]] <- as.character(NULL)
#         else {
#             if (.gcall_noninteractive(gtrack.info, src)$type !=
#                 "array")
#                 stop(sprintf("Track %s: only array tracks can be used as a source",
#                              src), call. = F)
#             colnames[[length(colnames) + 1]] <- names(.gtrack.array.get_colnames(src))
#         }
#     }
#     trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.",
#                                                            "/", trackstr), sep = "/"))
#     direxisted <- file.exists(trackdir)
#     if (!is.na(match(trackstr, get("GTRACKS"))))
#         stop(sprintf("Track %s already exists", trackstr), call. = F)
#     .gconfirmtrackcreate(trackstr)
#     success <- FALSE
#     tryCatch({
#         colnames <- .gcall("garrays_import", trackstr, srcs,
#                            colnames, new.env(parent = parent.frame()), silent = TRUE)
#         gdb.reload()
#         .gtrack.array.set_colnames(trackstr, colnames, FALSE)
#         created.by <- sprintf("gtrack.array.import(\"%s\", description, src = c(\"%s\"))",
#                               trackstr, paste(srcs, collapse = "\", \""))
#         .gtrack.attr.set(trackstr, "created.by", created.by,
#                          T)
#         .gtrack.attr.set(trackstr, "created.date", date(), T)
#         .gtrack.attr.set(trackstr, "description", description,
#                          T)
#         success <- TRUE
#     }, finally = {
#         if (!success && !direxisted) {
#             unlink(trackdir, recursive = TRUE)
#             gdb.reload(rescan=FALSE)
#             warning('Please run gdb.reload() once finished uploading array tracks')
#         }
#     })
#     retv <- 0
# }

# gtrack.create <- function (track = NULL, description = NULL, expr = NULL, iterator = NULL,
#           band = NULL, rescan=TRUE){
#     if (is.null(substitute(track)) || is.null(description) ||
#         is.null(substitute(expr)))
#         stop("Usage: gtrack.create(track, description, expr, iterator = NULL, band = NULL)",
#              call. = F)
#     .gcheckroot()
#     trackstr <- do.call(.gexpr2str, list(substitute(track)),
#                         envir = parent.frame())
#     exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
#     .iterator <- do.call(.giterator, list(substitute(iterator)),
#                          envir = parent.frame())
#     trackdir <- sprintf("%s.track", paste(get("GWD"), gsub("\\.",
#                                                            "/", trackstr), sep = "/"))
#     direxisted <- file.exists(trackdir)
#     if (!is.na(match(trackstr, get("GTRACKS"))))
#         stop(sprintf("Track %s already exists", trackstr), call. = F)
#     .gconfirmtrackcreate(trackstr)
#     success <- FALSE
#     tryCatch({
#         if (.ggetOption("gmultitasking"))
#             .gcall("gtrackcreate_multitask", trackstr, exprstr,
#                    .iterator, band, new.env(parent = parent.frame()),
#                    silent = TRUE)
#         else .gcall("gtrackcreate", trackstr, exprstr, .iterator,
#                     band, new.env(parent = parent.frame()), silent = TRUE)
#         gdb.reload(rescan = rescan)
#         .gtrack.attr.set(trackstr, "created.by", sprintf("gtrack.create(%s, description, %s, iterator=%s)",
#                                                          trackstr, exprstr, deparse(substitute(iterator),
#                                                                                     width.cutoff = 500)[1]), T)
#         .gtrack.attr.set(trackstr, "created.date", date(), T)
#         .gtrack.attr.set(trackstr, "description", description,
#                          T)
#         success <- TRUE
#     }, finally = {
#         if (!success && !direxisted) {
#             unlink(trackdir, recursive = TRUE)
#             gdb.reload()
#         }
#     })
#     retv <- 0
# }






# .gintervals.apply <- function (chroms, intervals, intervals.set.out, FUN, ...)
# {
#     if (!is.null(intervals.set.out))
#         fullpath <- .gintervals.check_new_set(intervals.set.out)
#     if (is.data.frame(intervals))
#         intervals <- list(intervals)
#     chroms$size <- NULL
#     if ("chrom" %in% colnames(chroms))
#         chroms <- data.frame(chrom = chroms[with(chroms, order(chrom)),
#                                             ])
#     else chroms <- chroms[with(chroms, order(chrom1, chrom2)),
#                           ]
#     if (any(unlist(lapply(intervals, function(intervals) {
#         .gintervals.is_bigset(intervals) || .gintervals.needs_bigset(intervals)
#     })))) {
#         stats <- NULL
#         zeroline <- NULL
#         success <- FALSE
#         t <- Sys.time()
#         progress.percentage <- -1
#         tryCatch({
#             if (!is.null(intervals.set.out))
#                 dir.create(fullpath, recursive = T, mode = "0777")
#             if (.gintervals.is1d(intervals[[1]])) {
#                 mapply(function(chrom) {
#                     loaded_intervals <- lapply(intervals, function(intervals) {
#                         .gintervals.load_ext(intervals, chrom = chrom)
#                     })
#                     res <- do.call(FUN, list(loaded_intervals,
#                                              ...))
#                     if (!is.null(intervals.set.out) && !is.null(res) &&
#                         nrow(res) > 0) {
#                         zeroline <<- res[0, ]
#                         .gintervals.big.save(fullpath, res, chrom = chrom)
#                         stat <- .gcall("gintervals_stats", res, new.env(parent = parent.frame()))
#                         stats <<- rbind(stats, data.frame(chrom = chrom,
#                                                           stat))
#                     }
#                     if (as.integer(difftime(Sys.time(), t, units = "secs")) >
#                         3) {
#                         t <<- Sys.time()
#                         percentage <- as.integer(100 * match(chrom,
#                                                              chroms$chrom)/nrow(chroms))
#                         if (percentage < 100 && progress.percentage !=
#                             percentage) {
#                             cat(sprintf("%d%%...", percentage))
#                             progress.percentage <<- percentage
#                         }
#                     }
#                 }, chroms$chrom)
#             }
#             else {
#                 mapply(function(chrom1, chrom2) {
#                     loaded_intervals <- lapply(intervals, function(intervals) {
#                         .gintervals.load_ext(intervals, chrom1 = chrom1,
#                                              chrom2 = chrom2)
#                     })
#                     res <- do.call(FUN, list(loaded_intervals,
#                                              ...))
#                     if (!is.null(intervals.set.out) && !is.null(res) &&
#                         nrow(res) > 0) {
#                         zeroline <<- res[0, ]
#                         .gintervals.big.save(fullpath, res, chrom1 = chrom1,
#                                              chrom2 = chrom2)
#                         stat <- .gcall("gintervals_stats", res, new.env(parent = parent.frame()))
#                         stats <<- rbind(stats, data.frame(chrom1 = chrom1,
#                                                           chrom2 = chrom2, stat))
#                     }
#                     if (as.integer(difftime(Sys.time(), t, units = "secs")) >
#                         3) {
#                         t <<- Sys.time()
#                         percentage <- as.integer(100 * which(chroms$chrom1 ==
#                                                                  chrom1 & chroms$chrom2 == chrom2)/nrow(chroms))
#                         if (percentage < 100 && progress.percentage !=
#                             percentage) {
#                             cat(sprintf("%d%%...", percentage))
#                             progress.percentage <<- percentage
#                         }
#                     }
#                 }, chroms$chrom1, chroms$chrom2)
#             }
#             if (!is.null(intervals.set.out)) {
#                 if (is.null(stats))
#                     return(retv <- NULL)
#                 .gintervals.big.save_meta(fullpath, stats, zeroline)
#             }
#             if (progress.percentage >= 0)
#                 cat("100%\n")
#             success <- TRUE
#             if (!is.null(intervals.set.out) && !.gintervals.needs_bigset(intervals.set.out))
#                 .gintervals.big2small(intervals.set.out)
#         }, finally = {
#             if (!success && !is.null(intervals.set.out))
#                 unlink(fullpath, recursive = TRUE)
#         })
#     }
#     else {
#         loaded_intervals <- lapply(intervals, .gintervals.load_ext)
#         res <- do.call(FUN, list(loaded_intervals, ...))
#         if (!is.null(intervals.set.out) && !is.null(res) && nrow(res) >
#             0) {
#             if (.gintervals.is1d(res))
#                 res <- res[order(res$chrom), ]
#             else res <- res[order(res$chrom1, res$chrom2), ]
#             if (.gintervals.needs_bigset(res))
#                 .gintervals.small2big(intervals.set.out, res)
#             else .gintervals.save_file(fullpath, res)
#         }
#         else return(NULL)
#     }
#     if (!is.null(intervals.set.out))
#         gdb.reload(rescan = FALSE)
# }
