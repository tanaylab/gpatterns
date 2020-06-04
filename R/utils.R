# Utility functions ------------------------------------------------

#' @export
qq <- GetoptLong::qq

#' @export
qqv <- function(text) {GetoptLong::qq(text, envir=parent.frame(), collapse=FALSE)}

#' @export
comify <- scales::comma

#' @export
red_message <- function(msg, log=TRUE){  
    msg <- red(glue(msg, .envir = parent.frame(1)))
    if (log){
        loginfo(msg)    
    }
    message(msg)
}

#' @export
blue_message <- function(msg, log=TRUE){ 
    msg <- blue(glue(msg, .envir = parent.frame(1)))
    if (log){
        loginfo(msg)    
    }
    message(msg)
}


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

#' Remove gpatterns track
#' @param track track name
#' @export
gpatterns.rm <- function(track){
    for (suffix in c('meth', 'unmeth', 'cov', 'avg')){
        gtrack.rm(qq('@{track}.@{suffix}'))
    }
    system('rm -rf @{.gpatterns.base_dir(track)}')
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
            qual = 'numeric')#,
    #         h = 'numeric', 
    #         H = 'numeric', 
    #         x = 'numeric', 
    #         X = 'numeric')
    }
    return(colClasses)
}

#' @export
.gpatterns.get_tidy_cpgs_from_dir <- function(dir,
                                              intervals = NULL,
                                              uniq = TRUE){
    if (stringr::str_sub(dir, -1) != '/') {
        dir <- paste0(dir, '/') 
    }

    files <- list.files(gsub('/$', '', dir), full.names=TRUE, pattern='tcpgs.gz')

    .get_tcpgs <- function(files){
        classes <- .gpatterns.tcpgs_colClasses(uniq)
        res <- files %>%
            map_df(function(f) fread(
                qq("gzip -d -c @{f} | cut -d',' -f1-@{length(classes)}"),
                colClasses = classes) %>%
                    as_tibble() %>%
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


#' @export
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
pheatmap1 <- function(pmat, annotation=NULL, annotation_colors=NULL, annotation_col=NULL, annotation_row = NULL, id_column = TRUE, border_color=NA, clustering_callback = function(hc, ...){dendsort::dendsort(hc)}, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", min_vals_row=1, min_vals_col=1, na_col='gray', ...){
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

        if (!is.null(annotation_colors)){
            annot_cols <- plyr::dlply(annotation_colors,
                                      .(type),
                                      function(x){
                                        var_ord <- x$variable
                                        x %>%
                                            select(-type) %>%
                                            spread(variable, color) %>%
                                            select(one_of(var_ord)) %>% 
                                            unlist
                                      }) 
            annot_cols <- annot_cols[colnames(annots)]         
        } else {
            annot_cols <- NA
        }
    } else {
        annots <- NA
    }
    
    if (min_vals_row > 1){
        pmat <- pmat[apply(pmat, 1, function(x) sum(!is.na(x)) >= 100), ]
    }

    if (min_vals_col > 1){
        pmat <- pmat[, apply(pmat, 2, function(x) sum(!is.na(x)) >= 100)]   
    }

    if (!is.null(annotation_col)){
        annotation_col <- annots %>% filter(samp %in% colnames(pmat)) %>% select(one_of(c(annotation_col, 'samp')))
        rownames(annotation_col) <- annotation_col$samp
        annotation_col <- annotation_col %>% select(-samp)        
    } else {
        annotation_col <- NA
    }

    if (!is.null(annotation_row)){
        annotation_row <- annots %>% filter(samp %in% rownames(pmat)) %>% select(one_of(c(annotation_row, 'samp')))
        rownames(annotation_row) <- annotation_row$samp
        annotation_row <- annotation_row %>% select(-samp)        
    } else {
        annotation_row <- NA
    }
    
    pheatmap::pheatmap(pmat,
                            border_color=border_color,
                            annotation_col=annotation_col,                            
                            annotation_colors=annot_cols,
                            annotation_row=annotation_row,   
                            clustering_callback=clustering_callback, 
                            clustering_distance_rows=clustering_distance_rows,
                            clustering_distance_cols=clustering_distance_cols, 
                            na_col=na_col,               
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

#' @export
do.call_tibble <- function(f, tbl, additional_params = list(), additional_funcs=list(), use_glue=TRUE){   

    if (length(additional_funcs) == 0){
        funcs <- c(f)
    } else {
        funcs <- c(f, additional_funcs)
    }
    func_args <- map(funcs, ~ names(as.list(args(.)))) %>% flatten_chr() %>% unique()

    if (use_glue){      
        for (arg in func_args){
            if (has_name(tbl, arg) && is.character(tbl[[arg]]) && length(tbl[[arg]]) == 1){
                tbl[[arg]] <- tbl %>% glue::glue_data(tbl[[arg]])
            }
        }        
    }

    suppressWarnings(tbl <- tbl %>% select(one_of(func_args)))

    tbl <- as.list(tbl) %>% modify_if(is.tibble, as.list) %>% modify_if(is.list, ~ .x[[1]])
    tbl <- modify_if(tbl, is.tibble, ~ as.list(.x)[[1]])

    if (length(additional_params) > 0){
        additional_params <- additional_params[!(names(additional_params) %in% names(tbl))]
    
        additional_params <- additional_params[names(additional_params) %in% func_args]
    }        
    
    if (!is.null(names(tbl))){
        do.call(f, c(additional_params, tbl))
    } else {
        do.call(f, additional_params)
    }
}


do.call_vec <- function(func, args_vec){
    func_args <- names(args_vec)[names(args_vec) %in% names(as.list(args(func)))]   
    do.call(func, args_vec[func_args])
}


# build_command <- function(args_vec, funcs, func_name, extra_args=''){
#     func_args <- map(funcs, ~ names(as.list(args(.)))) %>% flatten_chr() %>% unique()   
#     func_args <- names(args_vec)[names(args_vec) %in% func_args]    
#     args_vec <-  modify_if(args_vec, ~ is.tibble(.x) && nrow(.x) == 1, ~ .x[[1]])    
#     args_vec <-  modify_if(args_vec, ~ is.character(.x) && length(.x) == 1, ~ sprintf('"%s"', .x))    
#     args_vec <-  modify_if(args_vec, ~ length(.x) > 1, ~ as.character(as_tibble(.x)))          
#     sprintf("%s(%s%s)", func_name, extra_args, map2_chr(names(args_vec[func_args]), args_vec[func_args], ~ glue('{.x} = {.y}')) %>% paste(collapse=', ') )
# }

# add_function_args <- function(config, args, funcs){ 
#     func_args <- map(funcs, ~ names(as.list(args(.)))) %>% flatten_chr() %>% unique()    
#     args <- args[names(args) %in% func_args]    
#     for (arg in names(args)){
#         if(!has_name(config, arg)){
#             if (length(args[[arg]]) > 1){                
#                 config[[arg]] <- rep(list(args[[arg]]), nrow(config))
#             } else {
#                 if (is.character(args[[arg]])){
#                     config[[arg]] <- glue::glue_data(config, args[[arg]])       
#                 } else {
#                     config[[arg]] <- args[[arg]]   
#                 }                
#             }
#         }
#     }   
#     return(config)
# }

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

get_qual_colors <- function(n){
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    return(col_vector[1:n])
}

untable <- function(df, column){
    df <- as_tibble(df)    
    return(df[rep(1:nrow(df), df[[column]]), ])
}