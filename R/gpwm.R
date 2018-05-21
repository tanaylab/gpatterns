########################################################################
#' @export
gpwm.extract_motif <- function(track, intervals, colname=NULL, use_cache=TRUE){
    gmax.processes <- getOption('gmax.processes')
    options(gmax.processes=1)

    if (is.null(colname)) {
        colname = track
    }

    intervals$chrom <- factor(intervals$chrom, levels=levels(ALLGENOME[[1]]$chrom))

    full_intervals <- intervals
    if (use_cache) {
        intervals <- .gpwm_cache_misses(track, intervals)
    }
    intervals <- intervals %>%
                 select(chrom, start, end) %>%
                 distinct()

    if (nrow(intervals) > 0) {
        vmotif <- qq('@{track}__vtrack_max__')
        gvtrack.create(vmotif, track, 'max')    
        intervals <- gextract(vmotif, intervals=intervals, iterator=intervals, colnames='gpwm_value__') 
        gvtrack.rm(vmotif)

        intervals <- intervals %>% select(-intervalID)
    }

    options(gmax.processes=gmax.processes)
    if (!use_cache) {
        return(full_intervals %>% 
               left_join(intervals, by=c('chrom', 'start', 'end')) %>%
               rename_(.dots=setNames('gpwm_value__', colname)))
    }

    if (nrow(intervals) > 0) {
        .gpwm_cache__[[track]] <- bind_rows(.gpwm_cache__[[track]], intervals) %>%
                                  arrange(chrom, start, end)
    }

    return(full_intervals %>% 
           left_join(.gpwm_cache__[[track]], by=c('chrom', 'start', 'end')) %>%
           rename_(.dots=setNames('gpwm_value__', colname)))
}


########################################################################
#' @export
gpwm.extract <- function(..., intervals=NULL, colnames=NULL, tidy=FALSE, 
                           parallel=NULL, 
                           opt.flags = "",
                           max.jobs = 400,                           
                           R = paste0(R.home(component='bin'), '/R'),
                           packages = NULL,
                           jobs_title = NULL,
                           job_names = NULL,                           
                           queue = NULL,
                           memory = NULL,
                           threads = NULL,
                           io_saturation = NULL,
                           queue_flag = '-q @{queue}',
                           memory_flag = '-l mem_free=@{memory}G',
                           threads_flag = '-pe threads @{threads}',
                           io_saturation_flag = '-l io_saturation=@{io_saturation}',
                           script =paste(Sys.getenv('ANALYSIS_HOME'), 'common', 'sgjob.sh', sep='/')){
    arg_intervals <- NULL
    arg_colnames <- NULL
    arg_parallel <- NULL

    tracks <- list(...)
    point <- which(!sapply(tracks, is.character))[1]
    if (point == length(tracks)) {
        arg_intervals <- tracks[[length(tracks)]]
        tracks <- tracks[1:(length(tracks)-1)]
    }
    else if (point == length(tracks)-1) {
        arg_intervals <- tracks[[length(tracks) - 1]]
        arg_colnames <-  tracks[[length(tracks)]]
        tracks <- tracks[1:(length(tracks)-2)]
    }
    else if (point == length(tracks)-2) {
        arg_intervals <- tracks[[length(tracks) - 2]]
        arg_colnames <- tracks[[length(tracks) - 1]]
        arg_parallel <- tracks[[length(tracks)]]
        tracks <- tracks[1:(length(tracks)-3)]        
    }
    else {
        stop('Too many positional arguments.')
    }
    
    if (!is.null(arg_intervals)) {
        if (!is.null(intervals)) {
            stop('formal argument "intervals" matched by multiple actual arguments')
        }
        intervals <- arg_intervals
    }
    if (!is.null(arg_colnames)) {
        if (!is.null(colnames)) {
            stop('formal argument "colnames" matched by multiple actual arguments')
        }
        colnames <- arg_colnames
    }
    if (!is.null(arg_parallel)) {
        if (!is.null(parallel)) {
            stop('formal argument "parallel" matched by multiple actual arguments')            
        }
        parallel <- arg_parallel;
    }
    
    tracks <- do.call(c, tracks)
        
    if (is.null(intervals)) {
        stop('argument "intervals" is missing')
    }
    
    if (is.null(colnames)) {
        colnames <- tracks
    }
    if (!is.character(colnames)) {
        stop('Invalid value for argument "colnames"')
    }
    if (length(tracks) != length(colnames)) {
        stop('Length of "tracks" and "colnames" do not match')
    }
    names(colnames) <- tracks
    
    if (is.null(parallel)) {
        parallel = TRUE
    }
    
    intervals$chrom <- factor(intervals$chrom, levels=levels(ALLGENOME[[1]]$chrom))    
    full_intervals <- intervals
    
    intervals <- intervals %>%
                 select(chrom, start, end) %>%
                 distinct()
    
    
    intervals <- lapply(tracks, function(track) {.gpwm_cache_misses(track, intervals) %>% select(chrom, start, end)})
    names(intervals) <- tracks
    intervals <- intervals[sapply(intervals, nrow) > 0]
    
    if (length(intervals) > 0) {
        commands <- sapply(1:length(intervals),
                           function(i) {qq("gpwm.extract_motif(names(intervals)[@{i}], intervals[[@{i}]], 'gpwm_value__', use_cache=FALSE)")})       

        if (!parallel) {
            parsed_cmds <- parse(text=commands)
            names(parsed_cmds) <- names(intervals)
            for (track in names(intervals)) {
                result <- eval(parsed_cmds[track])
                .gpwm_cache__[[track]] <- bind_rows(.gpwm_cache__[[track]], result) %>%
                                          arrange(chrom, start, end)
            }
        }
        else {
          
            results <- gcluster.run2(command_list=commands, max.jobs=max.jobs, R=R, packages=packages, jobs_title=jobs_title, job_names=job_names, collapse_results=F, memory=memory, threads=threads, io_saturation=io_saturation, memory_flag=memory_flag, threads_flag=threads_flag, io_saturation_flag=io_saturation_flag, script=script, queue=queue, queue_flag=queue_flag)            
            
            names(results) <- names(intervals)
            
            .gpwm.verify_jobs_successes(results, commands)
            
            for (track in names(intervals)) {
                .gpwm_cache__[[track]] <- bind_rows(.gpwm_cache__[[track]], results[[track]]$retv) %>%
                                          arrange(chrom, start, end)
            }
        }
    }
    
    for (track in tracks) {
        full_intervals <- full_intervals %>%
                          left_join(.gpwm_cache__[[track]], by=c('chrom', 'start', 'end')) %>%
                          rename_(.dots=setNames('gpwm_value__', colnames[track]))
    }
    
    if (tidy){        
        full_intervals <- full_intervals %>% gather('track', 'val', -(chrom:end))
    }
    return(full_intervals)
}



########################################################################
#' @export
gpwm.extract_all <- function(pattern, intervals, colname_prefix=NULL, tidy=FALSE, parallel=TRUE, ...){
    tracks <- gtrack.ls(qq('^@{pattern}'), perl=TRUE)
    if (length(tracks) == 0) {
        return(intervals)
    }
    
    colnames <- tracks
    if (!is.null(colname_prefix)) {
        colnames <- paste0(colname_prefix, colnames)
    }
    
    res <- gpwm.extract(tracks, intervals, colnames, tidy=tidy, parallel=parallel, ...)
    if (tidy){
        res <- res %>% mutate(track = gsub(qq('^@{pattern}\\.?'), '', track))
    }
    return(res)
}



########################################################################
gpwm.clear_cache <- function(){
    for (name in ls(envir=.gpwm_cache__, all.names=FALSE)) {
        rm(name, envir=.gpwm_cache__)
    }
}


########################################################################
gpwm.load_cache <- function(filename=NULL){
    if (is.null(filename)) {
        filename <- .gpwm_cache__$.filename__
    }
    if (is.null(filename)) {
        return(xxx <- NULL)
    }
    
    if (!file.exists(filename)) {
        .gpwm_cache__$.filename__ <- filename
        return(xxx <- NULL)
    }
    
    current <- .gpwm_cache__
    .gpwm_cache__ <<- new.env()
    
    load(filename, envir=.gpwm_cache__)
    .gpwm_cache__$.filename__ <- filename

    for (name in ls(envir=current, all.names=TRUE)) {
        if (is.null(.gpwm_cache__[[name]])) {
            .gpwm_cache__[[name]] <- current[[name]]
            next
        }

        if (is.data.frame(.gpwm_cache__[[name]])) {
            .gpwm_cache__[[name]] <- bind_rows(.gpwm_cache__[[name]],
                                               current[[name]] %>% anti_join(.gpwm_cache__[[name]], by=c('chrom', 'start', 'end'))) %>%
                                     arrange(chrom, start, end)
        }
    }
}


########################################################################
gpwm.save_cache <- function(filename=NULL){
    if (is.null(filename)) {
        filename = .gpwm_cache__$.filename__
    }
    else {
        .gpwm_cache__$filename__ = filename
    }
    
    if (is.null(filename)) {
        stop("No cache filename is available.")
    }

    save(list=ls(envir=.gpwm_cache__), file=filename, envir=.gpwm_cache__)
}


########################################################################
.gpwm_cache_misses <- function(track, intervals){
    if (!is.null(.gpwm_cache__[[track]]) && (nrow(.gpwm_cache__[[track]]) > 0)) {
        return(intervals %>% 
            anti_join(.gpwm_cache__[[track]], by=c('chrom', 'start', 'end')) %>% 
            arrange(chrom, start, end))
    }
    
    if (is.null(.gpwm_cache__[[track]]) && gtrack.exists(track)) {
        .gpwm_cache__[[track]] <- data.frame()
    }

    return(intervals)
}

########################################################################
.gpwm.verify_jobs_successes <- function(results, commands){
    mask <- sapply(results, function(x) {x$exit.status != 'success'})
    if (any(mask)) {
        failed <- paste0("'", commands[mask], "': Bad exit.status", collapse='\n')
        stop(qq("gcluster.run() failed running the following commands:\n@{failed}"))
    }

    mask <- sapply(results, function(x) {class(x$retv) == 'try-error'}) 
    if (any(mask)) {
        failed <- paste0("'", commands[mask], ": ", results$retv[mask], collapse='\n')
        stop(qq("gcluster.run() failed running the following commands:\n@{failed}"))
    }
}



########################################################################
#' @export
gpwm.max_val_quantile <- function(track, size, quantiles = c(0.9, 0.95, 0.99, 
    0.995, 0.999)) {    
    gmax.processes <- getOption('gmax.processes')
    options(gmax.processes=1)    
    vtrack <- paste(track, "max", sep = "_")
    gvtrack.create(vtrack, track, "max")
    q <- gquantiles(vtrack, percentiles = quantiles, intervals = ALLGENOME[[1]], 
        iterator = size)
    options(gmax.processes=gmax.processes)
    return(q)
}

########################################################################
#' @export
gpwm.max_val_quantile_all <- function(pattern, size, quantiles = c(0.9, 
    0.95, 0.99, 0.995, 0.999),  parallel=TRUE, ...) {    
    motif_tracks <- gtrack.ls(qq('^@{pattern}'), perl=TRUE)
    quantiles_str <- paste("c(", paste(quantiles, collapse = ", "), ")")
        
    commands <- qqv('gpwm.max_val_quantile("@{motif_tracks}", @{size}, @{quantiles_str})')
    if (!parallel) {
        parsed_cmds <- parse(text=commands)
        names(parsed_cmds) <- motif_tracks    
        res <- map(motif_tracks, ~ eval(parsed_cmds[.x]))     
        res <-  map2_df(res, motif_tracks, ~ tibble(track = .y, quant=quantiles, value=.x$retv)) 
    } else {
        res <- gcluster.run2(command_list=commands, ...)                
        .gpwm.verify_jobs_successes(res, commands)
        res <-  map2_df(res, motif_tracks, ~ tibble(track = .y, quant=quantiles, value=.x$retv))      
    }       
    
    res <- res %>% 
        mutate(track = gsub(qq('^@{pattern}\\.?'), '', track),
               size = size)
    
    return(res)
}

########################################################################
#' @export
gpwm.get_global_quantiles <- function(pattern, size, quantiles=c(0.9, 
    0.95, 0.99, 0.995, 0.999), ...){
    global_quantiles_fn <- qq('@{.gpwm.base_dir(pattern)}/motif_max_val_quant_@{size}.csv')

    if (file.exists(global_quantiles_fn)){
        res <- read_csv(global_quantiles_fn, col_types=cols(track = col_character(),
                                                            quant = col_double(),
                                                            value = col_double(),
                                                            size = col_integer()))
    } else {
        res <- gpwm.max_val_quantile_all(pattern, size, quantiles, ...)
        write_csv(res, global_quantiles_fn)
    }
    return(res)
}

########################################################################
#' @export
gpwm.add_global_quantiles <- function(motif_intervals, global_quantiles=NULL, pattern=NULL, size=NULL, quantile_thresh = 0.99, ...){
    if (is.null(global_quantiles)){
        global_quantiles <- gpwm.get_global_quantiles(pattern=pattern, size=size, ...)    
    }
    global_quantiles <- global_quantiles %>% filter(quant == quantile_thresh) %>% rename(glob_val = value)
    motif_intervals %>% left_join(global_quantiles, by='track')    
}

########################################################################
#' @export
gpwm.motif_enrich <- function(fg, bg, global_quantiles=NULL, pattern=NULL, size=NULL, quantile_thresh=0.99, min_n_fg=4, min_n_bg=5){
    fg <- gpwm.add_global_quantiles(fg, global_quantiles=global_quantiles, pattern=pattern, size=size)
    bg <- gpwm.add_global_quantiles(bg, global_quantiles=global_quantiles, pattern=pattern, size=size)    

    fg_num <- fg %>% group_by(track) %>% summarise(n_fg = n(), n_fg_ok = sum(val >= glob_val, na.rm=TRUE))
    bg_num <- bg %>% group_by(track) %>% summarise(n_bg = n(), n_bg_ok = sum(val >= glob_val, na.rm=TRUE))

    counts <- fg_num %>% left_join(bg_num, by='track')    
    counts <- counts %>% purrrlyr::by_row(~ phyper(.x$n_fg_ok - 1, .x$n_bg_ok, .x$n_bg - .x$n_bg_ok, .x$n_fg, lower.tail=FALSE), .to='pval' ) %>% unnest(pval)

    counts <- counts %>% mutate(qval =p.adjust(pval))
    counts <- counts %>% mutate(rel_enrich = (n_fg_ok / n_fg) / (n_bg_ok / n_bg), 
                                abs_enrich = (n_fg_ok / n_fg) / (1 - quantile_thresh))
    return(counts)    
}

########################################################################
#' @export
gpwm.motif_enrich_per_cluster <- function(fg, bg, global_quantiles=NULL, pattern=NULL, size=NULL, parallel=TRUE, ...){
    if (is.null(global_quantiles)){
        global_quantiles <- gpwm.get_global_quantiles(pattern=pattern, size=size) 
    }    
    fg %>% plyr::ddply(plyr::.(clust), function(x) gpwm.motif_enrich(x, bg, global_quantiles=global_quantiles, ...), .parallel=parallel, .progress = 'text') %>% mutate(qval =p.adjust(pval)) %>% as_tibble()    
}

########################################################################
#' @export
gpwm.plot_cluster_motif_enrich <- function(motif_enrich, qval_thresh=0.05, fig_ofn=NULL, ncol=4, nrow=4, show_legend=FALSE){
    motifs <- motif_enrich %>% filter(qval <= qval_thresh) %>% .$track %>% unique
    motif_enrich <- motif_enrich %>% filter(track %in% motifs)

    .plot_motif <- function(motif){
        motif_enrich %>% filter(track == motif) %>% select(clust, rel_enrich, abs_enrich, qval) %>% gather('type', 'val', -clust, -qval) %>% mutate(lab = if_else(qval <= qval_thresh & type == 'rel_enrich', '*', '')) %>% ggplot(aes(x=factor(clust), y=val, fill=type, label=lab)) + geom_col(position='dodge') + xlab('Cluster') + ylab('Enrichment') + scale_fill_manual(name = '', values = c("darkred", "darkblue"),  guide=show_legend) + geom_text(size=10) + ggtitle(motif)
    }
    plots <- map(motifs, .plot_motif)
    p <- cowplot::plot_grid(plotlist=plots)

    if (!is.null(fig_ofn)){
        cowplot::save_plot(fig_ofn, p, nrow=nrow, ncol=ncol)
    }

    return(p)
}

########################################################################
#' @export
gpwm.create_lowres_motif_tracks <- function(pattern, new_pattern, resolution, ...){
    src_motifs <- gtrack.ls(pattern)
    motifs <- gsub(glue('{pattern}\\.'), '', src_motifs)
    out_motifs <- glue('{new_pattern}.{motifs}')
    vtracks <- glue('v_max_{motifs}')
    walk2(src_motifs, vtracks, ~ gvtrack.create(.y, .x, func='max'))
    on.exit(walk(vtracks, gvtrack.rm))
    commands <- glue('gtrack.create("{out_motifs}", "xxx", "{vtracks}", iterator={resolution})')
    
    gdir.create(new_pattern)
    res <- gcluster.run2(command_list=commands, ...)

    return(res)
}

########################################################################
.gpwm.base_dir <- function(track)
{
    map(track, ~  c(gdir.cwd(), strsplit(.x, '.', fixed=TRUE)[[1]])) %>%
        map(~ do.call(file.path, as.list(.x))) %>%
        as_vector %>%
        return()
}
    
########################################################################
if (!exists('.gpwm_cache__', envir=globalenv())) {
    .gpwm_cache__ <<- new.env()
}
