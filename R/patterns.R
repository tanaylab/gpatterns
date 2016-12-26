# Tidy CpGs Functions ------------------------------------------------

########################################################################
#' Get tidy cpgs from a track
#'
#' @param track name of track
#' @param intervals intervals
#'
#' @return
#' @export
#'
#' @examples
gpatterns.get_tidy_cpgs <- function(track,
                                    intervals = NULL){
    if (!gpatterns.track_exists(track)){
        stop('track does not exists')
    }
    files <- .gpatterns.tidy_cpgs_files(track)

    .get_tcpgs <- function(files){
        res <- files %>%
            map_df(function(f) fread(
                qq('gzip -d -c @{f}'),
                colClasses = c(
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
                    qual = 'numeric')) %>%
                    tbl_df %>%
                    mutate(cg_pos = cg_pos, meth = ifelse(meth == 'Z', 1, 0)))
        if (0 == nrow(res)){
            return(NULL)
        }
        return(res)
    }

    .intervals2files <- function(intervals, files) .gpatterns.intervals2files(intervals=intervals, files=files, dir = paste0(.gpatterns.base_dir(track), '/tidy_cpgs/'))


    if (!is.null(intervals)){
        tidy_intervals <- .gpatterns.get_tidy_cpgs_intervals(track)
        if (!is.character(intervals)){
            if (all(
                unite(intervals, 'coord', (chrom:end))$coord %in%
                unite(tidy_intervals, 'coord', (chrom:end))$coord)
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


########################################################################
#' Apply a function on tidy cpgs (internally separated by coordinates)
#'
#' @param track name of track
#' @param f function to apply (needs to return a data frame)
#' @param split_by_bin apply on each genomic bin separately
#' @param parallel execute parallely
#' @param use_sge use sge cluster for parallelization
#' @param max_jobs maximal number of jobs for sge
#' @param verbose print jobs status
#' @param ... additonal parameters for gcluster.run2
#'
#' @return a data frame with \code{f} applied to every genomic bin of tidy_cpgs
#' of \code{track}
#'
#' @export
#'
#' @examples
gpatterns.apply_tidy_cpgs <- function(track,
                                      f,
                                      split_by_bin = TRUE,
                                      parallel = getOption("gpatterns.parallel"),
                                      use_sge = FALSE,
                                      max_jobs = 300,
                                      verbose = FALSE,
                                      ...){
    # make sure that empty intervals would return NULL
    f_null <- function(track, intervals=NULL){
        tcpgs <- gpatterns.get_tidy_cpgs(track, intervals=intervals)
        if (is.null(tcpgs)){
            return(NULL)
        }
        return(f(tcpgs))
    }

    if (!split_by_bin){
        return(f_null(track))
    }

    coords <- .gpatterns.get_tidy_cpgs_intervals(track)

    if (use_sge){
        commands <- plyr::alply(coords, 1, function(cr)
            qq('f_null(track,
               intervals = gintervals("@{cr$chrom}", @{cr$start}, @{cr$end}))') %>%
                gsub('\n', ' ', .))

        res <- gcluster.run2(command_list = commands,
                             max.jobs = max_jobs,
                             debug = verbose,
                             packages = c('stringr', 'gpatterns'),
                             jobs_title = 'gpatterns.apply_tidy_cpgs',
                             collapse_results = FALSE,
                             ...)

        res <- map(res, function(x) x$retv)
    } else {
        res <- plyr::alply(coords, 1, function(cr) {
            f_null(track, intervals=cr)
        }, .parallel=parallel,
        .progress='text')
    }

    tryCatch(res <- res %>% compact %>% map_df(~ .),
             error=function(e) {
                 message('wasn\'t able to collate_results, returning raw output')
             })


    return(res)
}


########################################################################
#' @export
.gpatterns.get_tidy_cpgs_intervals <- function(track){
    .gpatterns.tidy_cpgs_files(track) %>%
        basename %>%
        gsub('\\.tcpgs\\.gz$', '', .) %>%
        strsplit('_') %>%
        map_df(function(x) tibble(chrom=x[1], start=as.numeric(x[2]), end=as.numeric(x[3])))
}

# Pattern generation Functions ------------------------------------------------

########################################################################
#' Calculate pattern coverage for each CpG
#'
#' @param track name of track
#' @param pat_len pattern length
#' @param max_span maximum span to look for CpGs
#' @param by_coord_bin execute separatly for each coordinates bin
#'        (for tracks that do not fit into memory)
#' @param parallel execute parallely
#'
#' @return
#' @export
#'
#' @examples
gpatterns.get_pat_cov <- function(track,
                                  pat_len,
                                  max_span=500,
                                  by_coord_bin=TRUE,
                                  parallel = getOption("gpatterns.parallel")){
    if (by_coord_bin){
        return(track %>%
                   gpatterns.apply_tidy_cpgs(
                       function(x) .gpatterns.tidy_cpgs_2_pat_cov(x,
                                                                  pat_len=pat_len,
                                                                  max_span=max_span),
                       parallel=FALSE))
    } else {
        return(gpatterns.get_tidy_cpgs(track) %>%
                   .gpatterns.tidy_cpgs_2_pat_cov(pat_len=pat_len,
                                                  max_span=max_span))
    }

}

########################################################################
#' Extracts patterns given pattern space
#'
#' @param track name of track
#' @param pattern_space pattern space data frame
#' (output of gpatterns.tracks_to_pat_space or gpatterns.intervs_to_pat_space)
#' needs to have the following fields: chrom,start,end,fid
#' @param pattern_space_file file that holds the pattern space
#' @param max_missing maximal number of missing data ('*') in pattern
#' @param min_cov minimal pattern coverage
#' @param by_coord_bin execute separatly for each coordinates bin
#'        (for tracks that do not fit into memory)
#' @param ... additional parameters to gpatterns.apply_tidy_cpgs
#' @return
#' @export
#'
#' @examples
gpatterns.tidy_cpgs_to_pats <- function(track,
                                        pattern_space = NULL,
                                        pattern_space_file = NULL,
                                        max_missing = Inf,
                                        min_cov = 1,
                                        by_coord_bin=TRUE,
                                        parallel = getOption("gpatterns.parallel"),
                                        ...){
    if (!is.null(pattern_space_file)){
        pattern_space <- fread(pattern_space_file)
    } else if (is.null(pattern_space)){
        stop('Please to provide either pattern_space or pattern_space_file')
    }

    pat_space <- pattern_space %>%
        group_by(fid) %>%
        mutate(pat_pos=1:n()) %>%
        select(fid, chrom, cg_pos=start, pat_pos) %>%
        ungroup


    max_pat_len <- max(pat_space$pat_pos)

    .cpgs_to_pats <- function(cpgs){
        cpgs <- cpgs %>% select(read_id, chrom, meth, cg_pos)
        pats <- cpgs %>%
            inner_join(pat_space, by=c('chrom', 'cg_pos')) %>%
            reshape2::dcast(read_id + fid ~ pat_pos, value.var='meth', fill='*') %>%
            unite_('pat', paste(1:max_pat_len), sep='') %>%
            mutate(n_na = stringr::str_count(pat, '\\*')) %>%
            filter(n_na <= max_missing) %>%
            select(-n_na)

        return(pats)
    }

    if (by_coord_bin){
        pats <- track %>% gpatterns.apply_tidy_cpgs(function(x) .cpgs_to_pats(x),
                                                    parallel=parallel, ...)
    } else {
        pats <- track %>% gpatterns.get_tidy_cpgs() %>% .cpgs_to_pats()
    }

    pats <- pats %>%
        group_by(fid) %>%
        filter(n() >= min_cov) %>%
        select(fid, pattern=pat, read_id) %>%
        ungroup %>%
        arrange(fid, pattern, read_id)

    return(pats)
}

########################################################################
#' Generate pattern space
#'
#' @param tracks tracks
#' @param pat_len pattern length
#' @param max_span  maximum span to look for CpGs
#' @param exclude_intervals intervals to exclude
#' @param exclude_cpgs CpGs to exclude
#' @param min_cov minimal coverage per CpG
#'
#' @return
#' @export
#'
#' @examples
gpatterns.tracks_to_pat_space <- function(tracks,
                                          pat_len, max_span=NULL,
                                          exclude_intervals=NULL,
                                          exclude_cpgs=NULL,
                                          min_cov=1){
    if (all(gtrack.exists(.gpatterns.pat_cov_track_name(tracks, pat_len)))){
        cov_tracks <- .gpatterns.pat_cov_track_name(tracks, pat_len)
    } else {
        cov_tracks <- .gpatterns.cov_track_name(tracks)
    }

    covered_cpgs <- gscreen(
        paste(sprintf("!is.na(%s)", cov_tracks), collapse=' | '),
        intervals = .gpatterns.genome_cpgs_intervals,
        iterator = .gpatterns.genome_cpgs_intervals)
    expr <- sprintf("sum_ignore_na_all(%s)", paste(cov_tracks, collapse=','))

    covs <- gextract(expr,
                     intervals = covered_cpgs,
                     iterator = covered_cpgs,
                     colnames = 'cov') %>%
        select(-intervalID) %>%
        filter(cov >= min_cov)

    if (!is.null(exclude_cpgs)){
        covs <- covs %>% anti_join(exclude_cpgs, by=c('chrom', 'start', 'end'))
    }

    if (!is.null(exclude_intervals)){
        covs <- covs %>%
            ungroup %>%
            gintervals.neighbors1(exclude_intervals) %>%
            filter(dist != 0) %>%
            select(-(chrom1:dist))
    }

    if (!is.null(max_span)){
        covs <- covs %>% filter(abs(end - lead(start, 1)) <= max_span)
    }

    space <- covs %>%
        arrange(chrom, start, end) %>%
        group_by(chrom) %>%
        filter(n() >= pat_len) %>%
        slice(1:(n() - n() %% pat_len))

    space <- space %>% ungroup %>% arrange(chrom, start, end) %>%
        mutate(fid = map(1:(nrow(space)/pat_len), ~ rep(., pat_len)) %>% combine)
    return(space %>% ungroup %>% select(chrom, start, end, fid))
}

########################################################################
#' Generate pattern space for specific genomic loci
#'
#' @param tracks tracks
#' @param intervals genomic loci to cover
#' @param pat_len pattern length
#' @param max_span maximum span to look for CpGs
#' @param contiguous generate patterns that do not skip a CpG
#' @param add_rest add rest of the genome (using gpatterns.tracks_to_pat_space)
#' @param min_cov minimal CpG coverage
#'
#' @return
#' @export
#'
#' @examples
gpatterns.intervs_to_pat_space <- function(tracks,
                                           intervals,
                                           pat_len,
                                           max_span=500,
                                           contiguous=TRUE,
                                           add_rest=FALSE,
                                           min_cov=1){
    if (all(gtrack.exists(.gpatterns.pat_cov_track_name(tracks, pat_len)))){
        cov_tracks <- .gpatterns.pat_cov_track_name(tracks, pat_len)
    } else {
        cov_tracks <- .gpatterns.cov_track_name(tracks)
    }
    expr <- sprintf('sum_ignore_na_all(%s)', paste(cov_tracks, collapse=','))
    covs <- gextract.left_join(expr,
                               intervals=intervals %>%
                                   mutate(start.orig=start, end.orig=end) %>%
                                   gintervals.expand(max_span),
                               iterator=.gpatterns.genome_cpgs_intervals,
                               colnames='cov')

    # each cpg belongs to only one interval (if more than 1 choose the first)
    covs <- covs %>% distinct(chrom, start, end, .keep_all=TRUE)
    if (contiguous){
        covs <- covs %>%
            group_by(chrom1, start1, end1) %>%
            filter(n() >= pat_len) %>%
            arrange(chrom1, start1, end1, chrom, start, end) %>%
            mutate(cg_num=1:n()) %>%
            filter(sum(chrom == chrom1 & start == start.orig & end == end.orig) > 0) %>%
            mutate(pos_cg_num=which(chrom == chrom1 & start == start.orig & end == end.orig)) %>%
            filter(abs(cg_num - pos_cg_num) < pat_len) %>%
            select(-cg_num, -pos_cg_num)
        space <- covs %>%
            group_by(chrom1, start1, end1) %>%
            arrange(chrom1, start1, end1, chrom, start, end) %>%
            mutate(m = rollmean(cov, pat_len, align='left', fill=NA)) %>%
            mutate(cg_num=1:n(), chosen_start=which.max(m)) %>%
            filter(cg_num %in% seq(chosen_start[1], chosen_start[1] + pat_len - 1)) %>%
            select(chrom, start, end, chrom1, start1, end1, start.orig, end.orig)
    } else {
        space <- covs %>%
            group_by(chrom1, start1, end1) %>%
            top_n(pat_len, cov) %>%
            ungroup
    }
    space <- space %>%
        ungroup %>%
        arrange(chrom, start, end) %>%
        mutate(fid = map(1:(nrow(space)/pat_len), ~ rep(., pat_len)) %>% combine) %>%
        ungroup %>%
        select(chrom, start, end, fid)


    if (add_rest){
        fid_intervs <- space %>%
            group_by(fid) %>%
            summarise(chrom=first(chrom), start = min(start), end=max(end)) %>%
            select(-fid) %>%
            ungroup
        space_rest <- gpatterns.tracks_to_pat_space(tracks,
                                                    pat_len,
                                                    max_span,
                                                    exclude_intervals=fid_intervs,
                                                    exclude_cpgs=space %>%
                                                        ungroup %>%
                                                        select(chrom, start, end)
                                                    )
        space_rest <- space_rest %>% mutate(fid = paste0(fid, 'R'))
        space <- space %>%
            mutate(fid = paste(fid)) %>%
            bind_rows(space_rest)
        space <- space %>%
            arrange(chrom, start, end) %>%
            mutate(fid = map(1:(nrow(space)/pat_len), ~ rep(., pat_len)) %>% combine)
    }


    return(space %>% ungroup)

}


########################################################################
.gpatterns.tidy_cpgs_2_pat_cov  <- function(calls, pat_len, max_span){
    calls <- calls %>%
        mutate(start = cg_pos, end = start+1) %>%
        select(read_id, chrom, start, end)

    cgs <- gextract.left_join(.gpatterns.genome_cpgs_track,
                              intervals= calls %>%
                                  distinct(chrom, start, end) %>%
                                  mutate(start.orig=start, end.orig=end, end=end+max_span) %>%
                                  as.data.frame %>%
                                  gintervals.force_range(),
                              iterator=.gpatterns.genome_cpgs_intervals,
                              colnames='CG') %>% select(-CG)

    cgs <- cgs %>%
        group_by(chrom1, start1, end1) %>%
        slice(1:pat_len) %>%
        mutate(pid=paste0(min(start), '_', max(end)))

    cgs <- cgs %>% left_join(calls, by=c('chrom', 'start', 'end'))
    valid_cg_ids <- cgs %>%  filter(start == start.orig, end == end.orig) %>% ungroup %>% select(pid, read_id)

    res <- cgs %>%
        inner_join(valid_cg_ids, by=c('pid', 'read_id')) %>%
        filter(!(start1 == start.orig & end == end.orig)) %>%
        summarise(pat_cov = n()) %>%
        rename(chrom=chrom1, start=start1, end=end1) %>% ungroup %>%
        mutate(end = end - max_span)

    res <- calls %>% distinct(chrom, start, end) %>% left_join(res, by=c('chrom', 'start', 'end')) %>% replace_na(list(pat_cov=0))
    return(res)
}

########################################################################
#' Generate CpG pattern frequency from tidy_cpgs
#'
#' @param calls tidy_cpgs
#' @param pat_length length of the pattern, e.g. for pat_lenth == 2 patterns would
#' be 00,01,10,11
#' @param min_cov minimal coverage
#' @param tidy tidy
#'
#' @return
#' @export
#'
#' @examples
gpatterns.tidy_cpgs_2_pat_freq <- function(calls, pat_length = 2, min_cov = 1, tidy=TRUE){

        #########################################################
        gen_pats <- function(calls, min_cov, pat_length=2){
            message(qq('cpg num: 1'))
            for (i in 1:(pat_length - 1)){
                message(qq('cpg num: @{i+1}'))
                calls <- calls %>%
                    mutate(next_meth = lead(meth, i)) %>%
                    group_by(read_id) %>%
                    mutate(next_pos = lead(start, i)) %>%
                    ungroup %>%
                    filter(!is.na(next_pos)) %>%
                    gintervals.neighbors1(.gpatterns.genome_next_cpg_intervals) %>%
                    filter(dist == 0) %>%
                    filter(next_pos == nextcg) %>%
                    select(-(next_pos:dist)) %>%
                    rename_(.dots = setNames('next_meth', paste0('next_', i, '_meth')))
            }

            pats_dist <- calls %>%
                unite('pat', ends_with('meth'), sep='') %>%
                group_by(chrom, start, end, pat) %>%
                summarise(n_pat = n()) %>%
                group_by(chrom, start, end) %>%
                filter(sum(n_pat) >= min_cov) %>%
                mutate(p_pat = n_pat / sum(n_pat)) %>%
                ungroup %>%
                mutate(pat = factor(pat))

            return(pats_dist)
        }

        #########################################################
        fill_columns <- function(pats_dist, pat_length){
            possible_pats <- apply(expand.grid(rep(list(0:1), pat_length)), 1,
                                   function(x) paste0(x, collapse='')) %>%
                paste0('pat_', .)

            for (column in possible_pats){
                if (!(column %in% colnames(pats_dist))){
                    if (nrow(pats_dist) == 0){
                        pats_dist[, column] <- numeric()
                    } else {
                        pats_dist[, column] <- 0
                    }
                }
            }
            return(pats_dist[, c('chrom', 'start', 'end', possible_pats)])
        }

        #########################################################
        calls <-  calls %>%
            mutate(start = cg_pos, end = start + 1) %>%
            select(chrom, start, end, read_id, meth)

        message('generating pats...')
        pats_dist <- gen_pats(calls, min_cov=min_cov, pat_length=pat_length)
        if (nrow(pats_dist) == 0){
            return(fill_columns(pats_dist[, c('chrom', 'start', 'end')], pat_length))
        }

        if (tidy){
            pats_dist <- k %>%
                tidyr::complete(nesting(chrom, start, end), pat, fill=list(p_pat=0, n_pat=0))
        } else {
            pats_dist <- pats_dist %>%
                select(-p_pat) %>%
                mutate(pat = paste0('pat_', pat)) %>% spread('pat', 'n_pat')
            pats_dist[is.na(pats_dist)] <- 0
            pats_dist <- fill_columns(pats_dist, pat_length)
        }

        return(pats_dist)
}




########################################################################
#' Generate pileup form tidy_cpgs
#'
#' @param calls tidy_cpgs
#' @param dsn downsampling n
#'
#' @return data frame with chrom,start,end,meth,unmeth,cov,avg
#' @export
#'
#' @examples
gpatterns.tidy_cpgs_2_pileup <- function(calls, dsn = NULL){
    pileup <- calls %>%
        select(chrom, start=cg_pos, call=meth) %>%
        mutate(call = ifelse(call == 1, 'meth', 'unmeth'))

    if (!is.null(dsn)) {
        pileup <- pileup %>%
            group_by(chrom, start) %>%
            filter(n() >= dsn) %>%
            sample_n(dsn)
    }

    pileup <- pileup %>%
        group_by(chrom, start, call) %>%
        summarise(n = n()) %>%
        mutate(end = start + 1) %>%
        select(chrom, start, end, call, n) %>%
        spread('call', 'n', fill=0)

    if (!('meth' %in% colnames(pileup))){
        pileup$meth <- 0
    }
    if (!('unmeth' %in% colnames(pileup))){
        pileup$unmeth <- 0
    }

    pileup <- pileup %>%
        mutate(cov = meth + unmeth, avg = meth / cov) %>%
        select(chrom, start, end, meth, unmeth, cov, avg)

    return(pileup)
}

# Pattern utility Functions ------------------------------------------------

########################################################################
#' Downsample patterns table
#'
#' @param patterns patterns table
#' @param dsn downsampling n
#'
#' @return
#' @export
#'
#' @examples
gpatterns.downsample_patterns <- function(patterns, dsn){
    patterns %>%
        group_by(fid) %>%
        filter(n() >= dsn) %>%
        sample_n(dsn) %>%
        ungroup
}



########################################################################
#' transforms patterns to summary statistics: n,n0,n1,nx,nc,meth and epipoly
#'
#' @param patterns_tab patterns table (needs to have 'fid' and 'pattern' fields)
#' @param noise_threshold thershold of a pattern to be considered 'noise'
#'
#' @return
#' @export
#'
#' @examples
gpatterns.frag_stats <- function(patterns_tab, noise_threshold=0.2){
    patterns_tab %>%
        mutate(cpgs = nchar(pattern), ones = cpgs-nchar(gsub('1', '', pattern))) %>%
        group_by(fid) %>%
        summarize(ncpg = min(cpgs),
                  n = n(),
                  n0 = sum(ones == 0),
                  n1 = sum(ones == cpgs),
                  nx = .gpatterns.count_noise(pattern, cpgs, ones, noise_threshold),
                  nc = n - n0 - n1 - nx,
                  pat_meth = sum(ones) / sum(cpgs),
                  epipoly = gpatterns.calc_epipoly(pattern)) %>%
        ungroup()
}

########################################################################
#' Extract patterns of a track
#'
#' @param track name of track
#' @param fids extract for specific fragments
#' @param tidy if TRUE: returns tidy format, else - a list
#' @param dsn extract downsampled patterns
#'
#' @return
#' @export
#'
#' @examples
gpatterns.extract_patterns <- function(track, fids = NULL, tidy = TRUE, dsn = NULL)
{
    if (!is.null(dsn)){
        track <- .gpatterns.downsampled_track_name(track, dsn)
    }

    if (!.gpatterns.patterns_exist(track)){
        stop(sprintf('no patterns available for %s', track))
    }

    patterns_tab <-  .gpatterns.load_patterns_tab(track)

    if(!is.null(fids)){
        patterns_tab <- patterns_tab %>% filter(fid %in% fids)
    }

    if (nrow(patterns_tab) == 0){
        warning('no patterns in fid')
    }
    if (tidy){
        return(patterns_tab %>% tbl_df)
    }

    patterns_tab <- patterns_tab %>%
        group_by(fid) %>%
        do(patterns=.$pattern) %>%
        ungroup()

    patterns_list <- patterns_tab$patterns
    names(patterns_list) <- paste0('fid', patterns_tab$fid)
    return(patterns_list)
}

########################################################################
#' Extract pattern data from tracks without specifying coordinates
#'
#' @param ... tracks (comma separated)
#' @param samples name for each track
#' @param elements elements to extract
#' @param add_bipolar_stats add bipolar stats
#' @param extract downsampled data
#' @param na.rm remove rows with missing data
#'
#' @return
#' @export
#'
#' @examples
gpatterns.extract_all <- function(...,
                                  samples=NULL,
                                  elements = c('fid', 'ncpg', 'n', 'n0', 'n1', 'nx', 'pat_meth', 'epipoly'),
                                  add_bipolar_stats = FALSE,
                                  dsn = NULL,
                                  na.rm = FALSE){
    if (!('fid' %in% elements)){
        elements <- c('fid', elements)
    }

    if (add_bipolar_stats){
        elements = c(elements, .gpatterns.bipolar_model_stats)
    }

    tracks <- sapply(list(...), as.character)

    samples <- samples %||% tracks

    if (!is.null(dsn)){
        tracks <- .gpatterns.downsampled_track_name(tracks, dsn)
    }

    tab <- map2_df(
        tracks,
        samples,
        function(track, name)
            .gpatterns.load_table(
                saved_name = .gpatterns.fids_tab_name(track),
                file = .gpatterns.fids_file_name(track))
        %>% mutate(track = name))

    if (any(elements %in% .gpatterns.bipolar_model_stats)){
        mix_tab <- map2_df(
            tracks,
            samples,
            function(track, name)
                .gpatterns.load_table(
                    saved_name = .gpatterns.bipolar_model_tab_name(track),
                    file = .gpatterns.bipolar_model_file_name(track))
            %>% mutate(track=name))
        tab <- tab %>%
            left_join(mix_tab, by='fid')
    }

    tab <- tab %>% select_(.dots=c('chrom', 'start', 'end', 'track', elements))
    if (na.rm){
        tab <- tab %>% na.omit()
    }
    return(tab)
}

########################################################################
#' Extract pattern data
#'
#' @param ... tracks to extract
#' @param intervals genomic scope to extract from intervals
#' @param colnames names for the tracks
#' @param elements elements to extract. default: fid,ncpg,n,n0,n1,nx,nc,pat_meth
#' @param add_bipolar_stats add stats from bipolar model (if computed)
#' @param add_patterns add patterns
#' @param add_cpg_content add centered cpg content
#' @param cpg_nhood scope for centered cpg content (in bp)
#' @param add_coordinates add chrom,start,end for fragments
#' @param tidy return tidy output
#' @param dsn extract downsampled data
#' @param na.rm remove rows with missing data
#'
#' @return
#' @export
#'
#' @examples
gpatterns.extract <- function(...,
                              intervals,
                              colnames=NULL,
                              elements = c('fid', 'ncpg', 'n', 'n0', 'n1', 'nx', 'nc', 'pat_meth'),
                              add_bipolar_stats = FALSE,
                              add_patterns = FALSE,
                              add_cpg_content = FALSE,
                              cpg_nhood=500,
                              add_coordinates = TRUE,
                              tidy = TRUE,
                              dsn = NULL,
                              na.rm = TRUE)
{

    if (add_coordinates && any(!(c('chrom', 'start', 'end') %in% elements))){
        elements <- c('chrom', 'start', 'end', elements)
    }
    if (!('fid' %in% elements)){
        elements <- c('fid', elements)
    }

    if (add_bipolar_stats){
        elements = c(elements, .gpatterns.bipolar_model_stats)
    }

    if (add_patterns){
        elements = c(elements, 'pattern')
    }

    if (add_cpg_content){
        elements = c(elements, 'cpg_content')
    }

    tracks <- sapply(list(...), as.character)
    colnames <- colnames %||% tracks

    if (!is.null(dsn)){
        tracks <- .gpatterns.downsampled_track_name(tracks, dsn)
    }

    add_cpg_content <- function(fids_tab, cpg_nhood){
        cg_cont <- fids_tab %>%
            distinct(chrom, start, end) %>%
            mutate(cpg_content = .gpatterns.centered_cpg_content(., cpg_nhood = cpg_nhood))
        return(fids_tab %>% left_join(cg_cont, by=c('chrom', 'start', 'end')))
    }

    tab <- map2_df(tracks, colnames, function(track, track_name){
        fids_scope <- gextract(.gpatterns.fid_track_name(track),
                               intervals = intervals, colnames = 'fid') %>%
            select(fid, intervalID)

        fids_tab <- .gpatterns.load_fids_tab(track) %>%
            inner_join(fids_scope, by = 'fid')

        if ('cpg_content' %in% elements){
            fids_tab <- fids_tab %>% add_cpg_content(cpg_nhood)
        }

        if ('pattern' %in% elements) {
            patterns_tab <- gpatterns.extract_patterns(track, fids=fids_tab$fid, tidy = TRUE)
            fids_tab <- fids_tab %>% left_join(patterns_tab, by = 'fid')
        }

        if (any(elements %in% .gpatterns.bipolar_model_stats)){
            bipolar_tab <- .gpatterns.load_bipolar_tab(track)
            fids_tab <- fids_tab %>% left_join(bipolar_tab, by = 'fid')
        }

        fids_tab <- fids_tab %>%
            mutate(samp = track_name) %>%
            select_(.dots=c('samp', elements, 'intervalID'))
    })

    if (na.rm){
        tab <- tab %>% na.omit()
    }

    if (tidy){
        return(tab)
    }

    spread_elements <- elements[which(!(elements %in% c('chrom', 'start', 'end', 'fid')))]
    untidy_tab <- tab %>%
        select(one_of(c('chrom', 'start', 'end', 'fid', 'samp', spread_elements[1]))) %>%
        spread_('samp', spread_elements[1])

    add_element_suffix <- function(tb, element){
        column_nums <- which(colnames(tb) %in% colnames)
        colnames(tb)[column_nums] <- paste0(colnames(tb)[column_nums], '.', element)
        return(tb)
    }

    untidy_tab <- untidy_tab %>% add_element_suffix(untidy_tab, spread_elements[1])

    for (element in spread_elements[-1]){
        untidy_tab1 <- tab %>%
            select(one_of(c('chrom', 'start', 'end', 'fid', 'samp', element))) %>%
            spread_('samp', element)
        untidy_tab1 <- untidy_tab %>% add_element_suffix(untidy_tab1, element)
        untidy_tab <- suppressMessages(untidy_tab %>% full_join(untidy_tab1))
    }

    return(untidy_tab)
}

########################################################################
.gpatterns.save_table <- function(x, saved_name, file)
{
    assign(saved_name, data.table::as.data.table(x))
    save(list=saved_name, file=file)
}


########################################################################
.gpatterns.load_table <- function(saved_name, file)
{
    load(file)
    return(get(saved_name) %>% tbl_df())
}

########################################################################
.gpatterns.load_fids_tab <- function(track){
    .gpatterns.load_table(saved_name=.gpatterns.fids_tab_name(track),
                          file=.gpatterns.fids_file_name(track))
}

########################################################################
.gpatterns.load_patterns_tab <- function(track){
    .gpatterns.load_table(saved_name=.gpatterns.patterns_tab_name(track),
                          file=.gpatterns.patterns_file_name(track))
}

########################################################################
.gpatterns.load_bipolar_tab <- function(track){
    .gpatterns.load_table(saved_name=.gpatterns.bipolar_model_tab_name(track),
                          file=.gpatterns.bipolar_model_file_name(track))
}

########################################################################
.gpatterns.patterns_exist <- function(tracks){
    map_lgl(tracks, function(track) file.exists(.gpatterns.patterns_file_name(track)))
}

########################################################################
.gpatterns.count_noise <- function(patterns, cpgs, ones, noise_threshold)
{
    noise_threshold <- noise_threshold * length(patterns)
    patterns <- patterns[(ones != 0) & (ones != cpgs)]
    tab      <- tabulate(as.factor(patterns))
    return(sum(tab[tab < noise_threshold]))
}

########################################################################
.gpatterns.remove.tracks <- function(tracks){
    for (tr in tracks){
        if(gtrack.exists(tr)){
            message(sprintf("removing %s", tr))
            gtrack.rm(tr, force=TRUE)
        }
    }
}

########################################################################
.gpatterns.centered_cpg_content <- function(intervals, cpg_nhood)
{
    intervals <- intervals %>%
        select(chrom, start, end) %>%
        mutate(start = floor((start+end-1)/2), end = start+1) %>%
        mutate(start= start - cpg_nhood, end = end + cpg_nhood) %>%
        as.data.frame %>%
        gintervals.force_range()

    intervals <- gvextract(.gpatterns.genome_cpgs_track,
                           intervals = intervals,
                           iterator = intervals,
                           colnames='cpg_content',
                           func = 'sum')

    return(intervals$cpg_content / (2*cpg_nhood + 1))
}


