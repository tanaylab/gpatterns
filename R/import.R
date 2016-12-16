# import Functions ------------------------------------------------

########################################################################
.step_invoke <- function(step, steps, f, ...){
    if (step %in% steps){
        message(qq('doing @{step}'))
        return(f(...))
    } else {
        message(qq('skipping @{step}'))
    }
}

########################################################################
#' Creates a track from bam files.
#'
#' @description Creates a track from bam files for a specific sample.
#' Use this methods only for small datasets. For large datasets please
#' the Snakemake pipeline.
#'
#' @param bams character vector with path of bam files
#' @param workdir directory in which the files would be saved
#' @param track name of the track to generate
#' @param description description of the track to generate
#' @param steps steps of the pipeline to do. Possible options are:
#' 'bam2tidy_cpgs', 'filter_dups', 'bind_tidy_cpgs', 'pileup', 'pat_freq', 'pat_cov'
#' @param conversion bisulfite conversion. could be 'ct' or 'ga'
#' @param paired_end bam files are paired end, with R1 and R2 interleaved
#' @param dsn downsampling n. Leave NULL for no downsampling
#' @param pat_cov_lens lengthes of patterns to calculate pattern coverage track for
#' @param max_span maximal span to look for patterns (usually the maximal insert length)
#' @param pat_freq_len lengthes of patterns to calculate pattern frequency track for
#' @param nbins number of genomic bins to separate the analysis.
#' @param groot root of misha genomic database to save the tracks
#' @param use_sge use sun grid engine for parallelization
#' @param max_jobs maximal number of jobs for sge parallelization
#' @param parallel parallelize using threads (number of threads is determined by gpatterns.set_parallel)
#'
#'
#' @return if 'stats' is one of the steps - data frame with statistics. Otherwise none.
#' @export
#'
#' @examples
gpatterns.import_from_bam <- function(bams,
                                      workdir,
                                      track,
                                      description,
                                      steps = c('bam2tidy_cpgs', 'filter_dups', 'bind_tidy_cpgs', 'pileup', 'pat_freq', 'pat_cov'),
                                      conversion = 'ct',
                                      paired_end = TRUE,
                                      dsn = NULL,
                                      pat_cov_lens = c(3,5,7),
                                      max_span = 500,
                                      pat_freq_len = 2,
                                      nbins = nrow(gintervals.all()),
                                      groot = GROOT,
                                      use_sge = FALSE,
                                      max_jobs = 400,
                                      parallel = getOption('gpatterns.parallel')){
    gsetroot(groot)
    all_steps <- c('bam2tidy_cpgs', 'filter_dups', 'bind_tidy_cpgs', 'pileup', 'pat_freq', 'pat_cov', 'stats')

    if (steps == 'all'){
        steps <- all_steps
    }

    stopifnot(all(steps %in% c(all_steps)))

    # get genomic bins
    genomic_bins <- gbin_intervals(intervals = gintervals.all(), nbins)

    # bam to tidy cpgs
    .step_invoke(
        'bam2tidy_cpgs',
        steps,
        .gpatterns.bam2tidy_cpgs,
        bams,
        tidy_cpgs_dir = qq('@{workdir}/tidy_cpgs'),
        stats_dir = qq('@{workdir}/tidy_cpgs/stats'),
        genomic_bins = genomic_bins,
        conversion = conversion,
        paired_end = paired_end,
        use_sge = use_sge,
        max_jobs = max_jobs,
        parallel = parallel)

    # filter dups
    .step_invoke(
        'filter_dups',
        steps,
        .gpatterns.filter_dups,
        tidy_cpgs_dir = qq('@{workdir}/tidy_cpgs'),
        stats_dir = qq('@{workdir}/tidy_cpgs_uniq/stats'),
        uniq_tidy_cpgs_dir = qq('@{workdir}/tidy_cpgs_uniq'),
        genomic_bins = genomic_bins,
        paired_end = paired_end,
        use_sge = use_sge,
        max_jobs = max_jobs,
        parallel = parallel)

    # bind tidy cpgs
    .step_invoke(
        'bind_tidy_cpgs',
        steps,
        .gpatterns.bind_tidy_cpgs,
        uniq_tidy_cpgs_dir = qq('@{workdir}/tidy_cpgs_uniq'),
        track = track)

    # pileup
    .step_invoke(
        'pileup',
        steps,
        .gpatterns.pileup,
        track = track,
        description = description,
        genomic_bins = genomic_bins,
        dsn = NULL,
        use_sge = use_sge,
        max_jobs = max_jobs,
        parallel = parallel)


    # pat_freq
    .step_invoke(
        'pat_freq',
        steps,
        .gpatterns.pat_freq,
        track = track,
        description = description,
        pat_freq_len = pat_freq_len,
        use_sge = use_sge,
        max_jobs = max_jobs,
        parallel = parallel)


    # pat_cov
    .step_invoke(
        'pat_cov',
        steps,
        .gpatterns.pat_cov,
        track = track,
        description = description,
        pat_cov_lens = pat_cov_lens,
        max_span = max_span,
        use_sge = use_sge,
        max_jobs = max_jobs,
        parallel = parallel)

    .step_invoke(
        'stats',
        steps,
        gpattens.get_pipeline_stats,
        track = track,
        tidy_cpgs_stats_dir = qq('@{workdir}/tidy_cpgs/stats'),
        uniq_tidy_cpgs_stats_dir = qq('@{workdir}/tidy_cpgs_uniq/stats'))
}

########################################################################
#' Get QC statistics for track
#'
#' @param track track name
#' @param tidy_cpgs_stats_dir directory with tidy_cpgs stats files
#' @param uniq_tidy_cpgs_stats_dir directory with filter_dups stats files
#' @param add_mapping_stats add full mapping statistics (singles, discordant etc.)
#'
#' @return
#' @export
#'
#' @examples
gpattens.get_pipeline_stats <- function(track,
                                        tidy_cpgs_stats_dir,
                                        uniq_tidy_cpgs_stats_dir,
                                        add_mapping_stats = FALSE){
    uniq_stats <- list.files(uniq_tidy_cpgs_stats_dir, full.names=T) %>%
        map_df(~ fread(.x)) %>%
        summarise(total_reads = sum(total_reads), uniq_reads = sum(uniq_reads)) %>%
        mutate(uniq_frac = uniq_reads / total_reads)
    tidy_cpgs_stats <- list.files(tidy_cpgs_stats_dir, full.names=T) %>%
        map_df(~ fread(.x)) %>%
        slice(1)
    stats <- tidy_cpgs_stats %>%
        summarise(total_reads = sum(.),
                  mapped_reads = good + single_R1 + single_R2,
                  mapped_frac = mapped_reads / total_reads) %>%
        bind_cols(uniq_stats %>% rename(good_reads = total_reads))
    sm <- gsummary(qq('@{track}.cov'))
    stats[['cg_num']] <- sm[1]
    stats[['meth_calls']] <- sm[5]
    stats[['global_avg_meth']] <- gsummary(qq('@{track}.avg'))[6]

    stats <- stats %>% bind_cols(gpatterns.apply_tidy_cpgs(track, function(x) x %>% summarise(insert_len = mean(abs(insert_len), na.rm=T), n = n())) %>% mutate(f = insert_len * n / sum(n)) %>% summarise(insert_len = sum(f)))

    if (add_mapping_stats){
        stats <- stats %>% bind_cols(tidy_cpgs_stats)
    }

    return(stats)
}

########################################################################
.gpatterns.bam2tidy_cpgs <- function(bams, tidy_cpgs_dir, stats_dir, genomic_bins, conversion, paired_end = TRUE, bin = .gpatterns.bam2tidy_cpgs_bin, ...){
    walk(c(tidy_cpgs_dir, stats_dir), ~ system(qq('mkdir -p @{.x}')))

    bam_prefix <- if (1 == length(bams)) 'cat' else 'samtools cat'
    single_end <- if (!paired_end) '--single_end' else ''
    commands <- genomic_bins %>% by_row( function(gbins){
        stats_fn <- qq('@{stats_dir}/@{gbins$chrom}_@{gbins$start}_@{gbins$end}.stats')
        output_fn <- qq('@{tidy_cpgs_dir}/@{gbins$chrom}_@{gbins$start}_@{gbins$end}.tcpgs.gz')
        qq('ulimit -u 16384;
           @{bam_prefix} @{paste(bams, collapse=\'\')} |
           @{bin} -i - -o - -s @{stats_fn} --conversion @{conversion}
           @{single_end} --chrom @{gbins$chrom}
           --genomic-range @{gbins$start} @{gbins$end} |
           awk \'NR==1; NR > 1 {print $0 | "sort --field-separator=, -k2,7 -k1 -k9"}\' |
           gzip -c > @{output_fn}') %>%
            gsub('\n', '', .)
        }, .collate =  'cols', .to = 'cmd')
    .gpatterns.run_command(commands, jobs_title = 'bam2tidy_cpgs', ...)
}

########################################################################
.gpatterns.filter_dups <- function(tidy_cpgs_dir, stats_dir, uniq_tidy_cpgs_dir, genomic_bins, paired_end = TRUE, bin = .gpatterns.filter_dups_bin, ...){
    walk(c(uniq_tidy_cpgs_dir, stats_dir), ~ system(qq('mkdir -p @{.x}')))
    all_files <- list.files(tidy_cpgs_dir, pattern='.tcpgs.gz', full.names = TRUE)
    files <- .gpatterns.intervals2files(genomic_bins, all_files, qq('@{tidy_cpgs_dir}/'))

    if (!all(all_files %in% files)){
        stop(sprintf('missing the following tidy cpgs files:\n\t%s',
                     paste(all_files[!(all_files %in% files)] ,collapse='\n\t')))
    }
    single_end <- if (!paired_end) '--only_R1' else ''

    commands <- genomic_bins %>% by_row(function(gbins){
        stats_fn <- qq('@{stats_dir}/@{gbins$chrom}_@{gbins$start}_@{gbins$end}.stats')
        input_fn <- qq('@{tidy_cpgs_dir}/@{gbins$chrom}_@{gbins$start}_@{gbins$end}.tcpgs.gz')
        output_fn <- qq('@{uniq_tidy_cpgs_dir}/@{gbins$chrom}_@{gbins$start}_@{gbins$end}.tcpgs.gz')
        qq('ulimit -u 16384; @{bin} -i @{input_fn} -o - -s @{stats_fn} @{single_end} | gzip -c > @{output_fn}')
        }, .collate =  'cols', .to = 'cmd')

    .gpatterns.run_command(commands, jobs_title = 'filter_dups', ...)

}

########################################################################
.gpatterns.bind_tidy_cpgs <- function(uniq_tidy_cpgs_dir, track = track, genomic_bins = genomic_bins){
    track_path <- .gpatterns.base_dir(track)
    system(qq('mkdir -p @{track_path}'))
    system(qq('ln -sf @{uniq_tidy_cpgs_dir} @{track_path}/tidy_cpgs'))
}

########################################################################
.gpatterns.pileup <- function(track, description, dsn = NULL, ...){
    message('calculating pileup...')
    pileup <- gpatterns.apply_tidy_cpgs(track, function(x) gpatterns.tidy_cpgs_2_pileup(x, dsn=dsn), ...)
    message('importing pileup to misha...')
    .gpattern.import_intervs_table(track, description, pileup, columns = c('meth', 'unmeth', 'cov', 'avg'))
}

########################################################################
.gpatterns.pat_freq <- function(track, description, pat_freq_len, ...){
    message(qq('calculating pattern frequency (@{pat_freq_len})...'))
    pat_freq <- gpatterns.apply_tidy_cpgs(track,
                                          function(x)
                                              gpatterns.tidy_cpgs_2_pat_freq(x,
                                                                             pat_length=pat_freq_len,
                                                                             tidy=FALSE),
                                          ...)
    if (nrow(pat_freq) == 0){
        warning('no patterns were found. Skipping creation of pattern frequency tracks')
    } else {
        message('importing pileup to misha...')
        .gpattern.import_intervs_table(track, description, pat_freq, columns = NULL)
    }
}

########################################################################
.gpatterns.pat_cov <- function(track, description, pat_cov_lens, max_span = 500, ...){
    message('calculating pattern coverage...')
    pat_cov <- pat_cov_lens %>%
        map_df(~ gpatterns.apply_tidy_cpgs(track,
                                           function(x) .gpatterns.tidy_cpgs_2_pat_cov(x,
                                                                                      pat_len=.x,
                                                                                      max_span=max_span),
                                           ...) %>%
                   mutate(pat_len = qq('pat_cov@{.x}')) )
    pat_cov <- pat_cov %>% spread(pat_len, pat_cov)

    message('importing pat_cov to misha...')
    .gpattern.import_intervs_table(track, description, pat_cov, columns = paste0('pat_cov', pat_cov_lens))
}

########################################################################
.gpattern.import_intervs_table <- function(track_pref, description, tab, columns = NULL){
    tab <- tbl_df(tab)
    columns <- columns %||% colnames(tab)[!(colnames(tab) %in% c('chrom', 'start', 'end'))]
    walk(columns, function(col){
        track_name <- qq('@{track_pref}.@{col}')
        if (gtrack.exists(track_name)){
            message(qq('removing @{track_name}'))
            gtrack.rm(track_name, force=TRUE)
        }
        message(qq('creating @{track_name}'))
        gtrack.create_sparse(track = track_name,
                             description = qq('@{description}: @{col} track'),
                             intervals = tab %>% select(chrom, start, end),
                             value = tab[[col]])
    })
}

########################################################################
.gpatterns.run_command <- function(commands, use_sge=FALSE, max_jobs=400, parallel = getOption('gpatterns.parallel'), jobs_title=''){
      if (use_sge){
        command_list <- 1:length(commands$cmd) %>% map(~ qq('system(commands$cmd[@{.x}])'))
        res <- gcluster.run2(command_list = command_list,
                             max.jobs = max_jobs,
                             packages = 'gpatterns',
                             jobs_title = jobs_title,
                             collapse_results = FALSE)
        codes <-  res %>% map_int(~ .x$retv)
    } else {
        res <- commands %>% plyr::alply(1, function(x) system(x$cmd), .parallel=parallel)
        codes <- res
    }
    walk2(codes, commands$cmd, function(code, cmd) if (code != 0) stop(qq('command "@{cmd}" failed')))
}



# Patterns import Functions ------------------------------------------------

########################################################################
#' Creates patterns attributes of track
#'
#' @param track track name
#' @param description a character string description
#' @param pat_space pattern space data frame
#' (output of gpatterns.tracks_to_pat_space or gpatterns.intervs_to_pat_space)
#' needs to have the following fields: chrom,start,end,fid
#' @param patterns_tab table with patterns (output of gpatterns.tidy_cpgs_to_pats)
#' needs to have the following fields: fid,pattern
#' @param add_read_id save the read_id together with the patterns
#' @param noise_threshold threshold to consider pattern as 'noise'
#' @param overwrite overwrite existing tracks
#' @param canonize convert to canonic form (see: \code{\link[gintervals.canonic]{misha}}).
#' warning: may lose data
#' @param add_biploar_stats run mixture model and add bipolaritly stats.
#' (warning: heavy computation)
#' @param ... additional parameters for gpatterns.calc_bipolarity
#'
#' @return
#' @export
#'
#' @examples
gpatterns.create_patterns_track <- function(track,
                                   description,
                                   pat_space,
                                   patterns_tab = NULL,
                                   add_read_id = TRUE,
                                   noise_threshold = 0.20,
                                   overwrite = TRUE,
                                   canonize = FALSE,
                                   add_biploar_stats = FALSE,
                                   ...){ #
    stopifnot( all(c('fid', 'chrom', 'start') %in% colnames(pat_space)) )

    if (is.null(patterns_tab)){
        message('extracting patterns from cpgs')
        patterns_tab <- gpatterns.tidy_cpgs_to_pats(track,
                                                    pat_space,
                                                    max_missing = 0,
                                                    min_cov = 1)
    }
    stopifnot( all(c('fid', 'pattern') %in% colnames(patterns_tab)) )


    # Create base dir for tracks
    dir.create(.gpatterns.base_dir(track), showWarnings=FALSE, recursive=TRUE)

    message("creating tables...")

    # Create a table mapping pattern positions to the fid
    loci_tab <- pat_space %>%
        group_by(fid, chrom) %>%
        summarize(start=min(start), end=max(end)+1) %>%
        ungroup() %>%
        select(chrom, start, end, fid)

    .gpatterns.create_tracks(
        track = track,
        description = description,
        patterns_tab = patterns_tab,
        loci_tab = loci_tab,
        overwrite = overwrite,
        add_read_id = add_read_id,
        noise_threshold = noise_threshold,
        add_biploar_stats = add_biploar_stats,
        canonize = canonize,
        ...)
}

########################################################################
gpatterns.create_downsampled_track <- function(track,
                                   dsn,
                                   description = NULL,
                                   patterns_tab = NULL,
                                   add_read_id = TRUE,
                                   noise_threshold = 0.20,
                                   overwrite = TRUE,
                                   canonize = FALSE,
                                   add_biploar_stats = FALSE,
                                   ...){

    stopifnot(gtrack.exists(.gpatterns.fid_track_name(track)))

    if (is.null(patterns_tab)){
        patterns_tab <- gpatterns.extract_patterns(track)
    }

    patterns_tab <- patterns_tab %>% gpatterns.downsample_patterns(dsn)


    track_ds <- gpatterns.downsampled_track_name(track, dsn)

    # Create base dir for tracks
    dir.create(.gpatterns.base_dir(track_ds), showWarnings=FALSE, recursive=TRUE)

    message("creating tables...")
    loci_tab <- .gpatterns.load_fids_tab(track) %>% select(chrom, start, end, fid)
    if (is.null(description)){
        description <- gtrack.attr.get(.gpatterns.fid_track_name(track), 'description')
    }

    .gpatterns.create_tracks(
        track = track_ds,
        description = description,
        patterns_tab = patterns_tab,
        loci_tab = loci_tab,
        overwrite = overwrite,
        add_read_id = add_read_id,
        noise_threshold = noise_threshold,
        add_biploar_stats = add_biploar_stats,
        canonize = canonize,
        ...)
}


########################################################################
.gpatterns.create_tracks <- function(track,
                                     description,
                                     patterns_tab,
                                     loci_tab,
                                     overwrite = FALSE,
                                     add_read_id = FALSE,
                                     add_biploar_stats = FALSE,
                                     noise_threshold = 0.2,
                                     canonize = FALSE,
                                     ...) {

     # Sort patterns table
    patterns_tab <- patterns_tab %>%
        mutate(cpgs=nchar(pattern), ones=cpgs-nchar(gsub('1', '', pattern))) %>%
        arrange(fid, cpgs, ones, pattern)

    if (add_read_id){
        .gpatterns.save_table(patterns_tab %>% select(fid, pattern, read_id),
                              saved_name=.gpatterns.patterns_tab_name(track),
                              file=.gpatterns.patterns_file_name(track))
    } else {
        .gpatterns.save_table(patterns_tab %>% select(fid, pattern),
                              saved_name=.gpatterns.patterns_tab_name(track),
                              file=.gpatterns.patterns_file_name(track))
    }


    # Create table tracks holding per fid statistics -
    # number of covered CpGs, read depth, pattern counts, avg methylation and noise
    fids_tab <- patterns_tab %>%
        gpatterns.frag_stats(noise_threshold = noise_threshold)

    fids_tab <- loci_tab %>%
        left_join(fids_tab, by='fid')

    .gpatterns.save_table(fids_tab,
                          saved_name=.gpatterns.fids_tab_name(track),
                          file=.gpatterns.fids_file_name(track))

    if (overwrite){
        tracks <- c(
            .gpatterns.fid_track_name(track),
            .gpatterns.ncpg_track_name(track),
            .gpatterns.n_track_name(track),
            .gpatterns.n0_track_name(track),
            .gpatterns.n1_track_name(track),
            .gpatterns.nx_track_name(track),
            .gpatterns.nc_track_name(track),
            .gpatterns.epipolymorphism_track_name(track),
            .gpatterns.pat_meth_track_name(track))
        .gpatterns.remove.tracks(tracks)
        if (file.exists(paste0(gdir.cwd(), '/', gsub('\\.', '/', track), '.track') )){
            gdir.rm(paste0(gdir.cwd(), '/', gsub('\\.', '/', track), '.track'), force=TRUE, recursive=TRUE)
        }
    }

    if (canonize){
        warning("Canonizing fragment intervals. May lose some regions")
        fids_tab <- fids_tab %>%
            gintervals.canonic %>%
            left_join(fids_tab %>%
                          distinct(chrom, start, end, .keep_all=TRUE),
                      by=c('chrom', 'start', 'end'))
    }

    message("creating tracks...")
    gtrack.create_sparse(.gpatterns.fid_track_name(track),  description, fids_tab, fids_tab$fid)
    gtrack.create_sparse(.gpatterns.ncpg_track_name(track), description, fids_tab, fids_tab$ncpg)
    gtrack.create_sparse(.gpatterns.n_track_name(track),    description, fids_tab, fids_tab$n)
    gtrack.create_sparse(.gpatterns.n0_track_name(track),   description, fids_tab, fids_tab$n0)
    gtrack.create_sparse(.gpatterns.n1_track_name(track),   description, fids_tab, fids_tab$n1)
    gtrack.create_sparse(.gpatterns.nx_track_name(track),   description, fids_tab, fids_tab$nx)
    gtrack.create_sparse(.gpatterns.nc_track_name(track),   description, fids_tab, fids_tab$nc)
    gtrack.create_sparse(.gpatterns.pat_meth_track_name(track), description, fids_tab, fids_tab$pat_meth)
    gtrack.create_sparse(.gpatterns.epipolymorphism_track_name(track), description, fids_tab, fids_tab$epipoly)

    if (add_biploar_stats){
        message("getting bipolarity stats")
        gpatterns.calc_bipolarity(track,
                                  overwrite=overwrite,
                                  ...)
    }
}

