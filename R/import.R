# standard tracks generation ------------------------------------------------
.gpatterns.gen_cg_track <- function(){

}

.get_nucs_df <- function(intervals, nucs=c('T', 'C', 'G', 'A')){   
    dinucs <- expand.grid(nucs, nucs)  %>% unite(Var1, Var2, col='dinucs', sep='') %>% .$dinucs 
    nuc_df <- tibble(start = seq(intervals$start, intervals$end - 1, 1), seq=(gseq.extract(intervals) %>% toupper %>% str_split('') %>% .[[1]])) %>% mutate(next_nuc = lead(seq)) %>% unite(seq, next_nuc, col='dinuc', sep='', remove=F) %>% mutate(chrom = intervals$chrom, end = start + 1) %>% select(chrom, start, end, seq, dinuc)
    for (nuc in nucs){
        nuc_df[nuc] <- as.numeric(nuc_df$seq == nuc)
    }
    for (nuc in dinucs){
        nuc_df[nuc] <- as.numeric(nuc_df$dinuc == nuc)
    }
    nuc_df <- nuc_df %>% select(-seq, -dinuc)
    return(nuc_df)        
}

.gpatterns.gen_next_CG_track <- function(){
    df <- map_df(gintervals.all()$chrom, ~ gintervals.load(.gpatterns.genome_cpgs_intervals, chrom=.x) %>%  mutate(nextcg = lead(start))) %>% arrange(chrom, start, end)
    gintervals.save(.gpatterns.genome_next_cpg_intervals, df)
}

.gpatterns.gen_re_fragments <- function(){

}

.gpatterns.gen_sticky_ends <- function(frags_intervals, intervals.set.out=NULL, file=NULL, frags_column='FID'){
    sticky_ends <- gintervals.load(frags_intervals) %>%
        group_by_(frags_column) %>%
        slice(1) %>%
        ungroup
    if (!is.null(intervals.set.out)){
        gintervals.save(intervals = sticky_ends, intervals.set.out = intervals.set.out)
    }
    if (!is.null(file)){
        fwrite(sticky_ends, file, sep=',')
    }
}



# Bissli Functions ------------------------------------------------

#' Align bisulfite converted reads to a reference genome using bowtie2
#'
#' @param r1_fastq a vector of FASTQ files to align (read1)
#' @param out_bam bam output file with aligned reads
#' @param genome_seq A directory holding the reference genome, one FASTA file per
#' chromosome
#' @param bissli2_idx The basename of the index to be searched. The index is created using
#' gpatterns.bissli2_build
#' @param r2_fastq a vector of FASTQ files to align (read2)
#' @param bissli2_bin path for bissli2-align
#' @param bowtie2 path for bowtie2
#' @param samtools path for samtools
#' @param maxins maximum fragment length (bowtie2 maxins argument)
#' @param threads number of threads to use
#' @param genome_type ct/ga/ct_ga. if 'ct' - Assume that the reads to be aligned underwent C->T conversion. If
#' paired-end reads are aligned, then assume read 1 underwent C->T
#' conversion while read 2 underwent G->A conversion. The --ct and --ga
#' options are mutually exclusive. if 'ga' - Assume that the reads to be aligned underwent G->A conversion. If
#' paired-end reads are aligned, then assume read 1 underwent G->A
#' conversion while read 2 underwent C->T conversion.
#' if 'ct_ga' reads would be aligned to both C->T and G->A, and the best match would
#' be chosen.
#'
#' @param tmp_dir Directory for storing temporary files
#' @param bissli2_params additional parameters to bissli2/bowtie2
#'
#' @return
#' @export
#'
#' @examples
gpatterns.bissli2 <- function(r1_fastq,
                              out_bam,
                              genome_seq,
                              bissli2_idx,
                              r2_fastq=NULL,
                              bissli2_bin=.gpatterns.bissli2_bin,
                              bowtie2='bowtie2',
                              samtools = 'samtools',
                              maxins=1000,
                              threads=10,
                              genome_type='ct',
                              tmp_dir = NULL,
                              bissli2_params = ''){
    tmp_dir <- tmp_dir %||% tempdir()
    r1_fastq <- paste(r1_fastq, collapse=',')
    if (genome_type == 'ct_ga'){
        genome_type <- 'ct --ga'
    }
    if (is.null(r2_fastq)){
        command <- qq('@{bissli2_bin} @{bissli2_params} --tmp-dir @{tmp_dir} --bowtie2 @{bowtie2} --@{genome_type} -g @{genome_seq} -x @{bissli2_idx} -U @{r1_fastq} --threads @{threads} | @{samtools} view -b -S -h -o @{out_bam} -')
    } else {
        r2_fastq <- paste(r2_fastq, collapse=',')
        command <- qq('@{bissli2_bin} @{bissli2_params} --tmp-dir @{tmp_dir} --bowtie2 @{bowtie2} --maxins @{maxins} --@{genome_type} -g @{genome_seq} -x @{bissli2_idx} -1 @{r1_fastq} -2 @{r2_fastq} --threads @{threads} | @{samtools} view -b -S -h -o @{out_bam} -')
    }
    system(command)
}

#' Create an bowtie2 index for a bisulfite converted genome
#'
#' @param reference a vector of FASTA file names holding the reference genome.
#' @param idx_base The basename of the index files that will be created.
#' @param bowtie2_build_bin The path of the bowtie2-build executable. [bowtie2-build]
#' @param bowtie2_options Any unknown options are passed transparently to bowtie2-build.
#' Note that bissli2-build can only handle FASTA files as input. Therefore the '-f' option is always passed to bowtie2-build
#' @param bissli2_build_bin binary of bissli2-build
#'
#' @return
#' @export
#'
#' @examples
gpatterns.bissli2_build <- function(reference,
                                    idx_base,
                                    bowtie2_options,
                                    bowtie2_build_bin='bowtie2-build',
                                    bissli2_build_bin=.gpatterns.bissli2_build_bin){
    reference <- paste(reference, collapse=',')
    command <- qq('@{bissli2_build_bin} --bowtie2-build @{bowtie2_build_bin} @{reference} @{idx_base} @{bowtie2_options}')
    system(command)
}

# import Functions ------------------------------------------------

.step_invoke <- function(step, steps, f, ...){
    if (step %in% steps){
        message(qq('doing @{step}'))
        return(f(...))
    } else {
        message(qq('skipping @{step}'))
    }
}




#' Create a track from tidy_cpgs files
#'
#'
#' @param tidy_cpgs tidy_cpgs data frame or a vector with directories of tidy_cpgs (use full path)
#' @param track name of the track to generate
#' @param description description of the track to generate
#' @param steps steps of the pipeline. Possible options are:
#' 'bind_tidy_cpgs', 'pileup', 'pat_freq', 'pat_cov'
#' @param overwrite overwrite existing tracks
#' @param cov_filt_cmd if numeric - maximal coverage for CpG. Else - command for filtering highly (or lowly) covered CpGs. string with the maximal coverage, where 'covs' can represent the command, e.g. 'max(500, quantile(covs, 0.95))'. 
#' @param dsn downsampling n. Leave NULL for no downsampling
#' @param pat_cov_lens lengthes of patterns to calculate pattern coverage track for
#' @param max_span maximal span to look for patterns (usually the maximal insert length)
#' @param pat_freq_len lengthes of patterns to calculate pattern frequency track
#' @param nbins number of genomic bins to separate the analysis.
#' @param groot root of misha genomic database to save the tracks
#' @param use_sge use sun grid engine for parallelization
#' @param max_jobs maximal number of jobs for sge parallelization
#' @param parallel parallelize using threads (number of threads is determined by gpatterns.set_parallel)
#'
#' @return
#' @export
#'
#' @examples
gpatterns.import_from_tidy_cpgs <- function(tidy_cpgs,
                                            track,
                                            description,
                                            steps = 'all',
                                            overwrite = TRUE,
                                            cov_filt_cmd = NULL,
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
    all_steps <- c('bind_tidy_cpgs', 'pileup', 'pat_freq', 'pat_cov')

    if (length(steps) == 1 && steps[1] == 'all'){
        steps <- all_steps
    }

    stopifnot(all(steps %in% c(all_steps)))

    if (.is_tidy_cpgs(tidy_cpgs)){
        genomic_bins <- gbin_intervals(intervals = gintervals.all(), nbins)
        .step_invoke(
            'bind_tidy_cpgs',
            steps,
            .gpatterns.tidy_cpgs_to_files,
            tidy_cpgs = tidy_cpgs,
            intervals = genomic_bins,
            track = track)
    } else {
        # bind tidy cpgs
        .step_invoke(
            'bind_tidy_cpgs',
            steps,
            .gpatterns.bind_tidy_cpgs,
            tidy_cpgs_dirs=tidy_cpgs,
            track = track)
    }

    # pileup
    .step_invoke(
        'pileup',
        steps,
        .gpatterns.pileup,
        track = track,
        description = description,
        overwrite = overwrite,
        genomic_bins = genomic_bins,
        dsn = dsn,
        use_sge = use_sge,
        max_jobs = max_jobs,
        parallel = parallel, 
        cov_filt_cmd = cov_filt_cmd)

    # pat_freq
    .step_invoke(
        'pat_freq',
        steps,
        .gpatterns.pat_freq,
        track = track,
        description = description,
        overwrite = overwrite,
        pat_freq_len = pat_freq_len,
        nbins = 350,
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
        overwrite = overwrite,
        pat_cov_lens = pat_cov_lens,
        max_span = max_span,
        use_sge = use_sge,
        max_jobs = max_jobs,
        parallel = parallel)

}

#' Create a track from bam files.
#'
#' Creates a track from bam files.
#'
#' @param bams character vector with path of bam files
#' @param workdir directory in which the files would be saved (please provide full path)
#' @param steps steps of the pipeline to do. Possible options are:
#' 'bam2tidy_cpgs', 'filter_dups', 'bind_tidy_cpgs', 'pileup', 'pat_freq', 'pat_cov'
#' @param paired_end bam files are paired end, with R1 and R2 interleaved
#' @param cgs_mask_file comma separated file with positions of cpgs to mask
#' (e.g. MSP1 sticky ends). Needs to have chrom and start fields with the
#' position of 'C' in the cpgs to mask
#' @param trim trim cpgs that are --trim bp from the beginning/end of the read
#' @param umi1_idx position of umi1 in index (0 based)
#' @param umi2_idx position of umi2 in index (0 based)
#' @param use_seq use UMI sequence (not only position) to filter duplicates
#' @param only_seq use only UMI sequence (without positions) to filter duplicates
#' @param frag_intervs intervals set of the fragments to change positions to.
#' @param maxdist maximal distance from fragments
#' @param rm_off_target if TRUE - remove reads with distance > maxdist from frag_intervs
#' if FALSE - those reads would be left unchanged
#' @param add_chr_prefix add "chr" prefix for chromosomes (in order to import to misha)
#' @param bismark bam was aligned using bismark
#' @param import_raw_tcpgs import raw tidy cpgs to misha (without filtering duplicates)
#' @param cmd_prefix prefix to run on 'system' commands (e.g. source ~/.bashrc)
#' @param run_per_interv split run of bam2tidy_cpgs scripts separatly for each interval.
#' @param ... gpatterns.import_from_tidy_cpgs parameters
#'
#' @inheritParams gpatterns.import_from_tidy_cpgs
#'
#' @return if 'stats' is one of the steps - data frame with statistics. Otherwise none.
#' @export
#'
#' @examples
gpatterns.import_from_bam <- function(bams,
                                      workdir = NULL,
                                      track = NULL,
                                      steps = 'all',
                                      paired_end = TRUE,
                                      cgs_mask_file = NULL,
                                      trim = NULL,
                                      umi1_idx = NULL,
                                      umi2_idx = NULL,
                                      use_seq = FALSE,
                                      only_seq = FALSE,
                                      frag_intervs = NULL,
                                      maxdist = 0,
                                      rm_off_target = TRUE,
                                      add_chr_prefix = FALSE,
                                      bismark = FALSE,
                                      nbins = nrow(gintervals.all()),
                                      groot = GROOT,
                                      import_raw_tcpgs = FALSE,
                                      use_sge = FALSE,
                                      max_jobs = 400,
                                      parallel = getOption('gpatterns.parallel'),
                                      cmd_prefix = '',
                                      run_per_interv = TRUE,
                                      ...){
    gsetroot(groot)
    all_steps <- c('bam2tidy_cpgs', 'filter_dups', 'bind_tidy_cpgs', 'pileup', 'pat_freq', 'pat_cov', 'stats')    

    if (length(steps) == 1 && steps == 'all'){
        steps <- all_steps
    }

    stopifnot(all(steps %in% c(all_steps)))

    if (is.null(workdir) && !is.null(track)){
        if ('bind_tidy_cpgs' %in% steps){
            workdir <- paste0(.gpatterns.base_dir(track), '/workdir')        
        } else {
            workdir <- .gpatterns.base_dir(track)         
        }        
    }

    if (is.null(workdir) && !is.null(track)){
        stop('need to supply either workdir or track')
    }

    workdir <- normalizePath(workdir)

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
        paired_end = paired_end,
        cgs_mask_file = cgs_mask_file,
        umi1_idx = umi1_idx,
        umi2_idx = umi2_idx,
        trim = trim,
        frag_intervs = frag_intervs,
        maxdist = maxdist,
        rm_off_target = rm_off_target,
        add_chr_prefix = add_chr_prefix,
        bismark = bismark,
        use_sge = use_sge,
        max_jobs = max_jobs,
        parallel = parallel,
        cmd_prefix = cmd_prefix,
        run_per_interv = run_per_interv)

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
        use_seq = use_seq,
        only_seq = only_seq,
        use_sge = use_sge,
        max_jobs = max_jobs,
        parallel = parallel,
        cmd_prefix = cmd_prefix)

    tidy_cpgs_steps <- c('bind_tidy_cpgs', 'pileup', 'pat_freq', 'pat_cov')

    tidy_cpgs_dirs <- if (import_raw_tcpgs) qq('@{workdir}/tidy_cpgs') else qq('@{workdir}/tidy_cpgs_uniq')

    if (any(steps %in% tidy_cpgs_steps)){
        gpatterns.import_from_tidy_cpgs(tidy_cpgs_dirs,
                                        track = track,
                                        steps = steps[steps %in% tidy_cpgs_steps],
                                        nbins = nbins,
                                        groot = groot,
                                        use_sge = use_sge,
                                        parallel = parallel,
                                        max_jobs = max_jobs,
                                        ...)
    }

    # .step_invoke(
    #     'stats',
    #     steps,
    #     gpatterns.get_pipeline_stats,
    #     track = track,
    #     tidy_cpgs_stats_dir = qq('@{workdir}/tidy_cpgs/stats'),
    #     uniq_tidy_cpgs_stats_dir = qq('@{workdir}/tidy_cpgs_uniq/stats'))
}

#' Separate a track by strands (for QC purposes)
#'
#' @inheritParams gpatterns::gpatterns.import_from_tidy_cpgs
#' @param out_track output track name if NULL \code{track} would be used
#' @param intervals intervals to extract from tidy_cpgs
#' @param minus_suffix suffix for the minus strand track
#' @param plus_suffix suffix for the plus strand track
#' @param ... additional parameters for \code{gpatterns.import_from_tidy_cpgs}
#'
#' @return none
#' @export
#'
#' @examples
gpatterns.separate_strands <- function(track, description, out_track=NULL, intervals=NULL, minus_suffix='.minus', plus_suffix='.plus', ...){
    out_track <- out_track %||% track
    tcpgs <- gpatterns.get_tidy_cpgs(track, intervals=intervals)
    tcpgs_plus <- tcpgs %>% filter(strand == '+')
    tcpgs_minus <- tcpgs %>% filter(strand == '-')

    message(qq('doing plus strand, writing to @{out_track}@{plus_suffix}'))
    gpatterns.import_from_tidy_cpgs(
        tcpgs_plus,
        track = qq('@{out_track}@{plus_suffix}'),
        description = qq('@{description} plus strand'),
        ...
    )
    message(qq('doing minus strand, writing to @{out_track}@{minus_suffix}'))
    gpatterns.import_from_tidy_cpgs(
        tcpgs_minus,
        track = qq('@{out_track}@{minus_suffix}'),
        description = qq('@{description} minus strand'),
        ...
    )
}




.gpatterns.bam2tidy_cpgs <- function(bams,
                                     tidy_cpgs_dir,
                                     stats_dir,
                                     genomic_bins,
                                     cgs_mask_file = NULL,
                                     trim = NULL,
                                     paired_end = TRUE,
                                     umi1_idx = NULL,
                                     umi2_idx = NULL,
                                     frag_intervs = NULL,
                                     maxdist = 0,
                                     rm_off_target = TRUE,
                                     sort_output = FALSE,
                                     only_seq = FALSE,
                                     bismark = FALSE,
                                     adjust_read_bin = .gpatterns.adjust_read_bin,
                                     bin = .gpatterns.bam2tidy_cpgs_bin,
                                     run_per_interv = TRUE,
                                     add_chr_prefix = FALSE,                                                    
                                     ...){
    walk(c(tidy_cpgs_dir, stats_dir), ~ system(qq('mkdir -p @{.x}')))
    bam_prefix <- if (1 == length(bams)) 'cat' else 'samtools cat'
    single_end <- if (!paired_end) '--single-end' else ''
    cgs_mask <- if (is.null(cgs_mask_file)) '' else qq('--cgs-mask @{cgs_mask_file}')
    trim_str <- if(is.null(trim)) '' else qq('--trim @{trim}')
    umi1_idx_str <- if (is.null(umi1_idx)) '' else qq('--umi1-idx @{umi1_idx}')
    umi2_idx_str <- if (is.null(umi2_idx)) '' else qq('--umi2-idx @{umi2_idx}')
    sort_fields <- if (only_seq) '6,7' else '2,7'
    rm_off_target_str <- if (rm_off_target) '--rm_off_target' else ''
    chr_prefix_str <- if(add_chr_prefix) '--add-chr-prefix' else ''
    bismark_str <- if(bismark) '--bismark' else ''

    if (!is.null(frag_intervs)){
        post_process_str <- qq(' | @{adjust_read_bin} @{rm_off_target_str} -f @{frag_intervs} --maxdist @{maxdist} --groot @{GROOT}')
    } else {
        post_process_str <- ''
    }
    if (sort_output){
        post_process_str <- qq('@{post_process_str} | awk \'NR==1; NR > 1 {print $0 | "sort --field-separator=, -k@{sort_fields} -k1 -k9"}\'')
    }

    if (run_per_interv){
        commands <- genomic_bins %>% purrrlyr::by_row( function(gbins){
            stats_fn <- qq('@{stats_dir}/@{gbins$chrom}_@{gbins$start}_@{gbins$end}.stats')
            output_fn <- qq('@{tidy_cpgs_dir}/@{gbins$chrom}_@{gbins$start}_@{gbins$end}.tcpgs.gz')
            qq('@{bam_prefix} @{paste(bams, collapse=\' \')} |
         @{bin} --no-progress -i - -o - -s @{stats_fn} @{bismark_str} @{umi1_idx_str} @{umi2_idx_str}
         --chrom @{gbins$chrom} --genomic-range @{gbins$start} @{gbins$end}
         @{chr_prefix_str} @{trim_str} @{single_end} @{cgs_mask} @{post_process_str} |
           gzip -c > @{output_fn}') %>%
                gsub('\n', '', .) %>% gsub('  ', ' ', .)
        }, .collate =  'cols', .to = 'cmd')

        .gpatterns.run_commands(commands, jobs_title = 'bam2tidy_cpgs', ...)
    } else {
        stats_fn <- qq('@{stats_dir}/all.stats')
        output_fn <- sprintf("%s.gz", tempfile())
        cmd <- qq('@{bam_prefix} @{paste(bams, collapse=\' \')} |
         @{bin} --no-progress -i - -o - -s @{stats_fn} @{bismark_str} @{umi1_idx_str} @{umi2_idx_str}
         @{chr_prefix_str} @{trim_str} @{single_end} @{cgs_mask} @{post_process_str} |
           gzip -c > @{output_fn}') %>%
            gsub('\n', '', .) %>% gsub('  ', ' ', .)

        system(cmd)
        tcpgs <- fread(qq('gzip -d -c @{output_fn}'), colClasses=.gpatterns.tcpgs_colClasses(uniq=FALSE)) %>% tbl_df
        .gpatterns.tidy_cpgs_to_files(tcpgs, genomic_bins, outdir=tidy_cpgs_dir)
    }

}

.gpatterns.filter_dups <- function(tidy_cpgs_dir,
                                   stats_dir,
                                   uniq_tidy_cpgs_dir,
                                   genomic_bins,
                                   paired_end = TRUE,
                                   use_seq = FALSE,
                                   only_seq = FALSE,
                                   sorted = FALSE,
                                   bin = .gpatterns.filter_dups_bin,
                                   ...) {

    walk(c(uniq_tidy_cpgs_dir, stats_dir), ~ system(qq('mkdir -p @{.x}')))
    all_files <- list.files(tidy_cpgs_dir, pattern='.tcpgs.gz', full.names = TRUE)
    files <- .gpatterns.intervals2files(genomic_bins, all_files, qq('@{tidy_cpgs_dir}/'))

    if (!all(all_files %in% files)){
        stop(sprintf('missing the following tidy cpgs files:\n\t%s',
                     paste(all_files[!(all_files %in% files)] ,collapse='\n\t')))
    }
    single_end <- if (!paired_end) '--only_R1' else ''
    use_seq <- if (use_seq) '--use-seq' else ''
    only_seq <- if (only_seq) '--only-seq' else ''
    sorted_str <- if (sorted) '--sorted' else ''

    commands <- genomic_bins %>% purrrlyr::by_row(function(gbins){
        stats_fn <- qq('@{stats_dir}/@{gbins$chrom}_@{gbins$start}_@{gbins$end}.stats')
        input_fn <- qq('@{tidy_cpgs_dir}/@{gbins$chrom}_@{gbins$start}_@{gbins$end}.tcpgs.gz')
        output_fn <- qq('@{uniq_tidy_cpgs_dir}/@{gbins$chrom}_@{gbins$start}_@{gbins$end}.tcpgs.gz')
        qq('@{bin} -i @{input_fn} -o - -s @{stats_fn} @{single_end} @{use_seq} @{only_seq} @{sorted_str} | gzip -c > @{output_fn}')
        }, .collate =  'cols', .to = 'cmd')
    .gpatterns.run_commands(commands, jobs_title = 'filter_dups', ...)

}


.gpatterns.join_filter_dups <- function(tidy_cpgs_dirs,
                                        tmp_dir = NULL,
                                        ...){
    tmp_dir <- tmp_dir %||% tempdir()
    genomic_bins <- .gpatterns.get_tidy_cpgs_intervals(path = tidy_cpgs_dirs[1])
    #TODO: add test that genomic bins are the same
    .gpatterns.bind_tidy_cpgs(tidy_cpgs_dirs = tidy_cpgs_dirs, path=tmp_dir)
    return(.gpatterns.filter_dups(tidy_cpgs_dir = file.path(tmp_dir, 'tidy_cpgs'), genomic_bins=genomic_bins, ...))
}


.gpatterns.bind_tidy_cpgs <- function(tidy_cpgs_dirs, track=NULL, path=NULL){
    .collate_gzips <- function(files, outfile) {
        if (length (files) == 1){
            return(system(qq('ln -sf @{files} @{outfile}')))
        } else {
            system(qq('cp @{files[1]} @{outfile}'))
            for (i in 2:length(files)){
                return(system(qq('cat @{files[i]} | gzip -d -c | tail -n +2 | gzip -c >> @{outfile}')))
            }
        }
    }

    stopifnot(!is.null(track) || !is.null(path))
    if (!is.null(track)){
        track_path <- .gpatterns.base_dir(track)
    } else {
        track_path <- path
    }

    system(qq('mkdir -p @{track_path}'))
    if (1 == length(tidy_cpgs_dirs)){
        system(qq('ln -sf @{tidy_cpgs_dirs} @{track_path}/tidy_cpgs'))
    } else {
        system(qq('mkdir -p @{track_path}/tidy_cpgs'))
        tidy_cpgs_dirs %>%
            map_df(~ tibble(lib=.x, fn=list.files(.x, pattern='.*\\.tcpgs\\.gz$'))) %>%
            mutate(fn1 = fn) %>%
            group_by(fn) %>%
            by_slice(~ .collate_gzips(paste0(.x$lib, '/', .x$fn1[1]),
                                      qq('@{track_path}/tidy_cpgs/@{.x$fn1[1]}')) )
    }
}


#' @export
.gpatterns.tidy_cpgs_to_files <- function(tidy_cpgs, intervals, track=NULL, outdir=NULL){
    if (!is.null(track)){
        outdir <- paste0(.gpatterns.base_dir(track), '/tidy_cpgs')
    }

    system(qq('mkdir -p @{outdir}'))
    tidy_cpgs <- tidy_cpgs %>%
        bind_cols(tidy_cpgs %>%
                      select(chrom, start=cg_pos) %>%
                      mutate(end = start + 1) %>%
                      gintervals.filter(intervals, bind_intervals2 = TRUE) %>%
                      select(chrom1, start1, end1)
                  )
    tidy_cpgs %>% unite('grp', chrom1, start1, end1, remove=F) %>% group_by(grp) %>% by_slice(
        function(x) {
        tmp <- tempfile()
        x %>%
            select(-chrom1, -start1, -end1) %>%
            mutate(meth = ifelse(meth == 1, 'Z', 'z')) %>%
            fwrite(tmp, sep=',', row.names=F, col.names=T)

        system(qq('cat @{tmp} | gzip -c > @{outdir}/@{x$chrom1[1]}_@{x$start1[1]}_@{x$end1[1]}.tcpgs.gz'))
        system('rm -f @{tmp}')
    } )
}


#' @export
.gpatterns.pileup <- function(track, description, dsn = NULL, columns = c('meth', 'unmeth', 'cov', 'avg'), overwrite=TRUE, cov_filt_cmd = NULL, ...){
    message('calculating pileup...')        
    pileup <- gpatterns.apply_tidy_cpgs(track, function(x) gpatterns.tidy_cpgs_2_pileup(x, dsn=dsn), ...) %>% ungroup

    if (!is.null(cov_filt_cmd)){
        message(qq('filtering using the following rule: cov <= @{cov_filt_cmd}'))        
        covs <- pileup[['cov']]
        max_cov <-  eval(parse(text=cov_filt_cmd) )
        pileup <- pileup %>% filter(cov <= max_cov)
    }
    message('importing pileup to misha...')    
    .gpatterns.import_intervs_table(track, description, pileup, columns=columns, overwrite = overwrite)
}

#'@export
.gpatterns.pat_freq <- function(track, description, pat_freq_len, nbins=NULL, split_by_bin = TRUE, overwrite=TRUE, ...){
    for (pat_freq_l in pat_freq_len){
        message(qq('calculating pattern frequency (pattern length: @{pat_freq_l})...'))

        if (!is.null(nbins)){
            intervals <- gbin_intervals(intervals = gintervals.all(), nbins)
        } else {
            intervals <- nbins
        }
        pat_freq <- gpatterns.apply_tidy_cpgs(
                track,
                function(x)
                    gpatterns.tidy_cpgs_2_pat_freq(x,
                                                   pat_length =
                                                       pat_freq_l,
                                                   tidy =
                                                       FALSE),
                intervals = intervals,
                split_by_bin = split_by_bin,
                ...)



        if (nrow(pat_freq) == 0){
            warning('no patterns were found. Skipping creation of pattern frequency tracks')
        } else {
            message('importing pat_freq to misha...')
            new_track <- qq('@{track}.pat@{pat_freq_l}')
            if (gtrack.exists(new_track)){
                if (overwrite){
                    gtrack.rm(new_track, force=TRUE)
                } else {
                    return(NULL)
                }
            }
            message(qq('creating @{new_track}'))
            gtrack.array.import_from_df(df = pat_freq, track=new_track, description=description)
        }
    }
}


.gpatterns.pat_cov <- function(track, description, pat_cov_lens, max_span = 500, overwrite=TRUE, ...){
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
    .gpatterns.import_intervs_table(track, description, pat_cov, columns = paste0('pat_cov', pat_cov_lens), overwrite = overwrite)
}

#' @export
.gpatterns.import_intervs_table <- function(track_pref, description, tab, columns = NULL, overwrite = FALSE){
    tab <- tbl_df(tab)
    tracks_dir <- .gpatterns.base_dir(track_pref)
    if (!dir.exists(tracks_dir)){
        system(qq('mkdir -p @{tracks_dir}'))
    }
    columns <- columns %||% colnames(tab)[!(colnames(tab) %in% c('chrom', 'start', 'end'))]
    walk(columns, function(col){
        track_name <- qq('@{track_pref}.@{col}')
        if (gtrack.exists(track_name)){
            if (overwrite){
                message(qq('removing @{track_name}'))
                gtrack.rm(track_name, force=TRUE)
            } else {
                warning(qq("did not overwrite track @{track_name}. To do that please set overwrite to TRUE"))
                return(NULL)
            }
        }
        message(qq('creating @{track_name}'))
        gtrack.create_sparse(track = track_name,
                             description = qq('@{description}: @{col} track'),
                             intervals = tab %>% select(chrom, start, end),
                             value = tab[[col]])
    })
}


.gpatterns.run_commands <- function(commands, use_sge=FALSE, max_jobs=400, parallel = getOption('gpatterns.parallel'), jobs_title='', cmd_prefix = '', ...){

      if (use_sge){
        command_list <- 1:length(commands$cmd) %>% map(~ qq('@{cmd_prefix} system(commands$cmd[@{.x}])'))
        res <- gcluster.run2(command_list = command_list,
                             max.jobs = max_jobs,
                             packages = 'gpatterns',
                             jobs_title = jobs_title,
                             collapse_results = FALSE, ....)
        codes <-  res %>% map_int(~ .x$retv)
    } else {
        res <- commands %>% plyr::alply(1, function(x)
            system(paste(qq('@{cmd_prefix} '), x$cmd)), .parallel=parallel)
        codes <- res
    }
    walk2(codes, commands$cmd, function(code, cmd) if (code != 0) stop(qq('command "@{cmd}" failed')))
}



# Patterns import Functions ------------------------------------------------


#' Creates patterns attributes of track
#'
#' @param track track name
#' @param description a character string description
#' @param pat_space pattern space data frame
#' (output of gpatterns.tracks_to_pat_space or gpatterns.intervs_to_pat_space)
#' needs to have the following fields: chrom,start,end,fid
#' @param patterns_tab table with patterns (output of gpatterns.tidy_cpgs_to_pats)
#' needs to have the following fields: fid,pattern. If NULL it would be generated
#' automatically.
#' @param add_read_id save the read_id together with the patterns
#' @param noise_threshold threshold to consider pattern as 'noise'
#' @param overwrite overwrite existing tracks
#' @param canonize convert to canonic form (see: \code{\link[misha]{gintervals.canonic}}).
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
    if (!dir.exists(.gpatterns.base_dir(track))){
        dir.create(.gpatterns.base_dir(track), showWarnings=FALSE, recursive=TRUE)    
    }    

    # save pattern_space
    if (gintervals.exists(.gpatterns.pat_space_intervs_name(track))){
        if (overwrite){
            gintervals.rm(.gpatterns.pat_space_intervs_name(track), force=TRUE)
        } else {
            stop('pattern space intervals exist. Please run with overwrite=TRUE')
        }
    }

    message('saving pattern space...')
    gintervals.save(intervals=as.data.frame(pat_space), intervals.set.out=.gpatterns.pat_space_intervs_name(track))

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
    if (!dir.exists(.gpatterns.base_dir(track))){
        dir.create(.gpatterns.base_dir(track_ds), showWarnings=FALSE, recursive=TRUE)
    }

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
                                  overwrite = overwrite,
                                  save_tab = TRUE,
                                  ...)
    }
}

# Import utils ------------------------------------------------


#' Change tidy_cpgs coordinates to fragments coordinates
#'
#' @description changes the coordinates of tidy_cpgs to the fragment coordinates.
#' This is necessary to avoid over-counting of reads with coordinates that are
#' shifted by several bps from the restriction site.
#'
#' @param calls tidy_cpgs data frame
#' @param frag_intervs intervals set of the fragment
#' @param maxdist maximal distance from fragments
#' @param rm_off_target if TRUE - remove reads with distance > maxdist from frag_intervs
#' if FALSE - those reads would be left unchanged
#'
#' @return tidy_cpgs data frame with chrom,start,end coordinates adjusted to fragment intervals
#' @export
#'
#' @examples
gpatterns.adjust_read_pos <- function(calls, frag_intervs, maxdist=0, rm_off_target=TRUE){
    calls <- calls %>%
        .gpatterns.force_chromosomes()
    new_coords <- calls %>%
        select(chrom, start, end) %>%
        mutate(end = ifelse(end == '-', start + 1, as.numeric(end))) %>%
        mutate(end1=end, start1=start, start=ifelse(end>start, start, end), end=ifelse(end>start, end, start)) %>%
        select(-end1, -start1) %>%
        mutate(end = ifelse(end == start, start+1, end))

    capture.output(new_coords <-  new_coords %>%
        gintervals.neighbors1(frag_intervs))
    on_target_lines <- abs(new_coords$dist) <= maxdist

    if (rm_off_target){
        new_coords <- new_coords %>% filter(on_target_lines)
        calls <- calls %>% filter(on_target_lines) %>% mutate(start = new_coords$start1, end = new_coords$end1)
    } else {
        calls <- calls %>% mutate(start = ifelse(on_target_lines, new_coords$start1, start), end = ifelse(on_target_lines, new_coords$end1, end))
    }

    return(calls)
}

#' Merge tracks
#'
#' @description Creates a track called \code{new_track} that is the sum of the methylation calls in \code{tracks}.
#' In addition, if \code{merge_tidy_cpgs} is TRUE, creates a union of tidy_cpgs of \code{tracks} and generates
#' pattern frequency and pattern coverage for the combined tidy_cpgs.
#' If patterns attributes are present, they will be combined and calculated for the new track.
#'
#' @param tracks tracks to merge
#' @param new_track new track name
#' @param description new track description
#' @param intervals intervals scope. default is all the genome (gintrevals.all())
#' @param iterator new track iterator. if left NULL, the iterator would be the union of all the covered CpGs of \code{tracks}
#' @param add_var calculate variance (of tracks '.avg'). Will be created in new_track.var
#' @param merge_tidy_cpgs merge also tidy cpgs and create pattern frequency and pattern coverage track.
#' Does not work for tracks without the single molecule data (e.g. arrays, processed data)
#' @param ... additional parameters for \code{\link[gpatterns]{gpatterns.import_from_tidy_cpgs}}
#'
#' @return name of the new_track
#' @export
#'
#' @examples
gpatterns.merge_tracks <- function(tracks, new_track, description, intervals=gintervals.all(), iterator=NULL, add_var=FALSE, merge_tidy_cpgs = TRUE, ...){
    if (is.null(iterator)){
        if (!all(gtrack.exists(qqv('@{tracks}.cov')))){
            stop('not all tracks exist')
        }
        message('calculating covered CpGs')
        expr <- paste(qqv('!is.na(@{tracks}.cov)'), collapse=' | ')
        iterator <- gscreen(expr, intervals=intervals, iterator=.gpatterns.genome_cpgs_intervals)
    }

    gdir.create(gsub('\\.', '/', new_track), showWarnings=FALSE)

    for (suffix in c('cov', 'meth', 'unmeth')){
        message(qq('doing @{suffix}'))
        if (!all(gtrack.exists(qqv('@{tracks}.@{suffix}')))){
            stop('not all tracks exist')
        }
        expr <- sprintf("sum(%s, na.rm=T)", paste(qqv('@{tracks}.@{suffix}'), collapse=', '))
        if (!gtrack.exists(qq('@{new_track}.@{suffix}'))){
            gtrack.create(qq('@{new_track}.@{suffix}'), description, expr, iterator=iterator)
        }
    }

    if (add_var){
        expr <-  sprintf("var(%s, na.rm=T)", paste(qqv('@{tracks}.avg'), collapse=', '))
        if (!gtrack.exists(qq('@{new_track}.var'))){
            gtrack.create(qq('@{new_track}.var'), description, expr, iterator=iterator)
        }
    }
    gdb.reload()

    message(qq('doing avg'))
    if (!gtrack.exists(qq('@{new_track}.avg'))){
        gtrack.create(qq('@{new_track}.avg'), description, qq('@{new_track}.meth / @{new_track}.cov'), iterator=qq('@{new_track}.meth'))    
    }
    

    if (merge_tidy_cpgs){
        tidy_dirs <- .gpatterns.tidy_cpgs_dir(tracks)
        do.call_ellipsis(gpatterns.import_from_tidy_cpgs, list(tidy_cpgs=tidy_dirs, track=new_track, description=description, steps=c('bind_tidy_cpgs', 'pat_freq', 'pat_cov')), ...)
    }

    if (all(.gpatterns.patterns_exist(tracks))){
         pat_space <- reduce(.gpatterns.pat_space_intervs_name(tracks), function(x, y) gintervals.load(x) %>% full_join(gintervals.load(y), by=c('chrom', 'start', 'end', 'fid'))) %>% tbl_df
         gpatterns.create_patterns_track(new_track, description, pat_space=pat_space, ...)
    }

    return(new_track)
}

#' Get fragment coverage
#'
#' @description get number of molecules per each fragment for a track. This can
#' give us some upper bound of the number of methylation patterns per fragment.
#'
#' @param track track
#' @param intervals fragment intervals. if NULL - fragments would be based chrom,start,end positions
#' of tidy_cpgs, assuming gpatterns.adjust_read_pos was called before (as in any
#' run of gpatterns.import_* with frag_intervs parameter).
#' @save_track save the result to \code{track}.fid_cov track
#' @param parallel parallel
#' @param ... other parameters of gpatterns.adjust_read_pos
#'
#' @return data frame with fragmets coordinates (chrom,start,end) and 'frag_cov'
#' field with the fragment coverage
#'
#' @export
gpatterns.get_fragment_cov <- function(track,
                                       intervals=NULL,
                                       save_track = FALSE,
                                       parallel=getOption('gpatterns.parallel'),
                                       ...){
    get_frag_cov <- function(tcpgs, intervals=NULL){
        if (is.null(intervals)){
            return(tcpgs %>% distinct(read_id, chrom, start, end) %>% count(chrom, start, end))
        } else {
            return(tcpgs %>% gpatterns.adjust_read_pos(frag_intervs = intervals, ...) %>% distinct(read_id, chrom, start, end) %>% count(chrom, start, end))
        }
    }
    res <- gpatterns.apply_tidy_cpgs(track, function(x) get_frag_cov(x), parallel=parallel) %>%  group_by(chrom, start, end) %>% summarise(n = sum(n, na.rm=T)) %>% rename(frag_cov=n)
    if (save_track){
        gtrack.create_sparse(paste0(track, '.frag_cov'), description = description, intervals = res %>% select(chrom, start, end), values = res$frag_cov)
    }
    return(res)
}

