
# gpatterns.track_stats <- function(track){
#     stats <- list()
#     sm <- gsummary(qq('@{track}.cov'))
#     stats[['cg_num']] <- sm[1]
#     stats[['meth_calls']] <- sm[5]
#     stats[['global_avg_meth']] <- gsummary(qq('@{track}.avg'))[6]
#     stats[['reads_per_umi']] <- .gpatterns.reads_per_umi(track)
#     return(stats)
# }

#' Calculate QC statistics for track
#'
#' @param track track name (or mutiple track names)
#' @param tidy_cpgs_stats_dir directory with tidy_cpgs stats files (cannot be used with multiple tracks)
#' @param uniq_tidy_cpgs_stats_dir directory with filter_dups stats files (cannot be used with multiple tracks)
#' @param add_mapping_stats add full mapping statistics (singles, discordant etc.)
#' @param add_insert_length_stats add average insert length
#'
#' @return data frame with 'track' field with the track name pipeline statistics 
#' @export
#'
#' @examples
gpatterns.get_pipeline_stats <- function(track,
                                        tidy_cpgs_stats_dir = NULL,
                                        uniq_tidy_cpgs_stats_dir = NULL,
                                        add_mapping_stats = FALSE,
                                        add_insert_length_stats = TRUE){
    if (length(track) > 1){
        if (!is.null(tidy_cpgs_stats_dir) || !is.null(uniq_tidy_cpgs_stats_dir)){
            stop('cannot use explicit tidy_cpgs_stats_dir or uniq_tidy_cpgs_stats_dir with multiple tracks')
        }
        return(map_df(track, gpatterns.get_pipeline_stats))
    }
    if (is.null(tidy_cpgs_stats_dir)){
        tidy_cpgs_stats_dir <- paste0(GROOT, '/tracks/', gsub('\\.', '/', track), '/workdir/tidy_cpgs/stats') 
    }
    if (is.null(uniq_tidy_cpgs_stats_dir)){
        uniq_tidy_cpgs_stats_dir <- paste0(GROOT, '/tracks/', gsub('\\.', '/', track), '/workdir/tidy_cpgs_uniq/stats')
    }
    stopifnot(dir.exists(tidy_cpgs_stats_dir))
    stopifnot(dir.exists(uniq_tidy_cpgs_stats_dir))

    s <- gpatterns.get_tcpgs_stats(tidy_cpgs_stats_dir, uniq_tidy_cpgs_stats_dir)
    stats <- s$stats
        
    sm <- gsummary(qq('@{track}.cov'))
    stats[['cg_num']] <- sm[1]
    stats[['meth_calls']] <- sm[5]
    stats[['global_avg_meth']] <- gsummary(qq('@{track}.avg'))[6]

    if (add_insert_length_stats){
        stats <- stats %>%
            bind_cols(gpatterns.apply_tidy_cpgs(track,
                    function(x) x %>%
                        summarise(insert_len = mean(abs(insert_len), na.rm=T), n = n())) %>%
                        mutate(f = insert_len * n / sum(n)) %>% 
                        summarise(insert_len = sum(f)))    
    }
    

    if (add_mapping_stats){        
        stats <- stats %>% bind_cols(s$mapping_stats)
    }

    stats <- stats %>% mutate(track = track) %>% select(track, everything()) %>% mutate(meth_call_per_read = meth_calls / total_reads)    

    return(stats)
}

#' Get QC statistics for track
#' @param track track name (or mutiple track names)
#' @inheritDotParams gpatterns.get_pipeline_stats
#' @export
gpatterns.track_stats <- function(track, ...){
    if (length(track) > 1){
        return(map_df(track, gpatterns.track_stats))        
    }
    stats_fn <- .gpatterns.stats_file_name(track)
    if (file.exists(stats_fn)){
        return(fread(stats_fn) %>% as.tibble())
    } else {
        stats <- gpatterns.get_pipeline_stats(track, ...)
        write_tsv(stats, stats_fn)
    }
    return(stats)
}

#' Get QC statistics from tcpgs dir
#' 
#' @inheritParams gpatterns.get_pipeline_stats
#' 
#' @export
gpatterns.get_tcpgs_stats <- function(tidy_cpgs_stats_dir, uniq_tidy_cpgs_stats_dir){ 
    stopifnot(dir.exists(tidy_cpgs_stats_dir))
    stopifnot(dir.exists(uniq_tidy_cpgs_stats_dir))
    
    uniq_stats <- list.files(uniq_tidy_cpgs_stats_dir, full.names=TRUE) %>%
        map_df(~ fread(.x)) %>%
        summarise(total_reads = sum(total_reads, na.rm=TRUE), uniq_reads = sum(uniq_reads, na.rm=TRUE)) %>%
        mutate(uniq_frac = uniq_reads / total_reads)
    tidy_cpgs_stats <- list.files(tidy_cpgs_stats_dir, full.names=TRUE) %>%
        map_df(~ fread(.x)) %>% mutate(CHH = mean(CHH, na.rm=TRUE), CpG = mean(CpG, na.rm=TRUE)) %>%
        slice(1)

    tidy_cpgs_stats <- tidy_cpgs_stats %>% select(-one_of("reads_with_cg"))
    # if ("reads_with_cg" %in% colnames(tidy_cpgs_stats)) {
    #     tidy_cpgs_stats$reads_with_cg <-  list.files(tidy_cpgs_stats_dir, full.names=TRUE) %>% map_df(~ fread(.x)) %>% pull(reads_with_cg) %>% sum(na.rm=TRUE)        
    # }
    for (f in c('good', 'single_R1', 'single_R2', 'bad_cigar', 'no_conv', 'unmapped', 'discordant')){
        if (!(f %in% colnames(tidy_cpgs_stats))){
            tidy_cpgs_stats[[f]] <- 0
        }
    }    
    
    stats <- tidy_cpgs_stats %>%
        mutate(total_reads = good + single_R1 + single_R2 + bad_cigar + no_conv + unmapped + discordant,
                  mapped_reads = good + single_R1 + single_R2,
                  mapped_frac = mapped_reads / total_reads)
    
    stats <- bind_cols(stats, uniq_stats %>% rename(reads_with_cg = total_reads, uniq_reads_with_cg = uniq_reads)) %>%
        select(one_of('total_reads', 'mapped_reads', 'mapped_frac', 'reads_with_cg', 'uniq_reads_with_cg', 'uniq_frac', 'CHH', 'CHG', 'CpG')) 
    
    return(list(stats=stats, mapping_stats=tidy_cpgs_stats))
}

# total reads     mapped.reads    mapped.uniq    perc mapped % unique overall     on target  perc on target   uniq_1      uniq_2     perc unique on target   perc unique off target   estimated on target complexity      regions not present     total methylation calls     CpGs per read   num of CpGs     average CpG cov 
#' Capture stats
#' @export
gpatterns.capture_stats <- function(track, 
                                    bam, 
                                    regions, 
                                    fastq=NULL, 
                                    use_read2=TRUE, 
                                    max_distance=200, 
                                    tidy_cpgs_stats_dir = NULL,
                                    uniq_tidy_cpgs_stats_dir = NULL,
                                    add_mapping_stats = FALSE, 
                                    dsn=NULL, 
                                    sample_raw = FALSE,
                                    sample_ontar = TRUE, 
                                    reads = NULL){    
    if (is.null(reads)){
        message('extrating reads...')
        reads <- .gpatterns.bam2reads(bams=bam)                
    } else {
        if (is.character(reads)){
            reads <- fread(reads) %>% as.tibble()
        }
    }

    if (is.character(regions)){        
        regions <- gintervals.load(regions)
    }

    if (!is.null(dsn)){
        if (sample_raw){
            message(sprintf('sampling %d raw reads...', dsn))
            bam_prefix <- if (1 == length(bam)) 'cat' else 'samtools cat'
            cmd <- qq('@{bam_prefix} @{paste0(bam, collapse=" ")} | samtools view - | sed -n 1~2p | cut -f1 | head -n @{dsn}')
            ids <- fread(cmd, 
                sep='\t', header=F)[,1]
            if (length(ids) != dsn){
                message(sprintf('skipping  - not enough reads (%d) for downsampling %d raw reads', nrow(reads), dsn))
                return(tibble(on_target=NA, frac_ontar=NA, frac_uniq_ontar=NA, frac_uniq_offtar=NA))
            }
            reads <- reads %>% filter(read_id %in% ids)

        } else {
            if (nrow(reads) < dsn){
                message(sprintf('skipping  - not enough reads (%d) for downsampling %d reads', nrow(reads), dsn))
                return(tibble(on_target=NA, frac_ontar=NA, frac_uniq_ontar=NA, frac_uniq_offtar=NA))
            }
            if (!sample_ontar){
                message(sprintf('sampling %d mapped reads...', dsn))
                reads <- reads %>% sample_n(dsn)    
            } 
        }        
    }

    message('filtering duplicates for all reads...')
    f_reads <- .gpatterns.filter_read_dups(reads, use_read2=use_read2)  
    conv_stats <- f_reads %>% summarise(h = sum(h), H = sum(H), x = sum(x), X = sum(X)) %>% mutate(CHH = H / (H + h), CHG = X / (x + X)) %>% select(CHH, CHG)

    message('identifying \'on target\' reads...')
    on_tar_reads <- .gpatterns.filter_ontar_reads(reads, regions, max_distance=max_distance) %>% as.tibble() 
    ontar_poss <- .gpatterns.filter_ontar_reads(reads, regions, max_distance=max_distance, return_orig=FALSE) %>% as.tibble() 
    regsion_not_present <- sum(regions %>% gintervals.neighbors1(ontar_poss) %>% pull(dist) %>% abs(.) >= max_distance)

    if (!is.null(dsn) && !sample_raw && sample_ontar){    
        if (nrow(on_tar_reads) >= dsn){
            message(sprintf('sampling %d \'on target\' reads...', dsn))
            on_tar_reads <- on_tar_reads %>% sample_n(dsn)  
        } else {
            message(sprintf('skipping  - not enough reads (%d) for downsampling %d on target reads', nrow(on_tar_reads), dsn))
            return(tibble(on_target=NA, frac_ontar=NA, frac_uniq_ontar=NA, frac_uniq_offtar=NA))
        }        
    }

    message('filtering duplicates for on traget...')    
    f_reads_ontar <- .gpatterns.filter_read_dups(on_tar_reads, use_read2=use_read2)    
    
    message('filtering duplicates for off traget...')
    offtar_reads <- reads %>% filter(!(read_id %in% on_tar_reads$read_id))
    f_reads_offtar <- .gpatterns.filter_read_dups(offtar_reads, use_read2=use_read2)
    
    message('getting other stats...')        
    stats <- gpatterns.get_pipeline_stats(track, tidy_cpgs_stats_dir=tidy_cpgs_stats_dir, uniq_tidy_cpgs_stats_dir=uniq_tidy_cpgs_stats_dir, add_mapping_stats=add_mapping_stats)
    stats[['on_target']] <- nrow(on_tar_reads)
    stats[['frac_ontar']] <- nrow(on_tar_reads) / nrow(reads)
    stats[['frac_uniq_ontar']] <- nrow(f_reads_ontar) / nrow(on_tar_reads)
    stats[['frac_uniq_offtar']] <- nrow(f_reads_offtar) / nrow(offtar_reads)
    stats[['regsion_not_present']] <- regsion_not_present

    stats <- stats %>% select(-one_of('CHH', 'CHG'))

    stats <- bind_cols(stats, conv_stats)

    if (!is.null(dsn)){
        stats <- stats %>% mutate(dsn = dsn)
    }
    
    return(stats)
}

.gpatterns.insert_length_reads_per_umi <- function(track, regions, dsn=NULL, reads_umi_breaks=c(0,1,4,10,500), break_levels = c('1' = '[0,1]', '2-4' = '(1,4]', '5-10' = '(4,10]', '>10' = '(10,500]')){
    tcpgs <-  gpatterns.get_tidy_cpgs(track)    
    tcpgs <- tcpgs %>% distinct(read_id, num, .keep_all=T) %>% filter(end != '-')

    message('identifying \'on target\' reads...')
    tcpgs_ontar <- .gpatterns.filter_ontar_reads(tcpgs, regions)
    res <- tcpgs_ontar %>%
        mutate(insert_len = abs(insert_len), 
               num = cut(num, reads_umi_breaks, include.lowest=T), 
               num = do.call(forcats::fct_recode, c(list(f = num), as.list(break_levels) ))) %>% 
        select(insert_len, reads=num)
    return(res)
}

#' filter on target reads
#' Filters on target reads
#' @export
.gpatterns.filter_ontar_reads <- function(reads, regions, max_distance=200, return_orig=TRUE){
    poss <- reads %>% mutate(end = as.integer(ifelse(end == '-', start + 1, end)), new_start = pmin(start, end), new_end=pmax(start, end), new_end = ifelse(new_end == new_start, new_start + 1, new_end)) %>% select(chrom, start=new_start, end=new_end, read_id) %>% .gpatterns.force_chromosomes() %>% gintervals.filter(regions, max_distance=max_distance)
    if (return_orig){
        return(reads %>% filter(read_id %in% poss$read_id))
    } else {
        return(poss)
    }    
}

.gpatterns.reads_per_region <- function(regions, bams=NULL, reads=NULL, max_distance=200){
    stopifnot(!is.null(bams) || !is.null(reads))
    if (is.null(reads)){
        message('extrating reads...')
        reads <- .gpatterns.bam2reads(bams=bams)        
    }
    reads <- .gpatterns.filter_read_dups(reads)
    reads <- reads %>% .gpatterns.filter_ontar_reads(regions, max_distance, return_orig=FALSE)   
    
    covs <- reads %>% 
        gintervals.neighbors1(regions) %>% 
        group_by(chrom1, start1, end1) %>% 
        summarise(n = n()) %>% 
        select(chrom=chrom1, start=start1, end=end1, cov=n)

    covs <- covs %>%
        full_join(regions) %>% 
        replace_na(replace=list(cov = 0))        

    return(covs)
}

#' reads per umi
#'
#' @param track
#' @param dsn
#'
#' @return
#'
#' @examples
.gpatterns.reads_per_umi <- function(track, dsn=NULL, regions=NULL){
    tcpgs <- gpatterns.get_tidy_cpgs(track) %>% distinct(read_id, num, .keep_all = TRUE)
    if (!is.null(regions)){        
        tcpgs <- .gpatterns.filter_ontar_reads(tcpgs, regions)
    }
    if (!is.null(dsn)){
        tcpgs_all <- tcpgs %>% untable('num') %>% select(-num)           
        tcpgs <- tcpgs_all %>% sample_n(dsn) %>% group_by(read_id) %>% summarise(num = n())
    }
    return(tcpgs %>% group_by(reads=num) %>% summarise(n=n()))
}


#dsns = c(2e5, 5e5, 7e5, 1e6, 1.5e6, 2e6, 3e6)

#' ds stats
#' downsampling stats
#'
#' @param track
#' @param dsns
#' @param bam
#' @param fastq
#' @param intervals
#' @param sort_rand
#'
#' @export
#' @return
#'
#' @examples
.gpatterns.ds_stats <- function(track, dsns, bam = NULL, fastq = NULL, tidy_cpgs_dir = NULL, intervals=NULL, sort_rand=FALSE){
    stats <- list()
    stats[['cg_num']] <- tibble(dsn=numeric(), cg_num=numeric())
    stats[['meth_calls']] <- tibble(dsn=numeric(), meth_calls=numeric())
    stats[['global_avg_meth']] <- tibble(dsn=numeric(), avg_meth=numeric())
    stats[['reads_per_umi']] <-  tibble(dsn=numeric(), reads=numeric(), n=numeric())
    stats[['umis']] <-  tibble(dsn=numeric(), umis=numeric())
    stats[['on_tar']] <-  tibble(dsn=numeric(), on_tar=numeric())

    sort_rand_str <- if(sort_rand) ' sort -R | ' else ''

    for (dsn in dsns){
        message(qq('doing @{scales::comma(dsn)}'))
        if (!is.null(tidy_cpgs_dir)){                  
            non_uniq_tcpgs <- .gpatterns.get_tidy_cpgs_from_dir(tidy_cpgs_dir, intervals=intervals, uniq=FALSE)
        }
        tcpgs <- gpatterns.get_tidy_cpgs(track, intervals=intervals, only_tcpgs=TRUE)
 
        all_tcpgs <- gpatterns.get_tidy_cpgs(track, only_tcpgs=TRUE)

        ids <- -1
        if (!is.null(fastq)){
            message('sampling from fastq')
            cmd <- qq('gzip -d -c @{paste(fastq, collapse=" ")} | sed -n 1~4p | @{sort_rand_str} head -n @{dsn}')
            ids <- fread(cmd, sep='\t', header=F)[,1] %>%
                gsub('^@', '', .) %>%
                str_split(' ') %>% 
                map_chr(~ .x[1])                
        } else if (!is.null(bam)){
            message('sampling from bam')
            cmd <- qq('samtools view @{bam} | sed -n 1~2p | cut -f1 | @{sort_rand_str} head -n @{dsn}')
            ids <- fread(cmd, 
                sep='\t', header=F)[,1]
        } else if (!is.null(tidy_cpgs_dir)){
            message('sampling from non unique tidy cpgs')            
            if (length(unique(non_uniq_tcpgs$read_id)) >= dsn){
                ids <- non_uniq_tcpgs %>% distinct(read_id) %>% sample_n(dsn) %>% .$read_id
            }                        
        } else {
            message('sampling from tidy cpgs')
            if (length(unique(tcpgs$read_id)) >= dsn){
                ids <- tcpgs %>% distinct(read_id, num) %>% sample_n(dsn) %>% .$read_id
            }            
        }

        if (length(ids) == dsn){            
            tcpgs <- tcpgs %>% filter(read_id %in% ids)
            all_tcpgs <- all_tcpgs %>% filter(read_id %in% ids)

            stats[['cg_num']] <- stats[['cg_num']] %>% bind_rows(tibble(dsn=dsn, cg_num = tcpgs %>% distinct(chrom, cg_pos) %>% nrow))
            stats[['meth_calls']] <- stats[['meth_calls']] %>% bind_rows(tibble(dsn=dsn, meth_calls = nrow(tcpgs)))
            stats[['global_avg_meth']] <- stats[['global_avg_meth']] %>% bind_rows(tibble(dsn=dsn, avg_meth = mean(tcpgs$meth)))
            stats[['reads_per_umi']] <- stats[['reads_per_umi']] %>% bind_rows(tcpgs %>% group_by(reads=num) %>% summarise(n=n()) %>% mutate(dsn = dsn))
            stats[['umis']] <- stats[['umis']] %>% bind_rows(tibble(dsn=dsn, umis = length(unique(tcpgs$read_id))))
            stats[['on_tar']] <- stats[['on_tar']] %>% bind_rows(tibble(dsn=dsn, on_tar = length(unique(tcpgs$read_id)) / length(unique(all_tcpgs$read_id))))
        } else {
            message(qq('skipping @{dsn} (not enough reads)'))
        }
    }

    stats <- stats %>% map(~ .x %>% mutate(track = track))

    return(stats)
}


#' bam 2 reads
#' @export
.gpatterns.bam2reads <- function(bams, paired_end=TRUE, umi1_idx=NULL, umi2_idx=NULL, add_chr_prefix=FALSE, reads_fn=NULL, bismark = FALSE, bin =  system.file("import", "tidy_cpgs.py", package="gpatterns")){
    return_reads <- FALSE

    if (is.null(reads_fn)){
        reads_fn <- tempfile()
        on.exit(system('rm -f @{reads_fn}'))
    }

    bam_prefix <- if (1 == length(bams)) 'cat' else 'samtools cat'
    single_end <- if (!paired_end) '--single-end' else ''        
    umi1_idx_str <- if (is.null(umi1_idx)) '' else qq('--umi1-idx @{umi1_idx}')
    umi2_idx_str <- if (is.null(umi2_idx)) '' else qq('--umi2-idx @{umi2_idx}')    
    
    chr_prefix_str <- if(add_chr_prefix) '--add-chr-prefix' else ''
    bismark_str <- if(bismark) '--bismark' else ''

    cmd <- qq('@{bam_prefix} @{paste(bams, collapse=\' \')} |
         @{bin} -i - -o /dev/null -s /dev/null -r @{reads_fn} @{bismark_str} @{umi1_idx_str} @{umi2_idx_str}
         @{chr_prefix_str} @{single_end}') %>%
            gsub('\n', '', .) %>% gsub('  ', ' ', .)
    
    system(cmd)    

    invisible(fread(reads_fn))    
}

#' @export
gpatterns.bam2uniq_reads <- function(bams, paired_end=TRUE, umi1_idx=NULL, umi2_idx=NULL, add_chr_prefix=FALSE, reads_fn=NULL, bismark = FALSE, use_read2=TRUE){
    reads <- .gpatterns.bam2reads(bams=bams, paired_end=paired_end, umi1_idx=umi1_idx, umi2_idx=umi2_idx, add_chr_prefix=add_chr_prefix, reads_fn=reads_fn, bismark=bismark)
    uniq_reads <- .gpatterns.filter_read_dups(reads, use_read2=use_read2)
    return(uniq_reads)
}

.gpatterns.filter_read_dups <- function(reads, use_read2=TRUE){
    if (use_read2){
        f_reads <- reads %>% 
            arrange(chrom, start, strand, desc(end)) %>% 
            distinct(chrom, start, end, strand, .keep_all=T) %>% 
            group_by(chrom, start, strand) %>% 
            mutate(n_double = sum(end != '-')) %>% 
            ungroup %>% 
            filter(n_double == 0 | end != '-') %>% 
            select(-n_double)    
            
    } else {
        f_reads <- reads %>% 
            distinct(chrom, start, strand, .keep_all=TRUE)   
            
    }        
    return(f_reads)
}

#' converstion stats
#' @export
gpatterns.conv_stats <- function(bams, paired_end=TRUE, bismark = FALSE, min_mapq=30, bin =  system.file("import", "tidy_cpgs.py", package="gpatterns")){
    bam_prefix <- if (1 == length(bams)) 'cat' else 'samtools cat'
    single_end <- if (!paired_end) '--single-end' else ''         

    bismark_str <- if(bismark) '--bismark' else ''

    stats_fn <- tempfile()

    cmd <- qq('@{bam_prefix} @{paste(bams, collapse=\' \')} |
         @{bin} -i - -o /dev/null -s @{stats_fn} @{bismark_str} @{single_end} --min_map @{min_mapq}') %>%
            gsub('\n', '', .) %>%
            gsub('  ', ' ', .)
    system(cmd)
    stats <- fread(stats_fn)
    return(stats)
}
