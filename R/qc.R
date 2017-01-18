
gpatterns.track_stats <- function(track){
    stats <- list()
    sm <- gsummary(qq('@{track}.cov'))
    stats[['cg_num']] <- sm[1]
    stats[['meth_calls']] <- sm[5]
    stats[['global_avg_meth']] <- gsummary(qq('@{track}.avg'))[6]
    stats[['reads_per_umi']] <- .gpatterns.reads_per_umi(track)
    return(stats)
}


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
gpatterns.get_pipeline_stats <- function(track,
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

    stats <- stats %>%
        bind_cols(
            gpatterns.apply_tidy_cpgs(track,
                                      function(x) x %>%
                                          summarise(insert_len = mean(abs(insert_len), na.rm=T),
                                                    n = n())) %>%
                mutate(f = insert_len * n / sum(n)) %>% summarise(insert_len = sum(f)))

    if (add_mapping_stats){
        stats <- stats %>% bind_cols(tidy_cpgs_stats)
    }

    return(stats)
}


#' reads per umi
#'
#' @param track
#' @param dsn
#'
#' @return
#'
#' @examples
.gpatterns.reads_per_umi <- function(track, dsn=NULL){
    tcpgs <- gpatterns.get_tidy_cpgs(track) %>% distinct(read_id, num)
    if (!is.null(dsn)){
        tcpgs <- tcpgs %>% sample_n(dsn)
    }
    tcpgs %>% group_by(reads=num) %>% summarise(n=n()) %>% return()
}


#dsns = c(2e5, 5e5, 7e5, 1e6, 1.5e6, 2e6, 3e6)

#' ds stats
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
.gpatterns.ds_stats <- function(track, dsns, bam = NULL, fastq = NULL, intervals=NULL, sort_rand=FALSE){
    stats <- list()
    stats[['cg_num']] <- tibble(dsn=numeric(), cg_num=numeric())
    stats[['meth_calls']] <- tibble(dsn=numeric(), meth_calls=numeric())
    stats[['global_avg_meth']] <- tibble(dsn=numeric(), avg_meth=numeric())
    stats[['reads_per_umi']] <-  tibble(dsn=numeric(), reads=numeric(), n=numeric())
    stats[['umis']] <-  tibble(dsn=numeric(), umis=numeric())
    stats[['on_tar']] <-  tibble(dsn=numeric(), on_tar=numeric())

    sort_rand_str <- if(sort_rand) ' sort -R | ' else ''

    for (dsn in dsns){
        message(qq('doing @{dsn}'))
        tcpgs <- gpatterns.get_tidy_cpgs(track, intervals=intervals)
        all_tcpgs <- gpatterns.get_tidy_cpgs(track)

        if (!is.null(fastq)){
            message('sampling from fastq')
            cmd <- qq('gzip -d -c @{paste(fastq, collapse=" ")} | sed -n 1~4p | @{sort_rand_str} head -n @{dsn}')
            ids <- fread(cmd, sep='\t', header=F)[,1] %>%
                gsub('^@', '', .)
        } else if (!is.null(bam)){
            message('sampling from bam')
            cmd <- qq('samtools view @{bam} | sed -n 1~2p | cut -f1 | @{sort_rand_str} head -n @{dsn}')
            ids <- fread(cmd, sep='\t', header=F)[,1]
        } else {
            message('sampling from tidy cpgs')
            ids <- tcpgs %>% distinct(read_id, num) %>% sample_n(dsn) %>% .$read_id
        }

        if (length(ids) == dsn){
            tcpgs <- tcpgs %>% filter(read_id %in% ids)
            all_tcpgs <- all_tcpgs %>% filter(read_id %in% ids)
            stats[['cg_num']] <- stats[['cg_num']] %>% bind_rows(tibble(dsn=dsn, cg_num = tcpgs %>% distinct(chrom, cg_pos) %>% nrow))
            stats[['meth_calls']] <- stats[['meth_calls']] %>% bind_rows(tibble(dsn=dsn, meth_calls = nrow(tcpgs)))
            stats[['global_avg_meth']] <- stats[['global_avg_meth']] %>% bind_rows(tibble(dsn=dsn, avg_meth = mean(tcpgs$meth)))
            stats[['reads_per_umi']] <-  stats[['reads_per_umi']] %>% bind_rows(tcpgs %>% group_by(reads=num) %>% summarise(n=n()) %>% mutate(dsn = dsn))
            stats[['umis']] <- stats[['umis']] %>% bind_rows(tibble(dsn=dsn, umis = length(unique(tcpgs$read_id))))
            stats[['on_tar']] <- stats[['on_tar']] %>% bind_rows(tibble(dsn=dsn, on_tar = length(unique(tcpgs$read_id)) / length(unique(all_tcpgs$read_id))))
        } else {
            message(qq('skipping @{dsn} (not enough reads)'))
        }
    }

    stats <- stats %>% map(~ .x %>% mutate(track = track))

    return(stats)
}
