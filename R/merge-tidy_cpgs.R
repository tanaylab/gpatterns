#' merge tidy cpgs directories
#' @export
gpatterns.merge_tidy_cpgs <- function(dirs, out_dir=tempdir(), nbins=nrow(gintervals.all()), paired_end = TRUE, stats_dir = paste0(out_dir, "/stats"), filter_dups_bin=system.file("import", "filter_dups_cpgs.py", package="gpatterns")){
    genomic_bins <- gbin_intervals(intervals = gintervals.all(), nbins)

    fn_df <- map_dfr(dirs, ~ tibble(fn = list.files(.x, pattern="tcpgs.gz", full.names=TRUE)) %>% mutate(name = gsub(".tcpgs.gz$", "", basename(fn))) %>% separate(name, c("chrom", "start", "end"), sep="_") %>% mutate(start = as.numeric(start), end = as.numeric(end))) %>% select(chrom, start, end, everything())

    bins_df <-  genomic_bins %>% gintervals.neighbors1(fn_df, maxneighbors=nrow(fn_df), mindist=0, maxdist=0) %>% filter(dist == 0) %>% distinct(chrom, start, end, chrom1, start1, end1, fn)

    dir.create(out_dir, showWarnings = FALSE)
    dir.create(stats_dir, showWarnings = FALSE)

    bins_df_cmd <- bins_df %>% 
        group_by(chrom, start, end) %>%
        summarise(
                end_cond = glue("$4 != \"-\" && $4 <= $3 && $4 <= {end[1]} && $4 >= {start[1]}"), 
                start_cond = glue("$3 != \"-\" && $3 < $4 && $3 <= {end[1]} && $3 >= {start[1]}"), 
                awk_cmd = glue("awk -F',' 'NR==1 || ({start_cond}) || ({end_cond})'"), 
                sort_cmd = glue("awk 'NR==1; NR > 1 {{print $0 | \"sort --field-separator=, -k2,7 -k1 -k9\"}}'"),  
                stats_fn = glue("{stats_dir}/{chrom[1]}_{start[1]}_{end[1]}.stats"),
                out_fn = glue("{out_dir}/{chrom[1]}_{start[1]}_{end[1]}.tcpgs.gz"),
                fns = paste(fn, collapse=" "),
                filter_dups_cmd = glue("{filter_dups_bin} -i - -s {stats_fn} --sorted"),
                filter_dups_cmd = ifelse(paired_end, filter_dups_cmd, paste(filter_dups_cmd, "--only_R1")),
                cmd = glue("gzip -d -c {fns} | {awk_cmd} | {sort_cmd} | {filter_dups_cmd} | gzip -c > {out_fn}"))

    plyr::l_ply(bins_df_cmd$cmd, function(x) system(x), .parallel=TRUE)
}