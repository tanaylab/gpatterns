# Epipolymorphism Functions ------------------------------------------------

########################################################################
#' Calculate epipolymorphism of samples
#'
#' @param patterns patterns
#'
#' @return
#' @export
#'
#' @examples
gpatterns.calc_epipoly <- function(patterns){
    p <- table(patterns) / length(patterns)
    epipoly <- 1 - sum(p^2)
    return(epipoly)
}

########################################################################
#' Returns TRUE if pattern is a 'border' pattern
#'
#' @param pattern pattern
#'
#' @return
#' @export
#'
#' @examples
gpatterns.border_pattern <- function(pattern){
    !stringr::str_detect(pattern, '01+0|10+1') & !gpatterns.clean_pattern(pattern)
}

########################################################################
#' Returns TRUE if pattern is fully methylated / unmethylated and FALSE if not
#'
#' @param pattern pattern
#'
#' @return
#' @export
#'
#' @examples
gpatterns.clean_pattern <- function(pattern){
    stringr::str_count(pattern, '1') == 0 | stringr::str_count(pattern, '0') == 0
}

########################################################################
#' Plots epipolymorphism
#'
#' @param tracks track names
#' @param intervals genomic scope for which the function is applied
#' @param colnames column names for the samples (useful when track names are too long
#' @param pat_len pattern length
#' @param meth_binsize size of bin for methyltion in x axis
#' @param median_line_binsize number of bins to smooth in median line
#' @param fig_ofn figure file name. if NULL plotting would occur on the current device
#' @param width figure width in inches
#' @param height figure width in inches
#' @param ncolumns number of columns
#' @param separate_figs plot each track in a separate figure
#'
#' @return plot. If seperate.figs=TRUE list of plots
#'
#' @examples
#'
#' @export
gpatterns.epipoly_plot <- function(
    tracks,
    intervals,
    colnames = NULL,
    pat_len = 5,
    median_line_binsize = 6,
    meth_binsize = 0.01,
    fig_ofn = NULL,
    width = NULL,
    height = NULL,
    ncolumns = 3,
    separate_figs = FALSE) {

    tab <- gpatterns.extract(tracks,
                             intervals = intervals,
                             colnames = colnames,
                             elements = c("meth", "epipoly"),
                             tidy = TRUE) %>%
            filter(!is.na(fid))

    m <- seq(0, 1, by = meth_binsize)

    median_line = tab %>%
        group_by(pat_meth = cut(pat_meth, breaks = m, include.lowest = TRUE), samp) %>%
        summarize(epipoly = median(epipoly)) %>%
        ungroup
    levels(median_line$pat_meth) = (m + meth_binsize/2)[-length(m)]
    median_line$pat_meth = as.numeric(as.character(median_line$pat_meth))
    median_line <- median_line %>%
        group_by(samp) %>%
        mutate(epipoly = zoo::rollmean(epipoly,
                                       k = median_line_binsize*2+1,
                                       fill = list(rep(NA, median_line_binsize), NULL, rep(NA, median_line_binsize))))

    num_obs = tab %>%
        group_by(samp) %>%
        summarise(label = sprintf("n = %s",
                                  scales::comma(sum(!is.na(epipoly)))))

    min_poly = tibble(pat_meth = m, epipoly = 2 * m * (1 - m))
    max_poly = tibble(pat_meth = m, epipoly = 1 - ((1 - 2 * m + 2 * m * m)^pat_len))

    if (separate_figs) {
        plotlist <- plyr::dlply(tab, .(samp), function(x) {
            p <- x %>% ggplot(aes(pat_meth, epipoly)) +
                geom_point(size = 0.4, color = "pink", alpha = 0.5) +
                geom_line(data = min_poly, color = "grey", lty = 2) +
                geom_line(data = max_poly, color = "grey") +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) +
                geom_line(data = median_line, color = "red") +
                xlab("average methylation") +
                ylab("epipolymorphism") +
                geom_text(data = num_obs, aes(x = 0.5, y = 0.1, label = label),
                          color = "black",
                          inherit.aes = FALSE,
                          parse = FALSE) +
                theme(strip.background = element_blank(),
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                      panel.margin = unit(2, "lines")) +
                ggtitle(x$samp[1])

            if (!is.null(fig_ofn)) {
                width <- width %||% 4
                height <- height %||% 4
                p <- p + ggsave(sprintf("%s_%s.png", fig_ofn, x$samp[1]), width = width,
                  height = height)
            }
            return(p)
        })
        return(plotlist)

    } else {
        p <- tab %>%
            ggplot(aes(pat_meth, epipoly)) +
            geom_point(size = 0.4, color = "pink", alpha=0.5) +
            geom_line(data = min_poly, color = "grey", lty = 2) +
            geom_line(data = max_poly, color = "grey") +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            geom_line(data = median_line, color = "red") +
            xlab("average methylation") +
            ylab("epipolymorphism") +
            facet_wrap(~samp, ncol = ncolumns) +
            geom_text(data = num_obs,
                      aes(x = 0.5, y = 0.1, label = label),
                      color = "black",
                      inherit.aes = FALSE,
                      parse = FALSE) +
            theme(strip.background = element_blank(),
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                  panel.margin = unit(2,"lines"))

        if (!is.null(fig_ofn)) {
            width <- width %||% 1.6 * ncolumns
            width <- width %||% 1.6 * 2.3 * ceiling(length(unique(tab$samp))/ncolumns)

            p <- p + ggsave(fig_ofn, width = width, height = height)
        }
        return(p)

    }
}

# Bipolar model Functions ------------------------------------------------

#' Run bipolar model
#'
#' @param track           track
#' @param uniform_mix     uniform_mix
#' @param max_sampling_n  max_sampling_n
#' @param min_sampling_n  min_sampling_n
#' @param init_num        init_num
#' @param min_pat_cov     min_pat_cov
#' @param save_tab        save_tab
#' @param parallel        parallel
#' @param thread_num      thread_num
#' @param verbose         verbose
#' @param use_sge         use_sge
#' @param max_jobs        max_jobs
#' @param ...             ...
#'
#' @return
#' @export
#'
#' @examples
gpatterns.run_bipolar_model <- function(track,
                            uniform_mix = 0.05,
                            max_sampling_n = 1000,
                            min_sampling_n = 100,
                            init_num = 3,
                            min_pat_cov = 5,
                            save_tab = FALSE,
                            parallel = getOption('gpatterns.parallel'),
                            thread_num = getOption('gpatterns.parallel.thread_num'),
                            verbose = FALSE,
                            use_sge = FALSE,
                            max_jobs = 300,
                            ...){

    if (use_sge){
        tmp_dirname <- tempfile(pattern = "", tmpdir = paste(get("GROOT"),
            "/tmp", sep = ""))
        if (!dir.create(tmp_dirname, recursive = T, mode = "0777")){
            stop(sprintf("Failed to create a directory %s", tmp_dirname),
                call. = F)
        }

        fn <- tempfile(tmpdir = tmp_dirname)
        group_num <- max_jobs
    } else{
        fn <- tempfile()
        group_num <- thread_num
    }

    patterns <- gpatterns.extract_patterns(track) %>%
        group_by(fid) %>%
        mutate(n = n()) %>%
        filter(n >= min_pat_cov) %>%
        ungroup %>%
        arrange(fid)
    patterns %>%
        write.table(fn, sep='\t', col.names=F, row.names=F, quote=F)
    fid_groups <- patterns %>%
        distinct(fid) %>%
        mutate(grp = ntile(fid, group_num)) %>%
        group_by(grp) %>%
        summarise(min_fid=min(fid), max_fid=max(fid))

    if (use_sge){
        commands <- plyr::alply(fid_groups, 1, function(x)
            qq('.gpatterns.run_bipolar_model("@{fn}",
                                             min_frag_id=@{x$min_fid},
                                             max_frag_id=@{x$max_fid},
                                             uniform_mix=@{uniform_mix},
                                             max_sampling_n=@{max_sampling_n},
                                             init_num=@{init_num},
                                             ret_tbl=TRUE,
                                             verbose=@{verbose},
                                             bipolar_model_bin="@{.gpatterns.bipolar_model_bin}")') %>%
                gsub('\n', ' ', .)
            )

        res <- gcluster.run2(command_list = commands,
                             max.jobs = max_jobs,
                             debug = verbose,
                             packages = 'gpatterns',
                             jobs_title = 'gpatterns.bipolar_model',
                             collapse_results = TRUE,
                             ...)

        unlink(tmp_dirname, recursive = TRUE)

    } else{
        res <- plyr::adply(fid_groups, 1, function(x)
            .gpatterns.run_bipolar_model(fn,
                                         min_frag_id=x$min_fid,
                                         max_frag_id=x$max_fid,
                                         uniform_mix=uniform_mix,
                                         max_sampling_n=max_sampling_n,
                                         init_num=init_num,
                                         ret_tbl=TRUE,
                                         verbose=verbose),
            .parallel=parallel)
        res <- res %>% select(-min_fid, -max_fid, -grp)
    }


    res <- res %>% mutate(qval = p.adjust(pval)) %>%
        tbl_df()

    return(res)
}

########################################################################
.gpatterns.run_bipolar_model <- function(fn,
                                         min_frag_id,
                                         max_frag_id,
                                         uniform_mix=0.05,
                                         max_sampling_n=1000,
                                         min_sampling_n=100,
                                         init_num=3,
                                         ofn = NULL,
                                         ret_tbl = FALSE,
                                         verbose = FALSE,
                                         bipolar_model_bin = .gpatterns.bipolar_model_bin){
    if (is.null(ofn)){
        ofn <- tempfile()
        ret_tbl <- TRUE
    }
    print(bipolar_model_bin)
    command <- paste(bipolar_model_bin,
                     fn,
                     qq('-min_frag_id=@{min_frag_id}'),
                     qq('-max_frag_id=@{max_frag_id}'),
                     qq('-num_of_inits=@{init_num}'),
                     qq('-K=@{max_sampling_n}'),
                     qq('-min_k=@{min_sampling_n}'),
                     qq('-unimodal_uniform_mix=@{uniform_mix}'),
                     qq('-file_format=1 -mode=2 > @{ofn}'))

    system(command, ignore.stderr = !verbose)

    if (ret_tbl){
        res <- read.table(ofn,
                     header = TRUE,
                     colClasses=c(frag_id = 'numeric',
                                  pval = 'numeric',
                                  z_score = 'numeric',
                                  ll_ratio = 'numeric',
                                  m_ratio = 'numeric',
                                  sd_ratio = 'numeric',
                                  mix1 = 'numeric',
                                  center1 = 'character',
                                  gain1 = 'numeric',
                                  loss1 = 'numeric',
                                  mix2 = 'numeric',
                                  center2 = 'character',
                                  gain2 = 'numeric',
                                  loss2 = 'numeric',
                                  uni.mix = 'numeric',
                                  mix.uni = 'numeric',
                                  center.uni = 'character',
                                  gain.uni = 'numeric',
                                  loss.uni = 'numeric',
                                  uni.mix2 = 'numeric')) %>%
            rename(fid = frag_id) %>% tbl_df %>%
            .gpatterns.arrange_biploar_model_output()


        return(res)
    }
}

########################################################################
.gpatterns.arrange_biploar_model_output <- function(tab, calc_qval = TRUE){
    if (calc_qval){
        tab <- tab %>% mutate(qval = p.adjust(pval))
    }

    orig_tab <- tab %>% mutate(mix = mix1 / (1 - uni.mix),
                               mix2 = mix2 / (1 - uni.mix2))
    tab <- tab %>% mutate(major = mix1 > mix2,
                          mix1 = ifelse(major, orig_tab$mix1, orig_tab$mix2),
                          mix2 = ifelse(major, orig_tab$mix2, orig_tab$mix1),
                          center1 = ifelse(major, orig_tab$center1, orig_tab$center2),
                          center2 = ifelse(major, orig_tab$center2, orig_tab$center1),
                          gain1 = ifelse(major, orig_tab$gain1, orig_tab$gain2),
                          gain2 = ifelse(major, orig_tab$gain2, orig_tab$gain1),
                          loss1 = ifelse(major, orig_tab$loss1, orig_tab$loss2),
                          loss2 = ifelse(major, orig_tab$loss2, orig_tab$loss1)) %>%
        select(-major)
    tab <- tab %>%
        mutate(cpgs=nchar(center1),
                          center1_m=cpgs-nchar(gsub('1', '', center1)),
                          center2_m=cpgs-nchar(gsub('1', '', center2))) %>%
        mutate(meth = center1_m >= center2_m,
               mix_meth = ifelse(meth, mix1, mix2),
               mix_unmeth = ifelse(meth, mix2, mix1),
               center_meth = ifelse(meth, center1, center2),
               center_unmeth = ifelse(meth, center2, center1),
               gain_meth = ifelse(meth, gain1, gain2),
               gain_unmeth = ifelse(meth, gain2, gain1),
               loss_meth = ifelse(meth, loss1, loss2),
               loss_unmeth = ifelse(meth, loss2, loss1),
               center_meth.ones = ifelse(meth, center1_m, center2_m),
               center_unmeth.ones = ifelse(meth, center2_m, center1_m),
               center_uni.ones = cpgs-nchar(gsub('1', '', center.uni))) %>%
        rename(mix_uni=uni.mix,
               center_uni=center.uni,
               gain_uni=gain.uni,
               loss_uni=loss.uni,
               uni_mix2=uni.mix2) %>%
        select(-meth, -center1_m, -center2_m, -cpgs, -mix.uni)

    return(tab)


}

# other Functions ------------------------------------------------

########################################################################
gpatterns.theta_distance <- function(track1, track2, meth_range=NULL, min_cov=NULL, similarity=FALSE){

    pats1 <- gpatterns.extract_patterns(track1, tabular=T) %>%
        group_by(fid, pattern) %>%
        summarise(n1=n()) %>%
        group_by(fid) %>%
        mutate(p1 = n1 / sum(n1)) %>%
        ungroup()

    pats2 <- gpatterns.extract_patterns(track2, tabular=T) %>%
        group_by(fid, pattern) %>%
        summarise(n2=n()) %>%
        group_by(fid) %>%
        mutate(p2 = n2 / sum(n2)) %>%
        ungroup()

    if (!is.null(min_cov)){
        pats1 <- pats1 %>%
            group_by(fid) %>%
            mutate(n = n()) %>%
            filter(n >= min_cov) %>%
            ungroup()
        pats2 <- pats2 %>%
            group_by(fid) %>%
            mutate(n = n()) %>%
            filter(n >= min_cov) %>%
            ungroup()
    }

    fids <- pats1$fid[pats1$fid %in% pats2$fid]

    pats <- full_join(pats1, pats2, by=c('fid', 'pattern') ) %>%
        filter(fid %in% fids) %>%
        mutate_each(funs(na2zero = ifelse(is.na(.), 0, .)))

    if (!is.null(meth_range)){
        fids1 <- .gpatterns.load_fids_tab(track1)
        fids2 <- .gpatterns.load_fids_tab(track2)
        pats <- pats %>%
            left_join(fids1 %>% select(fid, meth1=meth), by='fid') %>%
            left_join(fids2 %>% select(fid, meth2=meth), by='fid') %>%
            filter(between(meth1, meth_range[1], meth_range[2]), between(meth2, meth_range[1], meth_range[2]))
    }

    if (similarity){
        res <- pats %>%
            group_by(fid) %>%
            summarise(inter = sum(p1*p2),
                      intra = sqrt( (sum(p1^2)) * (sum(p2^2)) ),
                      theta=inter / intra) %>%
            filter(!is.na(theta))
    } else {
        res <- pats %>%
            group_by(fid) %>%
            summarise(inter = 1 - sum(p1*p2),
                      intra = sqrt( (1 - sum(p1^2)) * (1 - sum(p2^2)) ),
                      theta=inter / intra) %>%
            filter(!is.na(theta))
    }

    return(res)
}

########################################################################
gpatterns.theta_distance_sampling <- function(track1, track2, sampling_n, meth_range=NULL, rm_n0=FALSE, rm_n1=FALSE, rm_borders=FALSE, similarity=FALSE, replace=FALSE, min_cov=NULL, add_frag_stats=FALSE, add_intra=FALSE){
    if (!replace && is.null(min_cov)){
        min_cov <- sampling_n
    }

    pats1 <- gpatterns.extract_patterns(track1, tabular=T)
    pats1 <- .gpatterns_filter_pats(pats1, min_cov=min_cov, meth_range=meth_range, rm_n0=rm_n0, rm_n1=rm_n1, rm_borders=rm_borders, track=track1)
    if (add_frag_stats){
       frag_stats1 <- .gpatterns.frag_stats(pats1) %>%
            rename(ncpg_s1=ncpg, n_s1=n, n0_s1=n0, n1_s1=n1, nx_s1=nx, nc_s1=nc, meth_s1=meth, epipoly_s1=epipoly)
    }

    if (track1 == track2){
        res <- .calc_epipoly_sampling_intra(pats1, sampling_n, replace=replace, similarity=similarity)
        if (add_frag_stats){
            res <- res %>% left_join(frag_stats1, by='fid')
        }
        return(res)
    }

    pats2 <- gpatterns.extract_patterns(track2, tabular=T)
    pats2 <- .gpatterns_filter_pats(pats2, min_cov=min_cov, meth_range=meth_range, rm_n0=rm_n0, rm_n1=rm_n1, rm_borders=rm_borders, track=track2)

    fids <- pats1$fid[pats1$fid %in% pats2$fid]
    pats1 <- pats1 %>%
        filter(fid %in% fids)
    pats2 <- pats2 %>%
        filter(fid %in% fids)

    res <- .calc_epipoly_sampling(pats1, pats2, sampling_n, replace=replace, similarity=similarity)

    if (add_intra){
        intra1 <- .calc_epipoly_sampling_intra(pats1, sampling_n, replace=replace, similarity=similarity)
        intra2 <- .calc_epipoly_sampling_intra(pats2, sampling_n, replace=replace, similarity=similarity)
        res <- bind_cols(
            res %>% select(fid, inter=epipoly),
            intra1 %>% select(intra1=epipoly),
            intra2 %>% select(intra2=epipoly))
    }


    if (add_frag_stats){
        frag_stats2 <- .gpatterns.frag_stats(pats2) %>%
            rename(ncpg_s2=ncpg, n_s2=n, n0_s2=n0, n1_s2=n1, nx_s2=nx, nc_s2=nc, meth_s2=meth, epipoly_s2=epipoly)
        res <- res %>% left_join(frag_stats1, by='fid') %>% left_join(frag_stats2, by='fid')
    }

    return(res)
}


########################################################################
.gpatterns_filter_pats <- function(pats, min_cov=NULL, meth_range=NULL, rm_n0=FALSE, rm_n1=FALSE, rm_borders=FALSE, track=NULL){
    pats <- pats %>%
        mutate(
            cpgs = nchar(pattern),
            ones=cpgs-nchar(gsub('1', '', pattern)))

    if (rm_n0){
        pats <- pats %>% filter(ones != 0)
    }
    if (rm_n1){
        pats <- pats %>% filter(ones != cpgs)
    }
    if (rm_borders){
        pats <- pats %>% filter(!.gpatterns.border_pattern(pattern))
    }
    if (!is.null(meth_range)){
        if (is.null(track)){
            stop('no track was supplied for methylation range filtering')
        }
        fids <- .gpatterns.load_fids_tab(track)
        pats <- pats %>% left_join(fids %>% select(fid, meth), by='fid') %>% filter(between(meth, meth_range[1], meth_range[2])) %>% select(-meth)
    }

    if (!is.null(min_cov)){
        pats <- pats %>%
            group_by(fid) %>%
            filter(n() >= min_cov) %>%
            ungroup()
    }

    return(pats)
}

########################################################################
.calc_epipoly_sampling <- function(pats1, pats2, sampling_n, replace=FALSE, similarity=FALSE){
    pats1 <- pats1 %>%
        arrange(fid) %>%
        group_by(fid) %>%
        sample_n(sampling_n, replace=replace) %>%
        select(fid, pattern1=pattern)
    pats2 <- pats2 %>%
        arrange(fid) %>%
        group_by(fid) %>%
        sample_n(sampling_n, replace=replace) %>%
        select(fid2=fid, pattern2=pattern)
    res <- bind_cols(pats1, pats2) %>%
        filter(fid == fid2) %>%
        group_by(fid) %>%
        summarise(epipoly = sum(pattern1 != pattern2) / n())

    if (similarity){
        res <- res %>%
            mutate(epipoly = 1 - epipoly)
    }
    return(res)
}

########################################################################
.calc_epipoly_sampling_intra <- function(pats, sampling_n, replace=FALSE, similarity=FALSE){
    pats_min2 <- pats %>% group_by(fid) %>% filter(n() >= 2) %>% ungroup
    res <- plyr::adply(1:sampling_n, 1, function(x)
        pats_min2 %>%
            group_by(fid) %>%
            sample_n(2, replace=replace) %>%
            summarise(dif = pattern[1] != pattern[2])) %>%
            group_by(fid) %>%
            summarise(epipoly = sum(dif, na.rm=T) / n()) %>%
            ungroup
    res <- count(pats, fid) %>%
        left_join(res, by='fid') %>%
        select(-n) %>%
        arrange_old(fid)

    if (similarity){
        res <- res %>%
            mutate(epipoly = 1 - epipoly)
    }
    return(res)
}

########################################################################
gpatterns.theta_matrix_sampling <- function(tracks, samples=NULL, parallel=FALSE, symetric=TRUE, ...){
    if (parallel){
        .progress='none'
    } else {
        .progress='text'
    }
    if (is.null(samples)){
        samples <- tracks
    }
    if (length(samples) != length(tracks)){
        stop('number of samples needs to be the same as the number of tracks')
    }

    combs <- t(combn(tracks, 2)) %>% as.data.frame %>%
            rename(track1=V1, track2=V2) %>%
            bind_rows(data.frame(track1=tracks, track2=tracks))

    theta <- plyr::adply(combs, 1, function(x) {
            a <- gpatterns.theta_distance_sampling(
                track1=as.character(x$track1),
                track2=as.character(x$track2),
                add_intra=FALSE,
                ...)
            data.frame(theta=mean(a$epipoly), loci=nrow(a))
        }, .progress=.progress, .parallel=parallel )

    if (symetric){
        #make theta[i,j] equals theta[j,i]
        theta <- theta %>%
                bind_cols(t(combn(samples, 2)) %>% as.data.frame %>% bind_rows(data.frame(V1=samples, V2=samples))) %>%
                rename(samp1=V1, samp2=V2) %>%
                select(samp1, samp2, theta, loci)
            theta1 <- data.frame(samp1=theta$samp2, samp2=theta$samp1, theta=theta$theta, loci=theta$loci)
            theta <- bind_rows(theta, theta1) %>%
                distinct(samp1, samp2, .keep_all=T)
    } else {
        theta <- theta %>%
            bind_cols(expand.grid(samples, samples)) %>%
            rename(samp1=Var1, samp2=Var2) %>%
            select(samp1, samp2, theta, loci)
    }

    return(theta)
}
