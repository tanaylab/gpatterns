#' Impute missing methylation calls
#'
#' @param avgs 'tidy' output of \code{gpatterns.get_avg_meth}.
#' @param k_reg k_reg
#' @param min_locus_cov minimal total coverage per locus. If not NULL, filters out loci with less
#' than \code{min_locus_cov} methylation calls.
#' @param tidy return 'tidy' output.
#'
#' @return if \code{tidy}, \code{avgs} with additional \code{meth_reg}, \code{cov_reg} and \code{avg_reg} fields with the imputed methylation calls, total number of calls and avergae methylation (respectively) for each sample.
#' If \code{tidy == FALSE}: intervals set with additional fields with regularized average methylation for each sample.
#' 
#' @export
#'
gpatterns.impute <- function(avgs, k_reg=3, min_locus_cov=NULL, tidy=TRUE){
    meth_tab <- reshape2::dcast(avgs, chrom + start + end  ~ samp, value.var='meth', fill=0) %>% as_tibble()
    cov_tab <- reshape2::dcast(avgs, chrom + start + end  ~ samp, value.var='cov', fill=0) %>% as_tibble()

    if (!is.null(min_locus_cov)){
        message(glue('taking only that are covered >= {min_locus_cov}'))
        cov_f <- apply(cov_tab %>% select(-(chrom:end)), 1, sum, na.rm=TRUE) >= min_locus_cov
        cov_tab <- cov_tab[cov_f, ]
        meth_tab <- meth_tab[cov_f, ]
        message(glue('{comify(sum(cov_f))} loci remaining ({scales::percent(sum(cov_f) / length(cov_f))} out of {comify(length(cov_f))})'))
    }

    avg_m <- rowSums(meth_tab[, -1:-3], na.rm=T) / rowSums(cov_tab[, -1:-3], na.rm=TRUE)

    meth_reg_tab <- bind_cols(meth_tab[, 1:3], meth_tab[, -1:-3] + k_reg * avg_m )
    cov_reg_tab <- bind_cols(cov_tab[, 1:3], cov_tab[, -1:-3] + k_reg )
    met_reg <- bind_cols(meth_reg_tab[, 1:3], meth_reg_tab[, -1:-3] / cov_reg_tab[, -1:-3])

    if (tidy){
        met_reg <- meth_reg_tab %>%
            gather('samp', 'meth_reg', -(chrom:end)) %>%
            inner_join(cov_reg_tab %>%
                gather('samp', 'cov_reg', -(chrom:end)),
                    by=c('chrom', 'start', 'end', 'samp')) %>%
            inner_join(met_reg %>%
                gather('samp', 'avg_reg', -(chrom:end)),
                    by=c('chrom', 'start', 'end', 'samp')) %>%
            inner_join(avgs, by=c('chrom', 'start', 'end', 'samp'))
    }
    return(met_reg)
}

#' Filter loci with low variance across samples
#'
#' Clusters the loci with kmeans with a relatively high k, calculates the \code{sd} of each cluster centers, 
#' and then removes loci that are within clusters with low variance. 
#'
#' @param avgs 'tidy' output of \code{gpatterns.get_avg_meth}. In order to filter using the regularized (imputed) values (output of \code{gpatterns.impute}), set \code{avg_col} to 'avg_reg'
#' @param k number of clusters to divide to. if NULL - k would be chosen as number of loci divided by 150, 
#' in order to have ~100 loci per cluster. 
#' @param center_sd_thresh minimal sd of the cluster centers. if NULL - \code{center_sd_thresh} would be set by taking
#' the sd of the cluster with the largest difference to the next cluster sd that still
#' leaves \code{min_loci_frac} loci. 
#' @param min_loci_frac minimal fraction of loci to return in case \code{center_sd_thesh == NULL}.
#' @param seed seed to use in \code{TGL_kmeans_tidy}
#' @param tidy return tidy output
#' @param avg_col column of average methylation in \code{avgs}. 
#' @param plot_cluster_sd plot the standard deviation of clusters.
#' @param ret_clust_sd return a list with the filtered loci (under 'avgs'), and the cluster sd plot (under 'clust_sd_p')
#'
#' @return see \code{ret_clust_sd}. if \code{tidy}: \code{avgs} data frame without the filtered loci. if \code{tidy == FALSE}: intervals set with additional fields with regularized average methylation for each sample.
#' 
#' @export
#'
gpatterns.filter_loci <- function(avgs, k=NULL, center_sd_thresh=NULL, min_loci_frac=0.3, seed=NULL, tidy=TRUE, avg_col='avg', plot_cluster_sd = FALSE, ret_clust_sd=FALSE){

    n_inital <- distinct(avgs, chrom, start, end) %>% nrow()
    met_reg <- avgs %>% select(one_of('chrom', 'start', 'end', 'samp', avg_col)) %>% spread_('samp', avg_col)

    if (is.null(k)){
        k <- floor(nrow(met_reg) / 150)
        message(glue('setting k to {k}'))
    }    

    km <- TGL_kmeans_tidy(met_reg %>% select(-(chrom:end)), k=k, id_column=FALSE, seed=seed)
    clust_sd <- sqrt(apply(km$centers[, -1] - rowMeans(km$centers[, -1]), 1, var))

    d <- km$size %>% mutate(csd = clust_sd) %>% arrange(csd) %>% mutate(cn = rev(cumsum(rev(n))))

    if (is.null(center_sd_thresh)){
        center_sd_thresh <- d %>% mutate(diff_sd = lead(csd) - csd) %>% filter(cn >= max(cn)*min_loci_frac) %>% filter(diff_sd == max(diff_sd, na.rm=T)) %>% pull(csd)
        message(glue('setting center sd thershold to: {center_sd_thresh}'))
    }

    d <- d %>% mutate(below_f = csd < center_sd_thresh)

    ggp <- d %>% ggplot(aes(x=cn, y=csd)) + geom_point() + xlab('Loci remaining') + ylab('Cluster sd') + geom_hline(yintercept=center_sd_thresh, linetype='dashed', color='red') + geom_vline(xintercept=sum(d$n[!d$below_f]), linetype='dashed', color='red') + scale_x_continuous(label=comify)

    if (plot_cluster_sd){        
        print(ggp)
    }

    var_clusters <- km$centers$clust[clust_sd >= center_sd_thresh]

    message(glue('taking only loci within clusters that have sd >= {center_sd_thresh}'))
    met_reg_f <- bind_cols(km$cluster %>% select(clust), met_reg) %>% filter(clust %in% var_clusters) %>% select(-clust)

    message(glue('{comify(nrow(met_reg_f))} loci remaining, ({scales::percent(nrow(met_reg_f) / n_inital)} out of {comify(n_inital)})'))

    if (tidy){
        met_reg_f <- met_reg_f %>%
            gather('samp', avg_col, -(chrom:end)) %>%
            inner_join(avgs, by=c('chrom', 'start', 'end', 'samp'))
    }

    if (ret_clust_sd){
        return(list(avgs=met_reg_f, clust_sd_p=ggp))
    }
    return(met_reg_f)
}

#' Get background trend expected methylation signal
#' 
#' Calculate the expected average methylation given the background (CpG content driven) trend.
#'
#' @param intervals genomic scope for which the function is applied
#' @param tracks gpatterns tracks
#' @param trend pre-computed global methylation trend (output of \code{gpatterns.global_meth_trend}).
#' @param cg_cont_track CpG content track
#' @param cg_cont_breaks breaks to determine the bins of \code{cg_cont_track}
#' @param ... additional paramters to \code{gpatterns.global_meth_trend}
#'
#' @return intervals set with the additional fields 'exp_meth' 'exp_cov' and 'exp_avg' with the 
#' expected methylation calls, coverage and average methylation (respectively)
#' 
#' @export
#'
#' @inheritDotParams gpatterns.global_meth_trend
gpatterns.get_bg_meth <- function(intervals, tracks=NULL, trend=NULL, cg_cont_track = .gpatterns.cg_cont_500_track, cg_cont_breaks = seq(0, 0.15, by= 0.002), ...){

    if (is.null(tracks)){
        if (!has_name(intervals, 'track')){
            stop('Please provide either tracks or intervals data frame with "track" field')
        } else {
            tracks <- unique(intervals$track)
        }
    }

    if (is.null(trend)){
        trend <- gpatterns.global_meth_trend(tracks=tracks, names=NULL, strat_track = cg_cont_track, strat_breaks=cg_cont_breaks, ...) %>% .$trend
    }

    intervs <- intervals %>% distinct(chrom, start, end)
  
    covs <- gextract.left_join1(c(paste0(tracks, '.cov'), cg_cont_track), intervals=intervs, iterator=.gpatterns.genome_cpgs_intervals, colnames=c(tracks, 'cg_cont')) %>% as_tibble()   

    exp_avgs <- covs %>%
        mutate(breaks = cut(cg_cont, breaks=cg_cont_breaks, include.lowest=TRUE),
            breaks = gsub(',', ', ', as.character(breaks))) %>%
        select(-(chrom:end)) %>%
        gather('track', 'cov', -(cg_cont:breaks)) %>%
        left_join(trend %>% select(track, breaks, exp_meth=meth), by=c('breaks', 'track')) %>%
        mutate(exp_meth = exp_meth * cov) %>%
        group_by(chrom1, start1, end1, track) %>%
        summarise(exp_meth = sum(exp_meth, na.rm=TRUE),
                  exp_cov = sum(cov, na.rm=TRUE),
                  cg_cont=breaks[1]) %>% 
        rename(chrom=chrom1, start=start1, end=end1) %>%
        mutate(exp_avg = exp_meth / exp_cov)        

    return(exp_avgs)
}

gpatterns.cluster <- function(avgs, min_locus_cov = 5, k_reg=3, center_sd_thresh=0.04, tracks=NULL, names=NULL, cluster_samples=TRUE, seed=NULL, hclust_loci=FALSE, bootstrap=FALSE, N_boot=100, boot_ratio=0.75, norm_by_locus_mean = TRUE, parallel=getOption("gpatterns.parallel"), plot_coclust_heatmap=TRUE,  ...){
    tracks <- tracks %||% unique(avgs$track)
    names <- names %||% unique(avgs$samp)
    
    message('Normalizing by background trend...')
    exp_avgs <- gpatterns.get_bg_meth(avgs, tracks=tracks, ...) %>% left_join(avgs, by=c('chrom', 'start', 'end', 'track'))
    
    message('Imputing missing data...')
    # impute expected (without bg trend)
    exp_avgs <- exp_avgs %>%
        select(chrom, start, end, samp, meth=exp_meth, cov) %>%
        gpatterns.impute(tidy=TRUE, k_reg=k_reg, min_locus_cov=min_locus_cov) %>%
        rename(exp_avg_reg = avg_reg, exp_meth=meth, exp_meth_reg=meth_reg, exp_cov_reg=cov_reg) %>%
        inner_join(exp_avgs %>% select(-exp_avg, -exp_meth),
                   by=c('chrom', 'start', 'end', 'samp', 'cov'))
    
    # impute observed (same exp_avgs object has 'meth', 'unmeth', 'cov', 'avg' fields of the observed)    
    met_reg <- gpatterns.impute(exp_avgs, k_reg=k_reg, tidy=TRUE, min_locus_cov=min_locus_cov)

    # calculate obs - exp
    met_reg <- met_reg %>%
        mutate(avg_diff = avg_reg - exp_avg_reg)
    
    # normalize rows
    diff_col <- 'avg_diff'
    if (norm_by_locus_mean){
        message('Dividing by locus mean...')
        met_reg <- met_reg %>% group_by(chrom, start, end) %>%
            mutate(norm_diff = avg_diff - mean(avg_diff, na.rm=TRUE)) %>%
            ungroup()    
        diff_col <- 'norm_diff'
    }     
    

    message('Filtering...')
    met_reg_f <- gpatterns.filter_loci(met_reg, seed=seed, avg_col=diff_col, tidy=TRUE, center_sd_thresh=center_sd_thresh, plot_cluster_sd=FALSE, ret_clust_sd=TRUE)
    clust_sd_p <- met_reg_f$clust_sd_p
    met_reg_f <- met_reg_f$avgs

    message('Clustering...')
    met_diff <- met_reg_f %>% reshape2::dcast(chrom + start + end ~ samp, value.var = diff_col, fill=0) %>% as_tibble()

    met_mat <- met_diff[, -1:-3]    
    bt <- bootclust(met_mat, N_boot=N_boot, boot_ratio=boot_ratio, id_column=FALSE, seed=seed)

    if (plot_coclust_heatmap){        
        tglkmeans::plot_coclust_mat(bt)
    }
    
    if (cluster_samples){
        cm_samp <- tgstat::tgs_cor(as.matrix(met_mat), pairwise.complete.obs = TRUE)
        # cm_samp <- cor(met_mat, use='pairwise.complete.obs')
        colnames(cm_samp) <- colnames(met_mat)
        rownames(cm_samp) <- colnames(met_mat)
        # hc_samp <- hclust(dist(cm_samp), 'ward.D2')
        hc_samp <- hclust(tgstat::tgs_dist(cm_samp), 'ward.D2')
    } else {
        cm_samp <- NULL
        hc_samp <- NULL
    }

    if (hclust_loci){
        message('calculating loci correlation')
        cm <- tgstat::tgs_cor(t(as.matrix(met_mat)), pairwise.complete.obs = TRUE)
        # cm <- cor(t(met_mat), use='pairwise.complete.obs')
        loc_names <- paste0(met_diff$chrom, '_', met_diff$start, '_', met_diff$end)
        colnames(cm) <- loc_names
        rownames(cm) <- loc_names
        message('clustering loci (hclust)')        
        # hc <- hclust(dist(cm), "ward.D2")
        hc <- hclust(tgstat::tgs_dist(cm), "ward.D2")
        
    } else {
        cm <- NULL
        hc <- NULL
    }    

    return(compact(list(avgs=met_reg_f, met_diff=met_diff, cm=cm, hc=hc, cm_samp=cm_samp, hc_samp=hc_samp, clust_sd_p=clust_sd_p, bt=bt)))
}


gpatterns.cutree <- function(avgs, k, min_coclust=0.5, reg_by_cluster = TRUE, k_reg_clust=5){
    met_reg_f <- avgs$avgs
    met_diff <- avgs$met_diff
    bt <- avgs$bt

    bt <- bt %>% cutree_bootclust(k=k, min_coclust=min_coclust, tidy=TRUE)

    met_reg_f <- met_reg_f %>% left_join(met_diff %>% select(chrom, start, end) %>% mutate(clust = bt$clust$clust), by=c('chrom', 'start','end'))

    if (reg_by_cluster){
        message('regularizing for plotting (by cluster)')

        met_clust_reg <- met_reg_f %>%
            group_by(samp, clust) %>%
            mutate(clust_avg = sum(meth_reg, na.rm=TRUE) / sum(cov_reg, na.rm=TRUE),
                   exp_clust_avg = sum(exp_meth_reg, na.rm=TRUE) / sum(exp_cov_reg, na.rm=TRUE)) %>%
            ungroup() %>%
            mutate(avg_reg_clust = (meth_reg + k_reg_clust * clust_avg) / (cov_reg + k_reg_clust),
                exp_avg_reg_clust = (exp_meth_reg + k_reg_clust * exp_clust_avg) / (exp_cov_reg + k_reg_clust)) %>%
            mutate(clust_reg_diff = avg_reg_clust - exp_avg_reg_clust) %>%
            group_by(chrom, start, end) %>%
            mutate(norm_diff_reg = clust_reg_diff - mean(clust_reg_diff, na.rm=TRUE)) %>%
            ungroup()
    } else {
        met_clust_reg <- NULL
    }

    met_diff <- met_diff %>% mutate(clust = bt$clust$clust) %>% select(chrom, start, end, clust, everything())

    avgs$avgs <- met_reg_f
    avgs$met_diff <- met_diff
    avgs$bt <- bt
    avgs$avgs_reg <- met_clust_reg

    return(avgs)
}

gpatterns.plot_clustering <- function(clust, fig_fn=NULL, width=7, height=14, device='pdf', samp_ord=NULL, top_annotation=NULL, bottom_annotation=NULL, annotation_colors=NULL, show_gene_names=FALSE, points=NULL, col=circlize::colorRamp2(c(-1,-0.25, 0, 0.25, 1), c("black", "#00688B", "white", "#FF413D", "yellow")), cluster_rows = TRUE, cluster_columns=TRUE, show_column_names=FALSE, show_row_dend = FALSE, use_raster=FALSE, row_names_gp = gpar(fontsize = 4), column_order=NULL, draw_hm=TRUE, split_by_cluster=TRUE, ...){

    if (tibble::has_name(clust, 'avgs_reg')){
        met_diff <- clust[['avgs_reg']] %>% reshape2::dcast(chrom + start + end + clust ~ samp, value.var = 'norm_diff_reg', fill=0) %>% as_tibble()
    } else {
        met_diff <- clust[['met_diff']]
    }

    met_mat <- met_diff %>% select(-(chrom:clust))

    if (!is.null(top_annotation)){
        top_col <- annotation_colors[colnames(top_annotation)[-1]]
        if (!is.null(points)){
            points <- tibble(samp = colnames(met_mat)) %>% left_join(points) %>% select(-samp) %>% pull(1)
            ha_top <- HeatmapAnnotation(df = as.data.frame(tibble(samp = colnames(met_mat)) %>% left_join(top_annotation, by='samp') %>% select(-samp)), col=top_col, points=anno_points(points, axis=TRUE), annotation_height = unit(c(10, 40), "points"))
        } else {
            ha_top <- HeatmapAnnotation(df = as.data.frame(tibble(samp = colnames(met_mat)) %>% left_join(top_annotation, by='samp') %>% select(-samp)), col=top_col)
        }
    } else {
        ha_top <- NULL
    }

    if (!is.null(bottom_annotation)){
        bottom_col <- annotation_colors[colnames(bottom_annotation)[-1]]
        ha_bottom <- HeatmapAnnotation(df = as.data.frame(tibble(samp = colnames(met_mat)) %>% left_join(bottom_annotation, by='samp') %>% select(-samp)), col=bottom_col)
    } else {
        ha_bottom <- NULL
    }

    if (show_gene_names){
        gene_names <- met_diff %>%
            gintervals.neighbors1('intervs.global.tss')

        subset <- which(abs(gene_names$dist) < 2000)

        # gene_names <- gene_names %>%
        #   mutate(name = as.character(name)) %>%
        #   mutate(name = ifelse(abs(dist) < 2000, name, '')) %>%
        #   pull(name)
        ha_gene_names <- rowAnnotation(text = row_anno_link(at=subset, labels=as.character(gene_names$name[subset]), labels_gp=gpar(fontsize = 4)))

        # ha_gene_names <- rowAnnotation(text = row_anno_text(gene_names, gp=gpar(fontsize = 4)))
    }

    if (!is.null(column_order)){
        column_order <- match(column_order, colnames(met_mat))
    }

    if (cluster_columns){
        cluster_columns <- clust$hc_samp
    }

    if (!is.null(fig_fn)){
        do.call(device, list(fig_fn, width=width, height=height))
    }

    if (split_by_cluster){
        split <- met_diff$clust
    } else {
        split <- NULL
    }

    met_hm <- Heatmap(met_mat, name='diff', heatmap_legend_param = list(color_bar = "continuous"), col=col, cluster_rows = cluster_rows, cluster_columns=cluster_columns, show_column_names=show_column_names, top_annotation=ha_top, bottom_annotation=ha_bottom, split=split, column_title=glue('{comify(nrow(met_mat))} loci'), show_row_dend = show_row_dend, use_raster=use_raster, row_names_gp = row_names_gp, column_order=column_order, ...)

    if (show_gene_names){
        res <- met_hm + ha_gene_names
    } else {
        res <- met_hm
    }

    if (draw_hm){
        draw(res)
    }

    if (!is.null(fig_fn)){
        dev.off()
    }
    invisible(res)
}
