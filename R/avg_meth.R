# Extraction Functions ------------------------------------------------


#' Extract average methylation
#'
#' @param tracks methylation tracks
#' @param intervals genomic scope for which the function is applied
#' @param iterator see iterator in \code{\link[gextract]{misha}}. if NULL iterator
#' would be set to CpGs
#' @param min_cov minimal coverage for iterator interval
#' @param mask_by_cov change loci with coverage < min_cov to NA. Not relevant in
#' pre_screen or tidy mode
#' @param use_cpgs use CpGs as iterator
#' @param min_samples minimal number of samples with cov >= min_cov. if min_cov
#' is NULL it would be set to 1.
#' @param min_cpgs minimal number of CpGs per iterator interval. note that the
#' intervalID column may be incorrect.
#' @param min_var minimal variance (across samples) per iterator interval
#' @param var_quantile minimal quantile of variance per iterator interval
#' @param names alternative names to tracks. similar to colnames in
#' \code{\link[gextract]{misha}} if tidy == FALSE. Note that names should be
#' shorter than the maximal length of R data frame column name
#' @param tidy if TRUE returns a tidy data frame with the following fields:
#' chrom, start, end, intervalID, samp, meth, unmeth, avg, cov.
#' if FALSE returns a data frame with average methylation,
#' similar to \code{\link[gextract]{misha}}'. Note that for a large number of
#' intervals tidy == FALSE may be the only memory feasable option.
#' @param pre screen for min_samples and min_cov (for large number of
#' tracks / large number of intervals). Note that the intervalID column may be incorrect
#' and if use_cpgs is TRUE, the intervals set would become the cpgs.
#' @use_disk for really big datasets - save intermediates on disk.
#' @param file save output to file (only in non tidy mode, would not filer by variance)
#' @param intervals.set.out save output big intervals set (only in tidy mode,
#' would not filter by variance)
#'
#'
#'
#' @return
#' @export
#'
#' @examples
gpatterns.get_avg_meth <- function(
    tracks,
    intervals,
    iterator = NULL,
    min_cov = NULL,
    mask_by_cov = FALSE,
    use_cpgs = FALSE,
    min_samples = NULL,
    min_cpgs = NULL,
    min_var = NULL,
    var_quantile = NULL,
    names = NULL,
    tidy = TRUE,
    pre_screen = FALSE,
    use_disk = FALSE,
    file = NULL,
    intervals.set.out = NULL) {

    .check_tracks_exist(tracks, c('meth', 'unmeth'))
    min_samples <- min_samples %||% 1

    if (is.null(iterator) && length(tracks) == 1){
        iterator <- .gpatterns.meth_track_name(tracks)
    }

    if ((is.null(iterator) && length(tracks) > 1) || use_cpgs) {
        message("Using cpgs as iterator")
        iterator <- .gpatterns.genome_cpgs_intervals
    } else if (!is.null(min_cpgs)) {
        message(qq('Taking only intervals with at least @{min_cpgs} CpGs'))
        intervals <- gpatterns.filter_cpgs(intervals, min_cpgs)
    }

    names <- names %||% tracks

    if (pre_screen || (!tidy & (!is.null(min_samples) || !is.null(min_cov)))){
        message(qq("Taking only intervals with coverage >= @{min_cov} in at least @{min_samples} samples"))
        if (use_disk){
            f_intervs <- .random_track_name()
        } else {
            f_intervs <- NULL
        }
        intervals <- gpatterns.screen_by_coverage(tracks = tracks,
                                                  intervals = intervals,
                                                  iterator = iterator,
                                                  min_cov = min_cov,
                                                  min_samples = min_samples,
                                                  intervals.set.out = f_intervs)
        if (use_disk){
            intervals <- f_intervs
        }

    }

    if (!tidy){
        return(
            .gpatterns.get_avg_meth_not_tidy(
                tracks,
                intervals,
                iterator = iterator,
                min_var = min_var,
                var_quantile = var_quantile,
                names = names,
                file = file,
                intervals.set.out = intervals.set.out,
                rm_intervals = !is.null(f_intervs))
        )
    }

    message('extracting...')
    avgs <- gvextract(
        c(qq('@{tracks}.meth', collapse=F), qq('@{tracks}.unmeth', collapse=F)),
        intervals = intervals,
        iterator = iterator,
        colnames = c(qq('@{names}.meth', collapse=F), qq('@{names}.unmeth', collapse=F)),
        func='sum')

    meth_cols <- grep('\\.meth$', colnames(avgs), value=T)
    meth <- avgs[, c('chrom', 'start', 'end', 'intervalID', meth_cols)] %>%
        gather('samp', 'meth', -(chrom:end), -intervalID)

    unmeth_cols <- grep('\\.unmeth$', colnames(avgs), value=T)
    unmeth <- avgs[, c('chrom', 'start', 'end', 'intervalID', unmeth_cols)] %>%
        gather('samp', 'unmeth', -(chrom:end), -intervalID)

    stopifnot(all(meth$start == unmeth$start) && all(meth$end == unmeth$end))

    avgs <- meth %>% bind_cols(unmeth %>% select(unmeth))

    avgs <- avgs %>% mutate(samp = gsub('\\.meth$', '', samp))

    avgs <- avgs %>%
        mutate(avg = meth / (meth + unmeth), cov = meth + unmeth)

    if ((!is.null(min_cov) || !is.null(min_samples)) && !pre_screen) {
        if (is.null(min_cov)){
            min_cov <- 1
            warning('Did not provide minimal coverage (min_cov). Setting to 1')
        }
        message(qq("Taking only intervals with coverage >= @{min_cov} in at least @{min_samples} samples"))
        avgs <- .gpatterns.filter_coverage(avgs,
                                           min_cov,
                                           min_samples,
                                           mask_by_cov)
    }

    if (!is.null(min_var) || !is.null(var_quantile)) {
        if (!is.null(var_quantile)) {
            vars <- avgs %>%
                group_by(chrom, start, end) %>%
                summarise(v = var(avg, na.rm=T)) %>% .$v
            min_var <- quantile(vars, probs = 1 - var_quantile, na.rm=T)
        }

        message(qq('Taking only intervals with variance >= @{round(min_var, digits = 2)}'))
        avgs <- avgs %>%
            group_by(chrom, start, end) %>%
            mutate(v = var(avg, na.rm=T)) %>%
            ungroup %>%
            filter(v >= min_var)
    }

    n_intervals <- distinct(avgs, chrom, start, end) %>% nrow
    message(qq('number of intervals: @{scales::comma(n_intervals)}'))


    return(avgs %>%
               select(chrom, start, end, intervalID, samp, meth, unmeth, avg, cov) %>%
               tbl_df())

}


.gpatterns.get_avg_meth_not_tidy <- function(tracks,
                                             intervals,
                                             iterator = NULL,
                                             min_var = NULL,
                                             var_quantile = NULL,
                                             names = NULL,
                                             file = NULL,
                                             intervals.set.out = NULL,
                                             rm_intervals = FALSE) {

    names <- names %||% tracks

    vtracks_pref <- .random_track_name()
    vtracks_meth <- paste0(vtracks_pref, '_', 1:length(tracks), '_meth')
    vtracks_unmeth <- paste0(vtracks_pref, '_', 1:length(tracks), '_unmeth')

    walk2(vtracks_meth, .gpatterns.meth_track_name(tracks), gvtrack.create, func='sum')
    walk2(vtracks_unmeth, .gpatterns.unmeth_track_name(tracks), gvtrack.create, func='sum')

    expr <- qqv('@{vtracks_meth} / ( @{vtracks_meth} + @{vtracks_unmeth} )')
    avgs <- gextract(expr, intervals=intervals, iterator=iterator, colnames=names, file=file, intervals.set.out = intervals.set.out)

    if (!is.null(file) || !is.null(intervals.set.out)){
        return(NULL)
    }

    if (!is.null(min_var) || !is.null(var_quantile)) {
        vars <- apply(avgs %>% select(-(chrom:end), -intervalID), 1, function(x) var(x, na.rm=T))
        if (!is.null(var_quantile)) {
            min_var <- quantile(vars, probs = 1 - var_quantile, na.rm=T)
        }

        message(qq('Taking only intervals with variance >= @{round(min_var, digits = 2)}'))
        vars <- avgs %>% filter(vars >= min_var)
    }

    n_intervals <- nrow(avgs)
    message(qq('number of intervals: @{scales::comma(n_intervals)}'))

    walk(c(vtracks_meth, vtracks_unmeth), gvtrack.rm)

    if (rm_intervals){
        gintervals.rm(intervals, force = TRUE)
    }

    return(avgs %>% tbl_df)
}

#' screen intervals by coverage in multiple tracks
#'
#' @param tracks methylation tracks
#' @param intervals genomic scope for which the function is applied
#' @param iterator see iterator in \code{\link[gextract]{misha}}. if NULL iterator
#' would be set to CpGs
#' @param min_cov minimal coverage for iterator interval
#' @param min_samples minimal number of samples with cov >= min_cov. if min_cov
#' is NULL it would be set to 1.
#' @param intervals.set.out intervals.set.out
#'
#' @return
#' @export
#'
#' @examples
gpatterns.screen_by_coverage <- function(tracks,
                                         intervals,
                                         iterator,
                                         min_cov,
                                         min_samples = NULL,
                                         intervals.set.out = NULL){
    cov_tracks <- .gpatterns.cov_track_name(tracks)
    min_samples <- min_samples %||% 1
    expr <- paste(qqv('(@{cov_tracks} >= min_cov)'), collapse = ', ')
    expr <- qq('sum(@{expr}, na.rm=T) >= @{min_samples}')
    intervs <- gscreen(expr,
                       intervals = intervals,
                       iterator = iterator,
                       intervals.set.out = intervals.set.out)
    return(intervs)
}

#' Filter intervals by number of CpGs
#'
#' @param intervals intervals to filter
#' @param min_cpgs minimal number of CpGs
#' @param cgs_track cpgs track (has '1' at every CpG in the genome and '0' otherwise)
#'
#' @return
#' @export
#'
#' @examples
gpatterns.filter_cpgs <- function(intervals,
                                  min_cpgs,
                                  cgs_track = .gpatterns.genome_cpgs_track) {
    gvtrack.create("cpg_count", cgs_track, func = "sum")
    intervals <- gscreen(qq('cpg_count >= @{min_cpgs}'),
                         intervals = intervals,
                         iterator = intervals)
    gvtrack.rm("cpg_count")
    return(intervals)
}


#' Filter intervals set by distance from another intervals set
#'
#' @param intervals intervals to filter
#' @param annot_intervs annotation intervals to filter by
#' @param annot_intervs_dist maximal distance from annot_intervs
#' @param orientation orientation (upstream, downstream or both)
#' @param use_tss use TSS as annot intervs
#' @param tss_dist annot_intervs_dist for TSS
#'
#' @return
#' @export
#'
#' @examples
gpatterns.filter_by_dist <- function(intervals,
                                     annot_intervs = NULL,
                                     annot_intervs_dist = 0,
                                     orientation = "upstream",
                                     use_tss = FALSE,
                                     tss_max_dist = 1000){

    stopifnot(orientation %in% c('both', 'upstream', 'downstream'))

    if (use_tss) {
        annot_intervs <- "intervs.global.tss"
        annot_intervs_dist <- tss_max_dist
    }

    if (!is.null(annot_intervs)) {
        nn <- gintervals.neighbors1(intervals, annot_intervs, na.if.notfound = T)
        if (orientation == "both") {
            if(is.null(nrow(annot_intervs))){
                message(sprintf("taking only intervals that are %s from %s",
                                annot_intervs_dist,
                                annot_intervs))
            }
            intervals <- nn %>%
                filter(abs(dist) <= annot_intervs_dist) %>%
                select_(.dots=colnames(intervals))
        } else if (orientation == "upstream") {
            if(is.null(nrow(annot_intervs))){
                message(sprintf("taking only intervals that are %s upstream to %s",
                                annot_intervs_dist,
                                annot_intervs))
            }
            intervals <- nn %>%
                filter(dist >= 0, abs(dist) <= annot_intervs_dist) %>%
                select_(.dots=colnames(intervals))

        } else if (orientation == "downstream") {
            if(is.null(nrow(annot_intervs))){
                message(sprintf("taking only intervals that are %s downstream to %s",
                                annot_intervs_dist,
                                annot_intervs))
            }
            intervals <- nn %>%
                filter(dist <= 0, abs(dist) <= annot_intervs_dist) %>%
                select_(.dots=colnames(intervals))
        }
    }

    return(intervals)
}

.gpatterns.filter_coverage <- function(avgs,
                                       min_cov,
                                       min_samples = NULL,
                                       mask_cpgs = FALSE){
    min_samples <- min_samples %||% 1
    avgs <- avgs %>%
        group_by(chrom, start, end) %>%
        mutate(n = sum(cov >= min_cov, na.rm=T)) %>%
        ungroup %>%
        filter(n >= min_samples)

    if (mask_cpgs){
        avgs <- avgs %>% mutate(avg = ifelse(cov >= min_cov, avg, NA))
    }
    return(avgs)
}

# Analysis Functions ------------------------------------------------


#' Cluster average methylation of multiple tracks
#'
#' @param avgs avgs pre-computed intervals set with average methylation, output of gpatterns.get_avg_meth
#' @param K number of clusters for kmeans++
#' @param verbose display kmeans++ messages
#' @param clust_order_func order clusters by clust_order_func [default: mean]
#' @param clust_columns cluster the columns (samples) using hclust
#' @param column_k cut the hclust of the columns and return an additional
#' 'samp_clust' column with the column cluster
#' @param ret_hclust return a list with avgs (the clustering result) and
#' hc (the hclust object)
#'
#' @return data frame with the following fields:
#' \itemize{
#'      \item{'chrom'}{}
#'      \item{'start'}{}
#'      \item{'end'}{}
#'      \item{'ord'}{order of intervals within the cluster}
#'      \item{'samp'}{sample id}
#'      \item{'samp_clust'}{id of the column cluster (hclust cutree)}
#'      \item{'clust'}{cluster id}
#'      \item{'original columns'}{original columns of avgs (meth,unmeth,avg,cov)}
#' }
#' @export
#'
#' @examples
gpatterns.cluster_avg_meth <- function(
    avgs,
    K,
    verbose = FALSE,
    clust_order_func = mean,
    clust_columns = TRUE,
    column_k = NULL,
    ret_hclust = FALSE) {

    avgs_mat <- avgs %>%
        select(chrom, start, end, samp, avg) %>%
        spread(samp, avg)
    clust <- avgs_mat %>%
        unite('id', chrom, start, end, sep='_') %>%
        TGL_kmeans(K=K,
                   verbose=verbose,
                   id_column=T) %>%
        separate(id, c('chrom', 'start', 'end')) %>%
        mutate_at(vars(start, end), as.numeric)

    avgs <- avgs %>%
        left_join(clust, by=c('chrom', 'start', 'end'))

    # order clusters
    avgs <- avgs %>%
        group_by(clust) %>%
        mutate(m = clust_order_func(avg)) %>%
        ungroup %>%
        group_by(chrom, start, end) %>%
        mutate(row_mean = mean(avg, na.rm=T)) %>%
        ungroup %>%
        arrange(m, row_mean) %>%
        select(-m, -row_mean)

    #order within cluster using hclust
    avgs <- avgs %>%
        distinct(clust) %>%
        mutate(new_clust = 1:(nrow(.))) %>%
        right_join(avgs, by='clust') %>%
        select(chrom, start, end, samp, clust=new_clust, meth, unmeth, avg, cov)
    avgs <- avgs %>%
        group_by(clust) %>%
        do(.hclust_order(.,
                         c('chrom', 'start', 'end'),
                         'samp',
                         'avg',
                         method='ward.D2')) %>%
        ungroup


    if (clust_columns) {
        avgs <- .gpatterns.cluster_columns(avgs, column_k = column_k)
    }

    if(ret_hclust){
        return(list(avgs=avgs, hc=hc))
    }

    return(avgs)
}

.gpatterns.cluster_columns <- function(avgs, column_k = NULL){
    avgs_mat <- avgs %>%
        select(chrom, start, end, clust, ord, samp, avg) %>%
        spread(samp, avg)
    dist_mat <-  avgs_mat %>%
        select(-(chrom:ord)) %>%
        as.matrix() %>%
        t %>%
        dist
    hc <- dist_mat %>%
        hclust(method = "ward.D2")
    hc <- gclus::reorder.hclust(hc, dist_mat)
    ord <- hc$order

    samples <- colnames(avgs_mat %>% select(-(chrom:ord)))

    if (!is.null(column_k)){
        ct <- cutree(hc, k=column_k)
        samp_clust <- data_frame(samp = as.character(names(ct)), samp_clust=ct) %>%
            left_join(data_frame(
                samp = samples[ord],
                samp_ord = 1:length(samples)), by='samp') %>%
            arrange(samp_ord)
        clust_ord <- samp_clust %>%
            distinct(samp_clust) %>%
            .$samp_clust
        samp_clust <- samp_clust %>%
            mutate(samp_clust = plyr::mapvalues(samp_clust, clust_ord, 1:max(clust_ord))) %>%
            select(-samp_ord)
        avgs <- avgs %>%
            left_join(samp_clust, by='samp') %>%
            select(chrom, start, end, ord, samp, samp_clust, clust, meth, unmeth, avg, cov)

    }
    avgs$samp <- factor(avgs$samp, levels=samples[ord])
    return(avgs)
}


# Plotting Functions ------------------------------------------------

#' Plots smoothScatter of average methylation of two tracks
#'
#' @param samples a vector of 2 sample names as present in 'samp' filed of avgs, or track names if
#' avgs is NULL
#' @param avgs pre-computed intervals set with average methylation, output of gpatterns.get_avg_meth.
#' if NULL, gpatterns.get_avg_meth(tracks=samples, ...) would be called. Note that
#' you using this option together with additional parameters to smoothScatter would
#' produce warnings.
#' @param fig_ofn figure filename. if NULL figure would be plotted on the current device
#' @param sample_names names for the two tracks to plot as the axis labels.
#'  if NULL the sample names would be used.
#' @param width plot width
#' @param height plot height
#' @param color.pal colo pallete
#' @param device function to open a device if fig.ofn is not NULL. default: 'png'
#' @param title_text title text
#' @param add_n add 'n = number_of_observations' to the title
#' @param ... other parameters for smoothScatter / gpatterns.get_avg_meth.
#'
#' @return None
#'
#' @examples
#'
#' @export
gpatterns.smoothScatter <- function(
    samples,
    avgs = NULL,
    fig_ofn = NULL,
    sample_names = NULL,
    width = 300,
    height = 300,
    color.pal = .smooth_scatter_pal2,
    device = "png",
    title_text = NULL,
    add_n = TRUE,
    ...) {

    avgs <- avgs %||% .do.call_ellipsis(gpatterns.get_avg_meth, list(tracks=samples, tidy=TRUE), ...)


    if (length(samples) != 2) {
        stop("Please provide 2 samples")
    }
    d <- avgs %>% filter(samp %in% samples) %>% select(chrom, start, end, samp, avg) %>% spread(samp, avg)

    sample_names <- sample_names %||% samples

    stopifnot(length(sample_names) == 2)

    if (!is.null(fig_ofn)) {
        do.call(device, list(fig.ofn, width = width, height = height))
    }

    if (add_n){
        n_str <- paste0("n = ", scales::comma(nrow(d)))
        if (!is.null(title_text)){
            title_text <- qq('@{title_text}\n@{n_str}')
        } else {
            title_text <- n_str
        }
    }

    smoothScatter(d[[samples[1]]], d[[samples[2]]], colramp = color.pal, xlab = sample_names[1],
                  ylab = sample_names[2], main = title_text, ...)

    abline(0,1,lty=2)
    grid()
    if (!is.null(fig_ofn)) {
        dev.off()
    }
}



# Utility Functions ------------------------------------------------

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

.do.call_ellipsis <- function(f, additonal_params, ...){
    f_args <- names(as.list(args(f)))
    elipsis <- list(...)
    elipsis <- map(names(elipsis), function(x) if (x %in% f_args) { return(elipsis[[x]]) } )
    do.call(f, c(additonal_params, elipsis))
}

.hclust_order <- function(d, keys, variable, value, ...){
    d_mat <- d %>%
        select(one_of(c(keys, variable, value))) %>%
        spread_(variable, value)
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
