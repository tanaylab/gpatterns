#' Define putative enhancer regions from chip-seq tracks
#'
#' @param min_tss_dist minimal (absolute) distance from TSS
#' @inheritParams gpatterns.chip2peaks
#'
#' @return intervals set with putative enhancers
#' @export
#'
#' @examples
gpatterns.putative_enhancers <- function(chip_tracks,
                                             quant_thresh=0.999,
                                             normalize=NULL,
                                             min_tss_dist=2000,
                                             logical_gate = '|'){
    enh <- gpatterns.chip2peaks(chip_tracks=chip_tracks, quant_thresh=quant_thresh, normalize=normalize, logical_gate=logical_gate)

    if (!is.null(min_tss_dist)){
        enh <- enh %>% gintervals.neighbors1('intervs.global.tss') %>% filter(abs(dist) >= min_tss_dist) %>% select(chrom, start, end)
    }

    return(enh)
}

#' Calculate peaks of chip-seq tracks
#'
#' @param chip_tracks names of chip-seq tracks
#' @param quant_thresh quantile of chip signal that would be considered as peak
#' @param normalize normalize all the regions to have the same size.
#' if NULL no normalization would be done.
#' @param logical_gate logical function to combine the different tracks, e.g. '|' or '&'
#'
#' @return intervals set with peaks
#' @export
#'
#' @examples
gpatterns.chip2peaks <- function(chip_tracks,
                                 quant_thresh=0.999,
                                 normalize=NULL,
                                 logical_gate = '|'){
    vtracks <- paste0('v_', chip_tracks)
    walk2(vtracks, chip_tracks, function(vt, t) gvtrack.create(vt, t, 'global.percentile.max'))
    on.exit(walk(vtracks, gvtrack.rm))

    expr <- paste(sprintf("%s >= %s", vtracks, quant_thresh), collapse=qq(' @{logical_gate} '))    
    peaks <- gscreen(expr)

    if (!is.null(normalize)){
        peaks <- gintervals.normalize(peaks, normalize)
    }


    return(as_tibble(peaks))
}

#' Plot chip-seq singal with quantile lines
#'
#' @param chip_track chip-seq track
#' @param intervals intervals to plot at
#' @param quantiles line quantiles to plot
#'
#' @export
gpatterns.plot_chip_quantiles <- function(chip_track,
                                          intervals,
                                          quantiles=c(0.8, 0.9, 0.95, 0.99, 0.995, 0.999)){
    vtrack <- paste0('v', chip_track)
    gvtrack.create(vtrack, chip_track, 'global.percentile.max')
    on.exit(gvtrack.rm(vtrack))
    quant_threshes <- gquantiles(vtrack, percentiles=quantiles)
    sig <- gextract(vtrack, intervals=intervals, colnames='signal')
    quant_threshes <- tibble(thresh=quant_threshes, quant=factor(quantiles))

    sig %>% ggplot(aes(x=start, y=-log2(1 - signal))) + geom_col() + geom_hline(data=quant_threshes, aes(yintercept=-log2(1-thresh), color=quant)) + scale_color_discrete(name='threshold\nquantile') + scale_x_continuous(labels=comify) + ylab(qq('-log2(1 - @{chip_track})'))
}
