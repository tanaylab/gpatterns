# Tracks and table names ------------------------------------------------
.gpatterns.downsampled_track_name <- function(track, dsn) {qqv('@{track}.ds@{dsn}') }
.gpatterns.cov_track_name <- function(track) { qqv('@{track}.cov') }
.gpatterns.avg_track_name <- function(track) { qqv('@{track}.avg') }
.gpatterns.meth_track_name <- function(track) { qqv('@{track}.meth') }
.gpatterns.unmeth_track_name <- function(track) { qqv('@{track}.unmeth') }
.gpatterns.pat_cov_track_name <- function(track, pat_len) { qqv('@{track}.pat@{pat_len}') }
.gpatterns.frag_cov_track_name     <- function(track) { qqv('@{track}.frag_cov')}
.gpatterns.fid_track_name     <- function(track) { qqv('@{track}.fid')}
.gpatterns.ncpg_track_name    <- function(track) { qqv('@{track}.ncpg')}
.gpatterns.n_track_name       <- function(track) { qqv('@{track}.n')}
.gpatterns.n0_track_name      <- function(track) { qqv('@{track}.n0')}
.gpatterns.n1_track_name      <- function(track) { qqv('@{track}.n1')}
.gpatterns.nx_track_name      <- function(track) { qqv('@{track}.nx')}
.gpatterns.nc_track_name      <- function(track) { qqv('@{track}.nc')}
.gpatterns.pat_meth_track_name    <- function(track) { qqv('@{track}.pat_meth')}
.gpatterns.pat_space_intervs_name    <- function(track) { qqv('@{track}.pat_space')}
.gpatterns.epipolymorphism_track_name <- function(track) { qqv('@{track}.epipoly')}

.gpatterns.patterns_track_names <- function(track) {paste0(track, '.', c('n', 'n0', 'n1', 'nx', 'nc', 'pat_meth', 'epipoly', 'fid'))}

.gpatterns.mix_qval_track_name <- function(track) { qqv('@{track}.mix_qval')}
.gpatterns.mix_meth_track_name <- function(track) { qqv('@{track}.mix_meth')}
.gpatterns.mix_unmeth_track_name <- function(track) { qqv('@{track}.mix_unmeth')}
.gpatterns.gain_meth_track_name <- function(track) { qqv('@{track}.gain_meth')}
.gpatterns.gain_unmeth_track_name <- function(track) { qqv('@{track}.gain_unmeth')}
.gpatterns.loss_meth_track_name <- function(track) { qqv('@{track}.loss_meth')}
.gpatterns.loss_unmeth_track_name <- function(track) { qqv('@{track}.loss_unmeth')}
.gpatterns.center_meth_track_name <- function(track) { qqv('@{track}.center_meth_ones')}
.gpatterns.center_unmeth_track_name <- function(track) { qqv('@{track}.center_unmeth_ones')}
.gpatterns.center_uni_track_name <- function(track) { qqv('@{track}.center_uni_ones')}
.gpatterns.gain_uni_track_name <- function(track) { qqv('@{track}.gain_uni')}
.gpatterns.loss_uni_track_name <- function(track) { qqv('@{track}.loss_uni')}

.gpatterns.patterns_tab_name  <- function(track) { 'patterns'}
.gpatterns.patterns_file_name <- function(track) { file.path(.gpatterns.base_dir(track), 'patterns.RData')}
.gpatterns.fids_tab_name      <- function(track) { 'fids' }
.gpatterns.fids_file_name     <- function(track) { file.path(.gpatterns.base_dir(track), 'fids.RData') }
.gpatterns.stats_file_name    <- function(track) { qq('@{.gpatterns.base_dir(track)}/stats.tsv') }

.gpatterns.tidy_cpgs_files <- function(track) {
    list.files(paste0(.gpatterns.base_dir(track), '/tidy_cpgs'), full.names=TRUE, pattern='tcpgs.gz')
}

.gpatterns.tidy_cpgs_dir <- function(track){
    paste0(.gpatterns.base_dir(track), '/tidy_cpgs')
}

.gpatterns.bipolar_model_tab_name     <- function(track) { 'mix' }
.gpatterns.bipolar_model_file_name     <- function(track) {
    file.path(.gpatterns.base_dir(track), 'bipolar.RData')
}
.gpatterns.bipolar_model_stats <- c('qval', 'mix.meth', 'mix.unmeth', 'center.meth', 'center.meth.ones', 'gain.meth', 'loss.meth', 'center.unmeth', 'center.unmeth.ones','gain.unmeth', 'loss.unmeth', 'center.uni', 'center.uni.ones', 'gain.uni', 'loss.uni')

.gpatterns.genome_cpgs_track <- 'seq.CG'
.gpatterns.genome_cpgs_intervals <- 'intervs.global.seq_CG'
.gpatterns.genome_next_cpg_intervals <- 'intervs.global.next_CG'
.gpatterns.cg_cont_500_track <- 'seq.CG_500_mean'



.gpatterns.special_intervals <- function(name){
    intervs_map <- list(tss = 'intervs.global.tss',
                        exon = 'intervs.global.exon',
                        utr3 = 'intervs.global.utr3',
                        intron = 'intervs.global.introns',
                        cgi = 'intervs.global.cgi_ucsc')
    if (name %in% names(intervs_map)){
        return(intervs_map[[name]])
    } else {
        stop(qq('interval @{name} does not exist'))
    }
}

.gpatterns.get_intervals <- function(intervals){
    if (is.character(intervals)){
        if (!gintervals.exists(intervals)){
            intervals <- .gpatterns.special_intervals(intervals)
        }
    }
    return(intervals)
}


# Color palettes ------------------------------------------------
# .blue_red_pal <- colorRampPalette(c("#87FFFF", "black", "#FF413D"))(1000)

#' @export
.blue_red_pal <- colorRampPalette(c("#00688B", "white", "#FF413D"))(1000)
.blue_black_red_yellow_pal <- colorRampPalette(c("white", "blue",  "black", "red", "yellow"))(1000)
.red_blue_pal <- rev(colorRampPalette(c("#87FFFF", "black", "#FF413D"))(1000))
.blue_red_yellow_pal <- colorRampPalette(c("white", "blue", "red", "yellow"))(1000)
.smooth_scatter_pal2 <- colorRampPalette(c("white", "white", "deepskyblue4", "gray",
    "darkgray", "black", "brown"))
.smooth_scatter_pal <- colorRampPalette(c("white", "white", "darkgrey", "black",
    "#FF413D", "yellow"))
.smooth_scatter_pal3 <- colorRampPalette(c("white", "blue", "red", "yellow", "black"))
