# import Functions ------------------------------------------------
#' Creates a track from bam files.
#'
#' @description Creates a track from bam files for a specific sample.
#' Using this method is discouraged: if you have large dataset use the
#' Snakemake pipeline instead.
#'
#' @param bams bams
#' @param workdir workdir
#' @param track track
#' @param description description
#' @param conversion conversion
#' @param paired_end paired_end
#' @param downsample downsample
#' @param dsn dsn
#' @param pat_cov_tracks pat_cov_tracks
#' @param pat_freq_len pat_freq_len
#'
#' @return
#' @export
#'
#' @examples
gpatterns.import_from_bam <- function(bams,
                                      workdir,
                                      track,
                                      description,
                                      conversion = 'CT',
                                      paired_end = TRUE,
                                      downsample = FALSE,
                                      dsn = NULL,
                                      pat_cov_tracks = c(3,5,7),
                                      pat_freq_len = 2){

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
#' needs to have the folloing fields: fid,pattern
#' @param add_read_id save the read_id together with the patterns
#' @param noise_threshold threshold to consider pattern as 'noise'
#' @param overwrite overwrite existing tracks
#' @param canonize convert to canonic form (see: misha::gintervals.canonic).
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


    # Create a table mapping pattern positions to the fid
    loci_tab <- pat_space %>%
        group_by(fid, chrom) %>%
        summarize(start=min(start), end=max(start)+1) %>%
        ungroup() %>%
        select(chrom, start, end, fid)

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


