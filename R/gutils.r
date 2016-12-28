########################################################################
.random_track_name <- function(len=12){
    stringi::stri_rand_strings(1, len, "[a-z]")
}

########################################################################
#' Finds neighbors between two sets of intervals (and does not return conflicting column names)
#'
#' @inheritParams misha::gintervals.neighbours
#' @param suffix suffix for conflicting colnames
#'
#' @return
#' @export
#'
#' @examples
gintervals.neighbors1 <- function(intervals1 = NULL,
                                  intervals2 = NULL,
                                  suffix = '1',
                                  maxneighbors = 1,
                                  mindist = -1e+09,
                                  maxdist = 1e+09,
                                  na.if.notfound = FALSE) {

    res <-
        gintervals.neighbors(
            intervals1 = intervals1,
            intervals2 = intervals2,
            maxneighbors = maxneighbors,
            mindist = mindist,
            maxdist = maxdist,
            na.if.notfound = na.if.notfound
        )

    intervs2_col_start <- which(colnames(res) == "chrom")[2]
    intervals1.cols <- 1:(intervs2_col_start - 1)
    intervals2.cols <- intervs2_col_start:ncol(res)
    conflict.names <- which(colnames(res)[intervals2.cols] %in% colnames(res)[intervals1.cols])
    colnames(res)[(intervs2_col_start - 1) + conflict.names] <- paste0(colnames(res)[(intervs2_col_start -
        1) + conflict.names], suffix)
    return(res %>% tbl_df)
}

########################################################################
#' Filter intervals set by distance to another.
#' A wrapper around gintervals.neighbours and filter(dist <= max_distance)
#'
#' @param intervals1 intervals set
#' @param intervals2 intervals set by which to filter intervals1
#' @param max_distance maximal distance of every interval in intervals1 from intervals2 (defualt 0)
#' @param abs_dist take the absolute distance
#' @param ... additional parameters to gintervals.neighbours1
#'
#' @return
#' @export
#'
#' @examples
gintervals.filter <- function(intervals1, intervals2, max_distance=0, abs_dist = TRUE, ...){
    intervals1_cols <- colnames(intervals1)
    res <- intervals1 %>% gintervals.neighbors1(intervals2, ...)
    if (abs_dist){
        res$dist <- abs(res$dist)
    }
    res <- res %>% filter(dist <= max_distance) %>% select(one_of(intervals1_cols))
    return(res)
}

########################################################################
gintervals.neighbors_y <- function(intervals1, intervals2, maxneighbors=1, mindist=-1e+09, maxdist=1e+09, na.if.notfound=FALSE)
{
    cols1 <- colnames(intervals1)
    cols2 <- colnames(intervals2)
    collisions <- cols2 %in% cols1
    cols2[collisions] <- paste0(cols2[collisions], '.y')
    cols12 <- c(cols1, cols2)

    if (nrow(intervals2) > 0) {
        neighbors <- gintervals.neighbors(intervals1, intervals2, maxneighbors=maxneighbors,
                                          mindist=mindist, maxdist=maxdist, na.if.notfound=na.if.notfound)
        colnames(neighbors)[1:length(cols12)] <- cols12
    }
    else {
        neighbors <- bind_cols(intervals1, data.frame(matrix(NA, nrow=nrow(intervals1), ncol=ncol(intervals2)+1)))
        colnames(neighbors) <- c(cols12, 'dist')
    }

    return(neighbors)
}

########################################################################
gintervals.left_join <- function(intervals1, intervals2, maxneighbors=1, mindist=-1e+09, maxdist=1e+09)
{
    joined <- gintervals.neighbors_y(intervals1, intervals2, maxneighbors=maxneighbors, mindist=mindist, maxdist=maxdist, na.if.notfound=TRUE)
    joined <- joined %>%
              group_by(chrom, start, end) %>%
              filter(is.na(dist) | (abs(dist) == min(abs(dist)))) %>%
              ungroup()

    return(joined)
}

########################################################################
gintervals.centers <- function(inv) {
    inv %>%
        mutate(center = floor((start + end) / 2), start=center, end=center+1) %>%
        select(-center)
}

########################################################################
gintervals.expand <- function(inv, expansion = 100) {
    inv %>%
        mutate(start = start - expansion, end = end + expansion) %>%
        as.data.frame %>%
        gintervals.force_range()
}

########################################################################
gintervals.normalize <- function(inv, size) {
    centers <- gintervals.centers(inv) %>%
        mutate(end = end - 1)
    return(gintervals.expand(centers, floor(size/2)))
}

########################################################################
gintervals.distance <- function(start1,end1,start2,end2) {
    pmax(pmax(start1, start2) - pmin(end1, end2), 0)
}

########################################################################
#' Creates a virtual track and runs gextract
#'
#' @param tracks tracks
#' @param intervals intervals
#' @param colnames colnames
#' @param iterator iterator
#' @param band band
#' @param file file
#' @param intervals.set.out intervals.set.out
#' @param func func
#' @param params params
#'
#' @return
#' @export
#'
#' @examples
gvextract <- function(tracks, intervals, colnames = NULL, iterator = NULL,
              band = NULL, file = NULL, intervals.set.out = NULL, func=NULL, params=NULL){
    vtracks_pref <- .random_track_name()
    vtracks <- paste0(vtracks_pref, '_', 1:length(tracks))

    walk2(vtracks, tracks, gvtrack.create, func=func, params=params)

    if (is.null(colnames)){
        colnames <- tracks
    }

    res <- gextract(vtracks,
                    intervals = intervals,
                    iterator = iterator,
                    colnames=colnames,
                    band=band,
                    file=file,
                    intervals.set.out=intervals.set.out)

    walk(vtracks, gvtrack.rm)

    return(res)
}


########################################################################
#' Returns the result of track expressions evaluation for each of the
#' iterator intervals, and cbinds the intervals (instead of intervalID)
#'
#' @inheritParams misha::gextract
#' @param suffix suffix for conflicting column names
#'
#' @return
#' @export
#'
#' @seealso \code{\link{misha::gextract}}
#' @examples
gextract.left_join <- function(expr, intervals = NULL, colnames = NULL, iterator = NULL, band = NULL, file = NULL, intervals.set.out = NULL, suffix='1'){
    d <- gextract(expr, intervals = intervals, colnames = colnames, iterator = iterator, band = band, file = file, intervals.set.out = intervals.set.out)
    conflict_names <- which(colnames(intervals) %in% colnames(d))
    colnames(intervals)[conflict_names] <- paste0(colnames(intervals)[conflict_names], suffix)
    intervals$intervalID <- 1:nrow(intervals)
    d <- d %>% arrange(intervalID) %>% left_join(intervals, by='intervalID') %>% select(-intervalID)
    return(d)
}

########################################################################
gintervals.which=function(expr, intervals, iterator=expr, which.func=max) {
  result.intervals = intervals
  f = function(x) {
    v <- which.func(x, na.rm=T)
    idx <- match(v, x)
    result.intervals[GAPPLY.INTERVID, ] <<- GAPPLY.INTERVALS[idx, ]
    return(v)
  }
  gintervals.apply("f(x1)", expr, intervals, enable.gapply.intervals=T, iterator=iterator)
  result.intervals
}

########################################################################
gmax.coord=function(chr) {
  x = ALLGENOME[[1]]
  x[x$chrom == chr,"end"]
}

########################################################################
gseq.extract_conv <- function(..., methylated=TRUE) {
    res <- toupper(gseq.extract(...))
    if (methylated) {
        return(gsub('C(?!G)', 'T', res, perl=T))
    } else {
        return(gsub('C', 'T', res))
    }
}

########################################################################
gseq.rev_comp <- function(s) {
    chartr('acgtACGT', 'tgcaTGCA', s) %>% stringi::stri_reverse()
}

########################################################################
sum_ignore_na <- purrr::partial(sum, na.rm=T)
sum_ignore_na_all <- sum_ignore_na

########################################################################
#' split genomic intervals to bins
#'
#' @param intervals intervals
#' @param nbins number of bins
#'
#' @return intervals set with \code{nbins} rows. if \code{nbins} < number of
#' intervals the original intervals would be returned
#'
#' @export
#'
#' @examples
gbin_intervals <- function(intervals, nbins){
    if (nbins < nrow(intervals)){
        warning('nbins is smaller than the number of intervals, would return intervals unchanged')
    }
    chrs <-  intervals %>% mutate(bins_chr=pmax(round(end / sum(end) * nbins), 1))
    bin_diff <- sum(chrs$bins_chr) - nbins
    if (bin_diff > 0){
        chrs <- chrs %>% arrange(-bins_chr)
        for (i in 1:abs(bin_diff)){
            chrs$bins_chr[i] <- max(1, chrs$bins_chr[i] - 1)
        }
    }
    if (bin_diff < 0){
        chrs <- chrs %>% arrange(bins_chr)
        for (i in 1:abs(bin_diff)){
            chrs$bins_chr[i] <- chrs$bins_chr[i] + 1
        }
    }
    binned_intervals <- plyr::adply(chrs, 1, function(x)
        tibble(start=round(seq(x$start, x$end, length.out=x$bins_chr + 1))) %>%
        mutate(end = ifelse(lead(start) == last(start), lead(start), lead(start) - 1)) %>%
        drop_na
        ) %>% select(-bins_chr)
    return(binned_intervals)
}


########################################################################
#' Return an alternate empirical cumulative distribution, counting only
#' values with a lower rank.
#'
#' For every value of x, return the number of values with a strictly
#' lower rank than x.
cume_dist_min <- function(x) {1 - cume_dist(-x)}


########################################################################
gintervals.mark_overlapping <- function(intervals, unify_touching_intervals=TRUE, var='overlap')
{
    canonic <- intervals %>%
               gintervals.force_range() %>%
               gintervals.canonic(unify_touching_intervals)
    return(intervals %>% mutate(overlap=attr(canonic, 'mapping')))
}

########################################################################
#' Runs R commands on a cluster that supports SGE
#'
#' @inheritParams misha::gcluster.run
#' @param command_list list of strings with R commands
#' @param packages packages to load for each command
#' @param jobs_title title of job names. names would have the title followed
#' by a serial number
#' @param job_names vector with the names of the jobs
#' @param collapse_results collapse return values of the jobs to a data frame.
#' if not possible - would return the usual list.
#'
#' @param memory memory requirments (would be called using \code{memory_flag})
#' @param threads threads requirments (would be called using \code{threads_flag})
#' @param io_saturation io_saturation requirments (would be called using \code{io_saturation_flag})
#' @param memory_flag flag for memory requirment (formatted as in qq)
#' @param threads_flag flag for threads requirment (formatted as in qq)
#' @param io_saturation_flag flag for io_saturation requirment (formatted as in qq)
#'
#' @return if collapse_results is TRUE: data frame with the results of all jobs (rbinded).
#' if collapse_results is FALSE returns the same as: \link[misha]{gcluster.run}
#'
#' @export
#' @seealso \code{\link{misha::gcluster.run}}
#' @examples
#'
gcluster.run2 <- function (...,
                           command_list = NULL,
                           opt.flags = "",
                           max.jobs = 400,
                           debug = FALSE,
                           R = paste0(R.home(component='bin'), '/R'),
                           packages = NULL,
                           jobs_title = NULL,
                           job_names = NULL,
                           collapse_results = FALSE,
                           memory = NULL,
                           threads = NULL,
                           io_saturation = NULL,
                           memory_flag = '-l mem_free=@{memory}G',
                           threads_flag = '-pe threads @{threads}',
                           io_saturation_flag = '-l io_saturation=@{io_saturation}',
                           script = paste(get(".GLIBDIR"), "exec", "sgjob.sh", sep = "/")){

    if (!is.null(command_list)){
        if (!is.null(packages)){
            packages_command <- paste(sprintf("library('%s');", packages), collapse=' ')
            commands <- purrr::map(command_list, function(x) paste('{ ', packages_command, x, ' }'))
        }
        commands <- purrr::map(command_list, function(x) parse(text=x))
    } else {
        commands <- as.list(substitute(list(...))[-1L])
    }

    if (!is.null(memory)){
        opt.flags <- paste(opt.flags, qq(memory_flag))
    }
    if (!is.null(threads)){
        opt.flags <- paste(opt.flags, qq(threads_flag))
    }
    if (!is.null(io_saturation)){
        opt.flags <- paste(opt.flags, qq(io_saturation_flag))
    }

    if (length(commands) < 1)
        stop("Usage: gculster.run2(..., command_list = NULL, opt.flags = \"\" max.jobs = 400, debug = FALSE)",
            call. = F)
    if (!length(system("which qsub", ignore.stderr = T, intern = T)))
        stop("gcluster.run2 must run on a host that supports Sun Grid Engine (qsub)",
            call. = F)
    .gcheckroot()
    tmp.dirname <- ""
    submitted.jobs <- c()
    tryCatch({
        tmp.dirname <- tempfile(pattern = "", tmpdir = paste(get("GROOT"),
            "/tmp", sep = ""))
        if (!dir.create(tmp.dirname, recursive = T, mode = "0777"))
            stop(sprintf("Failed to create a directory %s", tmp.dirname),
                call. = F)
        cat("Preparing for distribution...\n")
        save(.GLIBDIR, file = paste(tmp.dirname, "libdir", sep = "/"))
        vars <- ls(all.names = TRUE, envir = parent.frame())
        envir <- parent.frame()
        while (!identical(envir, .GlobalEnv)) {
            envir <- parent.env(envir)
            vars <- union(vars, ls(all.names = TRUE, envir = envir))
        }
        suppressWarnings(save(list = vars, file = paste(tmp.dirname, "envir",
            sep = "/"), envir = parent.frame()))
        .GSGECMD <- commands
        save(.GSGECMD, file = paste(tmp.dirname, "commands",
            sep = "/"))
        opts <- options()
        save(opts, file = paste(tmp.dirname, "opts", sep = "/"))
        cat("Running the commands...\n")
        completed.jobs <- c()
        progress <- -1
        repeat {
            num.running.jobs <- length(submitted.jobs) - length(completed.jobs)
            if (length(submitted.jobs) < length(commands) &&
                num.running.jobs < max.jobs) {
                istart <- length(submitted.jobs) + 1
                iend <- min(length(commands), istart + (max.jobs -
                  num.running.jobs) - 1)
                for (i in istart:iend) {
                  out.file <- sprintf("%s/%d.out", tmp.dirname,
                    i)
                  err.file <- sprintf("%s/%d.err", tmp.dirname,
                    i)
                  script <- paste(get(".GLIBDIR"), "exec", "sgjob.sh",
                    sep = "/")
                    if (!is.null(job_names)){
                        job.name <-job_names[i]
                    } else if (!is.null(jobs_title)){
                        job.name <- sprintf('%s_%s', jobs_title, i)
                    } else {
                        job.name <- sprintf('sgjob_%s', i)
                    }
                  command <- sprintf("qsub -terse -cwd -S /bin/bash -N %s -o %s -e %s -V %s %s %d '%s' '%s'",
                    job.name, out.file, err.file, opt.flags, script, i,
                    tmp.dirname, R)
                  jobid <- system(command, intern = TRUE)
                  if (length(jobid) != 1)
                    stop("Failed to run qsub", call. = FALSE)
                  if (debug)
                    cat(sprintf("\tSubmitted job %d (id: %s)\n",
                      i, jobid))
                  submitted.jobs <- c(submitted.jobs, jobid)
                }
            }
            Sys.sleep(3)
            running.jobs <- .gcluster.running.jobs(submitted.jobs)
            old.completed.jobs <- completed.jobs
            completed.jobs <- setdiff(submitted.jobs, running.jobs)
            if (debug) {
                delta.jobs <- setdiff(completed.jobs, old.completed.jobs)
                if (length(delta.jobs) > 0) {
                  for (jobid in delta.jobs) cat(sprintf("\tJob %d (id: %s) completed\n",
                    match(jobid, submitted.jobs), jobid))
                }
                if (!length(running.jobs) && length(submitted.jobs) ==
                  length(commands))
                  break
                new.progress <- length(completed.jobs)
                if (new.progress != progress) {
                  progress <- new.progress
                  cat(sprintf("\t%d job(s) still in progress\n",
                    length(commands) - progress))
                }
            }
            else {
                if (!length(running.jobs) && length(submitted.jobs) ==
                  length(commands))
                  break
                new.progress <- as.integer(100 * length(completed.jobs)/length(commands))
                if (new.progress != progress) {
                  progress <- new.progress
                  cat(sprintf("%d%%...", progress))
                }
                else cat(".")
            }
        }
        if (!debug && progress != -1 && progress != 100)
            cat("100%\n")
    }, interrupt = function(interrupt) {
        cat("\n")
        stop("Command interrupted!", call. = FALSE)
    }, finally = {
        if (length(submitted.jobs) > 0) {
            running.jobs <- .gcluster.running.jobs(submitted.jobs)
            answer <- c()
            for (i in 1:length(commands)) {
                res <- list()
                res$exit.status <- NA
                res$retv <- NA
                res$stdout <- NA
                res$stderr <- NA
                if (submitted.jobs[i] %in% running.jobs)
                  res$exit.status <- "interrupted"
                else {
                  fname <- sprintf("%s/%d.retv", tmp.dirname,
                    i)
                  if (file.exists(fname)) {
                    load(fname)
                    res$exit.status <- "success"
                    res$retv <- retv
                  }
                  else res$exit.status <- "failure"
                }
                out.file <- sprintf("%s/%d.out", tmp.dirname,
                  i)
                if (file.exists(out.file)) {
                  f <- file(out.file, "rc")
                  res$stdout <- readChar(f, 1000000)
                  close(f)
                }
                err.file <- sprintf("%s/%d.err", tmp.dirname,
                  i)
                if (file.exists(err.file)) {
                  f <- file(err.file, "rc")
                  res$stderr <- readChar(f, 1000000)
                  close(f)
                }
                answer[[i]] <- res
            }
            for (job in running.jobs) system(sprintf("qdel %s",
                job), ignore.stderr = T, intern = T)
            unlink(tmp.dirname, recursive = TRUE)

            if (collapse_results){
                canswer <- tryCatch(
                    purrr::map_df(answer, function(x) x$retv),
                    error = function(e) {
                        message('returning original output due to an error. collapse your reults manually (are all the parts data frames?)'
                        )
                        return(NULL)
                    })

                if (!is.null(canswer)) {
                    return(canswer)
                }
            }

            return(answer)
        }
        unlink(tmp.dirname, recursive = TRUE)
    })
}



########################################################################
#' Wrapper around gtrack.import_mappedseq for (multiple) bam files
#'
#' @param bam_files vector of bam (or sam) files to import
#' @param ... parameters of gtrack.import_mappedseq
#'
#' @return
#' @export
#' @seealso \code{\link{misha::gtrack.import_mappedseq}}
#' @examples
gtrack.import_mappedseq_bam <- function(bam_files, ...){
    cat_cmd <- 'cat'
    if (length(bam_files) > 1){
        cat_cmd <- 'samtools cat'
    }
    files <- paste(bam_files, collapse=' ')
    tmp_fifo <- tempfile()
    tryCatch({
        system(qq('mkfifo @{tmp_fifo}; @{cat_cmd} @{files} | samtools view -h > @{tmp_fifo}'), wait=FALSE)
        gtrack.import_mappedseq(file=tmp_fifo, ...)
        }, finally=system(qq('rm -f @{tmp_fifo}')))
}
