% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge-tidy_cpgs.R
\name{gpatterns.merge_tidy_cpgs}
\alias{gpatterns.merge_tidy_cpgs}
\title{merge tidy cpgs directories}
\usage{
gpatterns.merge_tidy_cpgs(
  dirs,
  out_dir = tempdir(),
  nbins = nrow(gintervals.all()),
  paired_end = TRUE,
  stats_dir = paste0(out_dir, "/stats"),
  filter_dups_bin = system.file("import", "filter_dups_cpgs.py", package = "gpatterns")
)
}
\arguments{
\item{dirs}{tidy_cpgs directories}

\item{out_dir}{output directory}

\item{nbins}{number of genomic bins for output tidy_cpgs}

\item{paired_end}{is the run paired end}

\item{stats_dir}{directory for tidy_cpgs stats}

\item{filter_dups_bin}{binary for filter_dups_cpgs python script}
}
\value{
NULL#'
}
\description{
merge tidy cpgs directories
}
