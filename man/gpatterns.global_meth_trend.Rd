% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/avg_meth.R
\name{gpatterns.global_meth_trend}
\alias{gpatterns.global_meth_trend}
\title{Plot global methylation stratified on other tracks}
\usage{
gpatterns.global_meth_trend(
  tracks,
  strat_track = .gpatterns.cg_cont_500_track,
  strat_breaks = seq(0, 0.08, by = 0.002),
  intervals = .gpatterns.genome_cpgs_intervals,
  iterator = .gpatterns.genome_cpgs_intervals,
  min_cov = NULL,
  min_cgs = NULL,
  names = NULL,
  groups = NULL,
  group_name = NULL,
  include.lowest = TRUE,
  ncol = 2,
  nrow = 2,
  width = 600,
  height = 560,
  fig_fn = NULL,
  xlab = strat_track,
  ylim = c(0, 1),
  title = "",
  legend = TRUE,
  colors = NULL,
  parallel = getOption("gpatterns.parallel")
)
}
\arguments{
\item{tracks}{tracks to plot}

\item{strat_track}{track to stratify average methylation by. default is CG content}

\item{strat_breaks}{breaks to determine the bins of strat_track}

\item{intervals}{genomic scope for which the function is applied}

\item{iterator}{track expression iterator (of both tracks and strat_track)}

\item{min_cov}{minimal coverage of each track}

\item{min_cgs}{minimal number of CpGs per bin}

\item{names}{alternative names for the track}

\item{groups}{a vector the same length of \code{tracks} with group for each track. Each group will on a different facet.}

\item{group_name}{name of the grouping variable (e.g. tumor, sample, patient, experiment)}

\item{include.lowest}{if 'TRUE', the lowest value of the range determined by breaks is included}

\item{ncol}{number of columns}

\item{nrow}{number of rows}

\item{width}{plot width (if fig_fn is not NULL)}

\item{height}{plot height (if fig_fn is not NULL)}

\item{fig_fn}{output filename for the figure (if NULL, figure would be returned)}

\item{xlab}{label for the x axis}

\item{ylim}{ylim of the plot}

\item{title}{title for the plot}

\item{legend}{add legend}

\item{colors}{custom colors}

\item{parallel}{get trends parallely}
}
\value{
list with trend data frame (under 'trend') and the plot (under 'p')
}
\description{
calculates the average methylation \code{(m / m + um)} in each
bin of \code{strat_track} and plots it. By default, plots the average methylation
in different bins of CpG content. This can be used as a sanity check for methylation
data - in general, methylation is high for regions with low CpG density,
and low for CpG dense regions (e.g. CpG islands).
}
\examples{

}
