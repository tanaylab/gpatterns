# Bissli Functions ------------------------------------------------

#' Align bisulfite converted reads to a reference genome using bowtie2
#'
#' @param r1_fastq a vector of FASTQ files to align (read1)
#' @param out_bam bam output file with aligned reads
#' @param genome_seq A directory holding the reference genome, one FASTA file per
#' chromosome
#' @param bissli2_idx The basename of the index to be searched. The index is created using
#' gpatterns.bissli2_build
#' @param r2_fastq a vector of FASTQ files to align (read2)
#' @param bissli2_bin path for bissli2-align
#' @param bowtie2 path for bowtie2
#' @param samtools path for samtools
#' @param maxins maximum fragment length (bowtie2 maxins argument)
#' @param threads number of threads to use
#' @param genome_type ct/ga/ct_ga. if 'ct' - Assume that the reads to be aligned underwent C->T conversion. If
#' paired-end reads are aligned, then assume read 1 underwent C->T
#' conversion while read 2 underwent G->A conversion. The --ct and --ga
#' options are mutually exclusive. if 'ga' - Assume that the reads to be aligned underwent G->A conversion. If
#' paired-end reads are aligned, then assume read 1 underwent G->A
#' conversion while read 2 underwent C->T conversion.
#' if 'ct_ga' reads would be aligned to both C->T and G->A, and the best match would
#' be chosen.
#'
#' @param tmp_dir Directory for storing temporary files
#' @param bissli2_params additional parameters to bissli2/bowtie2
#' @param log_file file for stderr / stdout of the bissli2 command
#' @param run_command run the bissli command
#' @param create_dirs create the bam file and log file directories
#'
#' @return if \code{run_command} is TRUE: bissli command exit status. if \code{run_command} is FALSE returns a string with the command
#' @export
#'
#' @examples
gpatterns.bissli2 <- function(r1_fastq,
                              out_bam,
                              genome_seq,
                              bissli2_idx,
                              r2_fastq=NULL,
                              bissli2_bin= system.file("bissli2", "bissli2-align.pl", package="gpatterns"),
                              bowtie2='bowtie2',
                              samtools = 'samtools',
                              maxins=1000,
                              threads=10,
                              genome_type='ct',
                              tmp_dir = NULL,
                              bissli2_params = '',
                              log_file = NULL, 
                              run_command = TRUE,
                              create_dirs = FALSE){
    if (create_dirs){
      system(qq('mkdir -p @{dirname(out_bam)}'))
      if (!is.null(log_file)){
        system(qq('mkdir -p @{dirname(log_file)}'))
      }
    }
    tmp_dir <- tmp_dir %||% tempdir()
    r1_fastq <- paste(r1_fastq, collapse=',')
    if (genome_type == 'ct_ga'){
        genome_type <- 'ct --ga'
    }

    log_str <- ''
    if (!is.null(log_file)){
      log_str <- qq(' 2> @{log_file}')
    }  

    if (is.null(r2_fastq)){
        command <- qq('@{bissli2_bin} @{bissli2_params} --tmp-dir @{tmp_dir} --bowtie2 @{bowtie2} --@{genome_type} -g @{genome_seq} -x @{bissli2_idx} -U @{r1_fastq} --threads @{threads} | @{samtools} view -b -S -h -o @{out_bam} - @{log_str}')
    } else {
        r2_fastq <- paste(r2_fastq, collapse=',')
        command <- qq('@{bissli2_bin} @{bissli2_params} --tmp-dir @{tmp_dir} --bowtie2 @{bowtie2} --maxins @{maxins} --@{genome_type} -g @{genome_seq} -x @{bissli2_idx} -1 @{r1_fastq} -2 @{r2_fastq} --threads @{threads} | @{samtools} view -b -S -h -o @{out_bam} -  @{log_str}')
    }

    if (!run_command){
      return(command)
    }
    
    system(command)
}

#' Create an bowtie2 index for a bisulfite converted genome
#'
#' @param reference a vector of FASTA file names holding the reference genome.
#' @param idx_base The basename of the index files that will be created.
#' @param bowtie2_build_bin The path of the bowtie2-build executable. [bowtie2-build]
#' @param bowtie2_options Any unknown options are passed transparently to bowtie2-build.
#' Note that bissli2-build can only handle FASTA files as input. Therefore the '-f' option is always passed to bowtie2-build
#' @param bissli2_build_bin binary of bissli2-build
#'
#' @return
#' @export
#'
#' @examples
gpatterns.bissli2_build <- function(reference,
                                    idx_base,
                                    bowtie2_options,
                                    bowtie2_build_bin='bowtie2-build',
                                    bissli2_build_bin= system.file("bissli2", "bissli2-build.pl", package="gpatterns")){
    reference <- paste(reference, collapse=',')
    command <- qq('@{bissli2_build_bin} --bowtie2-build @{bowtie2_build_bin} @{reference} @{idx_base} @{bowtie2_options}')
    system(command)
}

