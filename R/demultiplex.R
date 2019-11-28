
#' @export
gpatterns.demultiplex_fastqs <- function(config,
                                         workdir=NULL,
                                         fastq_dir = workdir,
                                         raw_reads_dir='{workdir}/{illumina_index}/raw', 
                                         split_dir = '{raw_reads_dir}/split',
                                         indexes_file='{split_dir}/indexes/indexes.tsv',
                                         raw_fastq_pattern = '.*_{read}_.*\\.fastq.gz', 
                                         R1_pattern = 'R1',
                                         R2_pattern = 'R2',                           
                                         indexes_tab=.gpatterns.all_indexes,
                                         idx1_pos = c(1,8),
                                         idx2_pos = c(1,8),
                                         umi1_pos=c(11, 19),
                                         umi2_pos=c(11, 19),
                                         read1_pos = c(20,80),
                                         read2_pos = c(20,80),
                                         hamming=1,
                                         paired_end=TRUE,
                                         commands_parallel=getOption('gpatterns.parallel'),
                                         use_sge=FALSE,
                                         stats_file=NULL,
                                         reads_per_file=4e6,
                                         illu_index=NULL,
                                         run_commands=TRUE,
                                         log_prefix=NULL,
                                         run_per_file=TRUE,
                                         step_file=NULL,
                                         ...){
    if (is.character(config)){
        config <- fread(config) %>% as.tibble()
    }


    if (!has_name(config, 'workdir')){
        if (is.null(workdir)){
            stop('Please supply a workdir (as an argument or as a column in config)')
        }
        config <- config %>% mutate(workdir = workdir)    
    } else {
        workdir <- config$workdir[1]
    }

    if (!is.null(illu_index)){
        config <- config %>% filter(illumina_index == illu_index)
    }

    if (!has_name(config, 'index1.seq')){
        config <- config %>% left_join(indexes_tab, by=c('index1' = 'index')) %>% rename(index1.seq = seq) %>% mutate(index1.seq = str_sub(start=idx1_pos[1], end=idx1_pos[2], string=index1.seq))
    }
    if (!has_name(config, 'index2.seq')){
        config <- config %>% left_join(indexes_tab, by=c('index2' = 'index')) %>% rename(index2.seq = seq) %>% mutate(index2.seq = str_sub(start=idx2_pos[1], end=idx2_pos[2], string=index2.seq))
    }

    indexes_config <- config %>% distinct(cell_id, row, plate_pos, column, index1, index2)

    fastq_files_pre <-  config %>% distinct(illumina_index, .keep_all=TRUE) %>% select(one_of('illumina_index', 'workdir', 'raw_reads_dir', 'split_dir', 'indexes_file')) %>% mutate(raw_reads_dir = glue(raw_reads_dir[1]), split_dir = glue(split_dir[1]))
    
    fastq_files <- get_raw_reads_files(fastq_files_pre, R1_pattern = R1_pattern, R2_pattern = R2_pattern, paired_end=paired_end)
    
    fastq_files <- fastq_files %>% mutate(raw_fastq_basename = gsub('R1_', '', gsub('\\.fastq\\.gz$', '', basename(raw_fastq_R1))))
    
    fastq_files <- fastq_files %>% mutate(indexes_file=glue(indexes_file[1]))

    config <- config %>% select(-one_of('raw_reads_dir', 'split_dir', 'indexes_file'))

    if (run_commands){
        res <- run_demultiplexing_commands(run_per_file, fastq_files, workdir, config, use_sge, log_prefix, idx1_pos, idx2_pos, umi1_pos, umi2_pos, read1_pos, read2_pos, hamming, reads_per_file, ...)
    }    
    
    config <- list_split_fastq_files(config, fastq_files, run_per_file)

    if (nrow(config) == 0){
        file.remove(step_file)
        logerror('No files were generated in demultiplexing. Please make sure that the indexes are correct.')
        stop('No files were generated in demultiplexing. Please make sure that the indexes are correct.')
    }

   if (!rlang::has_name(config, "empty")){
	config <- config %>% mutate(empty = TRUE)
   }

    missing_indexes <- indexes_config %>% left_join(config) %>% filter(is.na(r1_fastq), !empty)

    if (nrow(missing_indexes) > 0){
        missing_indexes <- missing_indexes %>% mutate(label = glue('{cell_id} ({plate_pos}: {index1}, {index2})')) %>% pull(label)
        red_message('The following indexes are missing:')
        walk(missing_indexes, ~ blue_message(.x))
        red_message('Please make sure that your indexes file is correct and run again with overwrite = TRUE')
    }

    # Prepare for mapping step    
    if (paired_end){
        config <- config %>% 
            mutate(r2_fastq = gsub(glue('good_{R1_pattern}'), glue('good_{R2_pattern}'), r1_fastq)) %>%
            mutate(index_exists = file.exists(r1_fastq) & file.exists(r2_fastq) & file.info(r1_fastq)$size != 0 & file.info(r2_fastq)$size != 0)
    } else {
        config <- config %>% mutate(index_exists = file.exists(r1_fastq) & file.info(r1_fastq)$size != 0)
    }

    config <- config %>% mutate(fastq_basename = gsub('\\.fastq.gz$', '', basename(r1_fastq)))

    if (!is.null(stats_file)){
        gpatterns.demultiplexing_stats(config)
    }

    return(config)  
}

list_split_fastq_files <- function(indexes_config, fastq_files, run_per_file){
    indexes_config <- indexes_config %>% select(-one_of('raw_reads_dir', 'split_dir', 'indexes_file'))

    if (run_per_file){
        # fastq_files <- fastq_files %>% distinct(illumina_index, raw_reads_dir, split_dir) %>% left_join(fastq_files %>% group_by(illumina_index) %>% nest(raw_fastq_R1, raw_fastq_R2, indexes_file, .key='raw_fastq'), by='illumina_index')
        fastq_files <- fastq_files %>% select(illumina_index, raw_reads_dir, split_dir, raw_fastq_R1, raw_fastq_R2, indexes_file) %>% group_by(illumina_index, raw_reads_dir, split_dir) %>% nest() %>% rename(raw_fastq = data)  %>% ungroup()
        indexes_config <- indexes_config %>% left_join(fastq_files, by='illumina_index')
    
        indexes_config <- indexes_config %>% plyr::adply(1, function(x) tibble(r1_fastq=list.files(x$split_dir, pattern=glue('.+_{x$lib}_good_R1(\\.\\d+)?\\.fastq.gz$'), full.names = TRUE) )) %>% as.tibble()    
    } else {
        # fastq_files <- fastq_files %>% distinct(illumina_index, raw_reads_dir, split_dir, indexes_file) %>% left_join(fastq_files %>% group_by(illumina_index) %>% nest(raw_fastq_R1, raw_fastq_R2, .key='raw_fastq'), by='illumina_index')
        fastq_files <- fastq_files %>% select(illumina_index, raw_reads_dir, split_dir, raw_fastq_R1, raw_fastq_R2) %>% group_by(illumina_index, raw_reads_dir, split_dir, indexes_file) %>% nest() %>% rename(raw_fastq = data)  %>% ungroup()
        

        indexes_config <- indexes_config %>% left_join(fastq_files, by='illumina_index')
    
        indexes_config <- indexes_config %>% plyr::adply(1, function(x) tibble(r1_fastq=list.files(x$split_dir, pattern=glue('{x$lib}_good_R1(\\.\\d+)?\\.fastq.gz$'), full.names = TRUE) )) %>% as.tibble()    
    }    
    return(indexes_config)
}

run_demultiplexing_commands <- function(run_per_file, fastq_files, workdir, config, use_sge, log_prefix, idx1_pos, idx2_pos, umi1_pos, umi2_pos, read1_pos, read2_pos, hamming, reads_per_file, ...){
    if (run_per_file){
        cmd_cfg <- fastq_files %>% distinct(raw_fastq_R1, .keep_all=TRUE) %>% mutate(illu_index=illumina_index, workdir=workdir, raw_fastq_basename = gsub('\\.fastq\\.gz', '', basename(raw_fastq_R1)))
        cmds <- glue('gpatterns:::do.call_tibble(gpatterns:::demultiplex_per_index, config, c(list(...), list(fastq_files=cmd_cfg[{1:nrow(cmd_cfg)}, ], config=config, run_per_file=run_per_file, idx1_pos=idx1_pos, idx2_pos=idx2_pos, umi1_pos=umi1_pos, umi2_pos=umi2_pos, read1_pos=read1_pos, read2_pos=read2_pos, hamming=hamming, reads_per_file=reads_per_file)))')   
    } else {
        cmd_cfg <- config %>% distinct(illumina_index, .keep_all=TRUE) %>% mutate(illu_index=illumina_index, workdir=workdir)    
        cmds <- glue('gpatterns:::do.call_tibble(gpatterns:::demultiplex_per_index, cmd_cfg[{1:nrow(cmd_cfg)}, ], c(list(...), list(fastq_files=fastq_files, config=config, run_per_file=run_per_file, idx1_pos=idx1_pos, idx2_pos=idx2_pos, umi1_pos=umi1_pos, umi2_pos=umi2_pos, read1_pos=read1_pos, read2_pos=read2_pos, hamming=hamming, reads_per_file=reads_per_file)))')       
    }

    loginfo('running %s commands', comify(length(cmds)))
    if (use_sge){       
        res <- gcluster.run2(command_list=cmds, ...)            
    } else {            
        res <- map(cmds, function(.x) eval(parse(text = .x)))  
    }        
    parse_commmands_res(res, cmds, cmd_cfg, use_sge=use_sge, log_prefix=log_prefix, ...)
    
    return(res)    
}

gpatterns.demultiplexing_stats <- function(config, R1_pattern = '.*_R1_.*\\.fastq.gz', R2_pattern = '.*_R2_.*\\.fastq.gz', paired_end=TRUE){
    fastq_files <- get_raw_reads_files(config)
    browser()
}

get_raw_reads_files <- function(config, raw_fastq_pattern = '.*_{read}_.*\\.fastq.gz', R1_pattern = 'R1', R2_pattern = 'R2', paired_end=TRUE){
    config <- config %>% select(one_of('illumina_index', 'raw_reads_dir', 'split_dir', 'indexes_file')) %>% distinct(illumina_index, raw_reads_dir, .keep_all = TRUE)
    
    fastq_files <- config %>% plyr::adply(1, function(x) 
            tibble(raw_fastq_R1 = list.files(x$raw_reads_dir, pattern=glue(raw_fastq_pattern, read=R1_pattern), full.names = TRUE)))
    if (paired_end){        
        fastq_files_R2 <- config %>% plyr::adply(1, function(x) 
            tibble(raw_fastq_R2 = list.files(x$raw_reads_dir, pattern=glue(raw_fastq_pattern, read=R2_pattern), full.names = TRUE))) %>% pull(raw_fastq_R2)

        fastq_files <- fastq_files %>% mutate(raw_fastq_R2 = gsub(R1_pattern, R2_pattern, raw_fastq_R1)) %>% filter(raw_fastq_R2 %in% fastq_files_R2)
    }
    return(as.tibble(fastq_files))
}


demultiplex_per_index <- function(fastq_files, config, illu_index=NULL, idx1_pos = c(1,8), idx2_pos = c(1,8), umi1_pos=c(11, 19), umi2_pos=c(11, 19), read1_pos = c(20,80), read2_pos = c(20,80), hamming=1, reads_per_file=4e6, run_per_file=TRUE){

    
    if (!is.null(illu_index)){
        fastq_files <- fastq_files %>% filter(illumina_index == illu_index)          
    }
    
    config <- config %>% filter(illumina_index %in% fastq_files$illumina_index) 

    if (length(unique(config$illumina_index)) > 1){
        stop('more than a single illumina index')
    } 

    split_dir <- paste0(fastq_files$split_dir[1], '/')
    indexes_file <- fastq_files$indexes_file[1]

    system(glue('mkdir -p {split_dir}'))
    system(glue('mkdir -p {dirname(indexes_file)}'))

    if (run_per_file){
        if (nrow(fastq_files) > 1){
           stop('more than a single file')
        }
        targets_pref <- paste0(split_dir, fastq_files$raw_fastq_basename[1], '_')
    } else {
        targets_pref <- split_dir
    }   
    
    conf2index(config, indexes_file, idx1_pos=idx1_pos, idx2_pos=idx2_pos, umi1_pos=umi1_pos, umi2_pos=umi2_pos, read1_pos =read1_pos, read2_pos = read2_pos, hamming=hamming[1], illu_index=illu_index, targets_pref=targets_pref)

    demultiplex_fastqs(R1_fastqs=fastq_files[['raw_fastq_R1']], R2_fastqs=fastq_files[['raw_fastq_R2']], indexes_file=indexes_file, reads_per_file=reads_per_file[1])
    
}

demultiplex_fastqs <- function(R1_fastqs, indexes_file, R2_fastqs=NULL,map_fastq_bin=.gpatterns.map_fastq_bin, reads_per_file=NULL){
    
    R1_str <- glue('-r1 <(gzip -d -c {paste(R1_fastqs, collapse=" ")})')
    R2_str <- ''
    if (!is.null(R2_fastqs)){
        R2_str <- glue('-r2 <(gzip -d -c {paste(R2_fastqs, collapse=" ")})')
    }
    reads_per_file_str <- ''
    if (!is.null(reads_per_file)){
        reads_per_file_str <- glue('-r {reads_per_file}')
    } 
    
    cmd <- glue('{map_fastq_bin} {R1_str} {R2_str} -i {indexes_file} {reads_per_file_str}')  
    loginfo('running %s', cmd)
        
    system(paste("/bin/bash -c", shQuote(cmd))) 
    
}

conf2index <- function(conf, out_fn, targets_pref='', idx1_pos = c(1,8), idx2_pos = c(1,8), umi1_pos=c(11, 19), umi2_pos=c(11, 19), read1_pos = c(20,80), read2_pos = c(20,80), hamming=1, illu_index=NULL, use_indexes=FALSE){
    
    if (!is.null(illu_index)){
        conf <- conf %>% filter(illumina_index == illu_index)
    }

    if (hamming > 0){        
        if (!is.null(idx1_pos)){
                idxs1 <- fill_hamming(conf, 'index1.seq', hamming=hamming, ext_columns='index1', rm_self=F)
                idxs1 <- rm_ambig_bcds(idxs1, 'index1.seq_hamming')  
                idxs1 <- idxs1 %>% full_join(conf %>% select(index1, index2, illumina_index, lib))                
        }
        if (!is.null(idx2_pos)){
                idxs2 <- fill_hamming(conf, 'index2.seq', hamming=hamming, ext_columns='index2', rm_self=F)
                idxs2 <- rm_ambig_bcds(idxs2, 'index2.seq_hamming')
                idxs2 <- idxs2 %>% full_join(conf %>% select(index1, index2, illumina_index, lib))   
        }
        
        if(is.null(idx1_pos)){
            idxs <- idxs2 %>% mutate(index1.seq_hamming = NA)
        } else if (is.null(idx2_pos)){
            idxs <- idxs1 %>% mutate(index2.seq_hamming = NA)           
        } else {            
            idxs <- idxs1 %>% full_join(idxs2)
        }
        
        conf <- idxs %>%                          
            distinct(index1.seq_hamming, index2.seq_hamming, lib, .keep_all=TRUE)
        
    } else {
        conf <- conf %>% mutate(index1.seq_hamming = index1.seq, index2.seq_hamming = index2.seq)
    print(conf)
    }
    
    conf <- add_file_path(conf, use_indexes, targets_pref=targets_pref)       
    
    
    if(is.null(idx1_pos)){
        conf <- conf %>% select(-index1.seq)
        colnames(conf) <- c('target',                
                sprintf("umi|r1:%s-%s", umi1_pos[1], umi1_pos[2]),
                sprintf("read:R1|r1:%s-%s", read1_pos[1], read1_pos[2]),
                sprintf("idx:index2|r2:%s-%s", idx2_pos[1], idx2_pos[2]),
                sprintf("umi|r2:%s-%s", umi2_pos[1], umi2_pos[2]),
                sprintf("read:R2|r2:%s-%s", read2_pos[1], read2_pos[2])
                )
    } else if (is.null(idx2_pos)){
        conf <- conf %>% select(-index2.seq)
        colnames(conf) <- c('target',
                sprintf("idx:index1|r1:%s-%s", idx1_pos[1], idx1_pos[2]),
                sprintf("umi|r1:%s-%s", umi1_pos[1], umi1_pos[2]),
                sprintf("read:R1|r1:%s-%s", read1_pos[1], read1_pos[2]),                
                sprintf("umi|r2:%s-%s", umi2_pos[1], umi2_pos[2]),
                sprintf("read:R2|r2:%s-%s", read2_pos[1], read2_pos[2])
                )        
    } else {            
        colnames(conf) <- c('target',
                sprintf("idx:index1|r1:%s-%s", idx1_pos[1], idx1_pos[2]),
                sprintf("umi|r1:%s-%s", umi1_pos[1], umi1_pos[2]),
                sprintf("read:R1|r1:%s-%s", read1_pos[1], read1_pos[2]),
                sprintf("idx:index2|r2:%s-%s", idx2_pos[1], idx2_pos[2]),
                sprintf("umi|r2:%s-%s", umi2_pos[1], umi2_pos[2]),
                sprintf("read:R2|r2:%s-%s", read2_pos[1], read2_pos[2])
                )
    }    
    
    loginfo('writing to %s', out_fn)  
    write_tsv(conf, out_fn)        
}


.split2char <- function(df, column, ext_columns=NULL){
    barcode_len <- str_length(df[1, column])
    if (!is.null(ext_columns)){
        df <- df[, c(column, ext_columns)]
    } else {
        df <- df %>% select_(column)
    }
    df %>% distinct_(column, .keep_all=TRUE) %>% mutate_(.dots = setNames(list(lazyeval::interp(~ as.character(a), a = as.name(column))), column)) %>% separate_(column, as.character(1:barcode_len), 1:(barcode_len-1), remove=F) %>% gather_('pos', 'char', paste(1:barcode_len)) %>% arrange_(column)
}

.unite_by_char <- function(df, column, new_column=NULL){
    barcode_len <- as.numeric(max(df$pos))
    if (is.null(new_column)){
        new_column <- paste0(column, '_hamming')
    }
    df %>% spread('pos', 'char') %>% unite_(new_column, paste(1:barcode_len), sep='')
}

.expand_hamming <- function(df, column, alphabet, barcode_len, new_column=NULL, ext_columns=NULL){
    df <- .split2char(df, column, ext_columns)
    plyr::adply(alphabet, 1, function(letter) plyr::adply(1:barcode_len, 1, function(p) df %>% mutate(char = ifelse(pos == p, letter, char)) %>% .unite_by_char(column, new_column=new_column)) %>% rename(pos = X1) %>% mutate(letter=letter)) %>% select(-X1) %>% select(-pos, -letter)
}


fill_hamming <- function(df, column, hamming=1, alphabet=c('A', 'C', 'T', 'G'), ext_columns=NULL, rm_self=TRUE){
    barcode_len <- str_length(df[1, column])
    
    loginfo('doing hamming %s', 1)
    hdf <- .expand_hamming(df, column, alphabet, barcode_len, ext_columns=ext_columns, new_column='h1')
    if (hamming > 1){        
        for (i in 2:hamming){
            loginfo('doing hamming %s', i)
            add_cols <- c(column, ext_columns)
            if (i > 2){
                add_cols <- c(add_cols, paste0('h', 1:(i-2)))
            }
            hdf <- .expand_hamming(hdf, paste0('h', i-1), alphabet, barcode_len, ext_columns=add_cols, new_column=paste0('h', i))   
        }
    }
    hdf <- hdf %>% gather_('hamming', paste0(column, '_hamming'), paste0('h', 1:hamming)) %>% select(-hamming) %>% distinct_(column, paste0(column, '_hamming'), .keep_all=TRUE)
    if (rm_self){
        hdf <- hdf %>% filter_(paste0(column, ' != ', column, '_hamming'))
    }
    return(hdf)
}


rm_ambig_bcds <- function(df, hamming_column){
    df <- df %>% group_by_(hamming_column) %>% mutate(n = n()) %>% ungroup
    loginfo("removed %s ambiguous barcodes", sum(df$n > 1))
    df %>% filter(n == 1) %>% select(-n)
}

add_file_path <- function(df, use_indexes=FALSE, targets_pref=''){    
    if (use_indexes){
            df <- df %>%
                filter(!is.na(lib)) %>%
                mutate(suf = ifelse(lib == 'empty', 'bad', 'good')) %>% 
                mutate(target = paste0(index1, '_', index2)) %>%  
                select(illumina_index, target, index1.seq=index1.seq_hamming, index2.seq=index2.seq_hamming, suf) %>%
                mutate(target = sprintf("%s%s_%s", targets_pref, target, suf)) %>%
                select(-illumina_index, -suf) %>% 
                distinct(target, index1.seq, index2.seq, .keep_all=TRUE) %>%
                mutate(umi1=NA, umi2=NA, read1=NA, read2=NA) %>%
                select(target, index1.seq, umi1, read1, index2.seq, umi2, read2)               
    } else {
        df <- df %>%
            filter(!is.na(lib)) %>%
            mutate(suf = ifelse(lib == 'empty', 'bad', 'good')) %>% 
            rename(target = lib) %>%            
            select(illumina_index, target, index1.seq=index1.seq_hamming, index2.seq=index2.seq_hamming, suf) %>%
            mutate(target = sprintf("%s%s_%s", targets_pref, target, suf)) %>%
            select(-illumina_index, -suf) %>% 
            distinct(target, index1.seq, index2.seq, .keep_all=TRUE) %>%
            mutate(umi1=NA, umi2=NA, read1=NA, read2=NA) %>%
            select(target, index1.seq, umi1, read1, index2.seq, umi2, read2)        
    }
    return(df)

}
