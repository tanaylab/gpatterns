
gpatterns.map <- function(config, workdir=NULL, out_bam='{workdir}/{illumina_index}/bam/{fastq_basename}.bam', log_file = '{workdir}/{illumina_index}/bam/log/{fastq_basename}.log', commands_parallel=getOption('gpatterns.parallel'), use_sge=FALSE, run_commands=TRUE, paired_end=TRUE, verify=TRUE, log_prefix=NULL, overwrite=TRUE, step_file=NULL, ...){
    # indexes_config <- config %>% distinct(cell_id, row, plate_pos, column, index1, index2)
    
    if (paired_end){
        config <- config %>% mutate(index_exists = file.exists(r1_fastq) & file.exists(r2_fastq) & file.info(r1_fastq)$size != 0 & file.info(r2_fastq)$size != 0)
    } else {
        config <- config %>% mutate(index_exists = file.exists(r1_fastq) & file.info(r1_fastq)$size != 0)
    }

    # fastq_config <- config

    config <- config %>% filter(index_exists)

    if (!has_name(config, 'workdir')){
        if (is.null(workdir)){
            stop('Please supply a workdir (as an argument or as a field in config)')
        }
        config <- config %>% mutate(workdir = workdir)    
    }   

    if (!is.null(out_bam)){        
        config <- config %>% select(-one_of('out_bam')) %>% mutate(out_bam = glue::glue_data(., out_bam))
    }

    if (!is.null(log_file)){        
        config <- config %>% select(-one_of('log_file')) %>% mutate(log_file = glue(log_file))
    }

    cmds <- glue('gpatterns:::do.call_tibble(gpatterns::gpatterns.bissli2, config[{1:nrow(config)}, ], c(list(...), create_dirs=TRUE))')

    if (has_name(config, 'overwrite')){
        overwrite <- config$overwrite[1]
    }

    if (!overwrite){        
        config <- config %>% mutate(run_map = !file.exists(out_bam))
        cmds <- cmds[which(!file.exists(config$out_bam))]
        loginfo(sprintf('skipping %d files that already exist. set overwrite to TRUE to override', nrow(config) - length(cmds)))
    }    

    if (run_commands && length(cmds) > 0){
        loginfo('running %s commands', comify(length(cmds)))
        if (use_sge){       
            res <- gcluster.run2(command_list=cmds, ...)
            
        } else {
             res <- map(cmds, function(.x) eval(parse(text = .x)))  
        }
        parse_commmands_res(res, cmds, config %>% filter(run_map), use_sge=use_sge, log_prefix=log_prefix, ...)
    }        
  

    config <- config %>% rename(bam_file = out_bam) %>% mutate(bam_exists = file.exists(bam_file)) 
    failed_bams <- config$bam_file[!config$bam_exists]
    if (length(failed_bams) > 0){
        red_message('The following bam files were not created:')
        walk(failed_bams, ~ blue_message(.x))
        logerror('some bams were not created: %s', paste(failed_bams, collapse = ', '))
    }
      
    
    return(config)
}

gpatterns.bam2tracks <- function(config, track_prefix='', track='{track_prefix}{experiment}.{lib}', description='experiment: {experiment}. library: {lib}', steps='all', commands_parallel=getOption('gpatterns.parallel'), use_sge=FALSE, run_commands=TRUE, log_prefix=NULL, ...){
    config <- config %>% mutate(bam_exists = file.exists(bam_file))

    orig_config <- config
    config <- config %>% distinct(experiment, lib, .keep_all=TRUE) %>% select(-bam_file)
    

    # nested_cfg <- orig_config %>% group_by(experiment, lib) %>% nest(bam_file, .key='bams')
    nested_cfg <- orig_config %>% select(experiment, lib, bam_file) %>% group_by(experiment, lib) %>% nest() %>% rename(bams = data) %>% mutate(bams = map(bams, ~ .x$bam_file))
    config <- config %>% left_join(nested_cfg, by=c('experiment', 'lib'))
    for( .x in c('r1_fastq', 'r2_fastq', 'fastq_basename')) { 
        if (has_name(config, .x)) { 
            nested_cfg <- orig_config %>% select(experiment, lib, !!.x) %>% group_by(experiment, lib) %>% nest() %>% rename(!!.x := data)  %>% ungroup()        
            config <- config %>% select(-one_of(.x)) %>% left_join(nested_cfg, by=c('experiment', 'lib'))
        }
    }

    if (!has_name(config, 'track_prefix')){
        config <- config %>% mutate(track_prefix = track_prefix)
    }

    if (!is.null(track)){        
        config <- config %>% select(-one_of('track')) %>% mutate(track = glue::glue_data(., track))
    }

    if (!is.null(description)){        
        config <- config %>% select(-one_of('description')) %>% mutate(description = glue::glue_data(., description))
    }

    run_cfg <- config %>% select(-one_of('workdir'))
    
    cmds <- glue('gpatterns:::do.call_tibble(gpatterns::gpatterns.import_from_bam, run_cfg[{1:nrow(run_cfg)}, ], additional_params=list(...), additional_funcs=c(gpatterns::gpatterns.import_from_tidy_cpgs))')        
    
    if (run_commands){
        loginfo('running %s commands', comify(length(cmds)))
        if (use_sge){       
            res <- gcluster.run2(command_list=cmds, ...)
        } else {
             res <- map(cmds, function(.x) eval(parse(text = .x)))  
        }
        parse_commmands_res(res, cmds, config, use_sge=use_sge, log_prefix=log_prefix, ...)
    }        
    
    
    if (track_prefix == ''){
        config <- config %>% select(-track_prefix)
    }
    return(config)

}