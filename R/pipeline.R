#' Run pipeline
#' @export
gpatterns.pipeline <- function(config_file, log_file=NULL, config_dir=NULL, run_commands=TRUE, run_dir = getwd(), overwrite=FALSE){
	logging::basicConfig()

    if (!is.null(log_file)){
        logging::addHandler(logging::writeToFile, file=log_file)
        logging::removeHandler('basic.stdout')
        on.exit(logging::removeHandler('logging::writeToFile'))
    }

	yaml <- yaml::yaml.load_file(config_file)

	exp_cfg <- yaml2cfg(yaml)	
	by_vars <- intersect(names(exp_cfg %>% unnest(config, names_repair = tidyr_legacy)), names(exp_cfg %>% select(-config, -pipeline_steps)))
    by_vars <- by_vars[!map_lgl(by_vars, ~ is.list(exp_cfg[[.x]]))]
	
	all_cfg <- exp_cfg %>% 
		unnest(config, names_repair = tidyr_legacy) %>%
		left_join(exp_cfg %>% select(-config, -pipeline_steps), by=by_vars)

	exp_steps <- all_cfg %>% distinct(experiment, pipeline_steps) %>% unnest(pipeline_steps)

	if (has_name(all_cfg, 'illumina.index')){
		all_cfg <- all_cfg %>% rename(illumina_index = illumina.index)
	}

	all_steps <- names(yaml$steps)

	config <- all_cfg

	
	for (step in all_steps){
		config$overwrite <- overwrite			
		
		step_file <- glue('{run_dir}/finished_{step}')			
		if (!file.exists(step_file) || overwrite){
			red_message('{step}')	
            
			config <- invoke_step(config, yaml$steps, step, exp_steps, step_file=step_file, run_commands=run_commands)
			if (run_commands){
				file.create(step_file)		
			}            
		} else {
			if (!run_commands){				
				red_message("skipping {step} (to run set run_commands to TRUE)")
			}
			if (!overwrite){				
				red_message("skipping {step} (to override set overwrite to TRUE)")
			}
			
			config <- invoke_step(config, yaml$steps, step, exp_steps, step_file=step_file, run_commands=FALSE)
		}       
        
		
		if (!is.null(config_dir)){			
            readr::write_rds(config, glue("{config_dir}/{step}_config.rds"))			
		}
	}	
		
	
}

invoke_step <- function(config, steps_cfg, step, exp_steps=NULL, step_file=NULL, run_commands=TRUE){
	if (!is.null(exp_steps)){
		exp_steps <- exp_steps %>% filter(pipeline_steps == step)
		step_config <- config %>% filter(experiment %in% exp_steps$experiment)	
	}

    if (nrow(step_config) > 1){
    	if (run_commands){
    		loginfo(steps_cfg[[step]]$description)
    		run_step(step_config, steps_cfg, step, step_file, run_commands) 
    	}     	
    }    
    
    # dry run to get config tibble
    return(run_step(config, steps_cfg, step, step_file, run_commands=FALSE))
    
}

run_step <- function(config, steps_cfg, step, step_file, run_commands){	
	func <- eval(parse(text=steps_cfg[[step]]$func))	
	args_vec <-  steps_cfg[[step]][['params']]

    if (run_commands){
    	loginfo('running `%s` with the following paramters:', steps_cfg[[step]]$func)
    	walk2(names(args_vec), args_vec, function(x, y) loginfo("%s: %s", x, y))	
    }    
    
    res <- do.call(func, c(list(config = config, step_file=step_file, run_commands=run_commands), args_vec))	
    if (run_commands){    	
    	red_message('finished {step}')
    }    
    return(res)
}

yaml2cfg <- function(conf){
	
	conf[['experiments']] <- map(conf[['experiments']], ~ apply_genome_conf(.x, conf))
	
	list_columns <- map(conf[['experiments']], names) %>% simplify() %>% unique() %>% keep(function(field) any(map_int(conf[['experiments']], ~ length(.x[[field]]) ) > 1))
	
	exp_cfg <- map2_dfr(conf[['experiments']],
					  names(conf[['experiments']]), 
					  ~ .x %>% 
					  	keep(~ !is.null(.x)) %>% 					  	
					  	modify_at(list_columns, ~ list(.x)) %>%
					  	as_tibble() %>%
					  	mutate(experiment=.y))

	# get indexes
	exp_cfg <- exp_cfg %>% mutate(config = map(config_file, ~ { 
			df <- fread(.x) %>% as_tibble() %>% mutate_if(is.numeric, as.character)
			if (!has_name(df, 'lib') && has_name(df, 'sample')){
				df <- df %>% unite('lib', contains('sample'), remove=FALSE)
			}
			return(df)
		} ))
	return(exp_cfg)

}


apply_genome_conf <- function(conf, main_conf){
    if ('genome' %in% names(conf) && 'genomes' %in% names(main_conf)){
        conf <- plyr::defaults(conf, main_conf$genomes[[conf$genome]])
    }
    return(conf)
}

parse_commmands_res <- function(res, cmds, config, use_sge=FALSE, log_prefix=NULL, ...){	
    if (use_sge){               
        exit_status <- map_chr(res, 'exit.status')
        if (any(exit_status != 'success')){
            err_id <- which(exit_status != 'success')
            logwarn('command %s failed', cmds[err_id])
        }
        
        if (!is.null(log_prefix)){        	
            config <- config %>% mutate(log_err = glue(log_prefix, '.err'), log_out =  glue(log_prefix, '.out'))
            walk(glue('mkdir -p {unique(dirname(config$log_err))}'), system)
            walk2(1:nrow(config), res, ~ readr::write_lines(.y$stdout, config %>% slice(.x) %>% .$log_out))
            walk2(1:nrow(config), res, ~ readr::write_lines(.y$stderr, config %>% slice(.x) %>% .$log_err))
        }
    } else {        
        exit_status <- map_int(res, ~ .x)
        if (!any(exit_status != 0)){
            err_id <- which(exit_status != 0)
            logwarn('command %s failed', cmds[err_id])   
        }
    }           
}

