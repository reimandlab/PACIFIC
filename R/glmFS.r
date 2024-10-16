#' Current Total Iterations
#'
#' Get the total number of iterations currently stored in a directory
#' @param dir WIP
#' @return WIP
#' @examples 
#' WIP
#' @export
current_total_iters <- function(dir){
    r <- list.files(dir, pattern = 'step1-results')
    r <- regmatches(r, regexec('-N-iters-\\s*(.*?)\\s*-timestamp-', r))
    N <- sum(as.numeric(sapply(r, '[[', 2)))
    return(N)
}

#' glmFS Step 1
#'
#' Run the first step of glmFS
#' @param data WIP
#' @param response WIP
#' @param event WIP
#' @param single_features WIP
#' @param interaction_features WIP
#' @param features_to_discretize WIP
#' @param discretization_method WIP
#' @param features_with_flexible_direction WIP
#' @param features_to_skip_sparsity_prefiltering WIP
#' @param features_to_skip_univariate_association_prefiltering WIP
#' @param sparsity_criterion WIP
#' @param univariate_p_cutoff WIP
#' @param subsampling_ratio WIP
#' @param num_iterations WIP
#' @param output_dir WIP
#' @param verbose WIP
#' @return WIP
#' @examples 
#' WIP
#' @export
glmFS_step1 <- function(data,
                        response,
                        event = NA,
                        single_features = NA,
                        interaction_features = NA,
                        features_to_discretize = NA,
                        discretization_method = 'median',
                        features_with_flexible_direction = NA,
                        features_to_skip_sparsity_prefiltering = NA,
                        features_to_skip_univariate_association_prefiltering = NA,
                        sparsity_criterion = '5_percent',
                        univariate_p_cutoff = 0.1,
                        subsampling_ratio = 0.8,
                        num_iterations = 10,
                        output_dir,
                        verbose = FALSE){

    S.TM <- Sys.time()
    
    if(verbose){ cat('----------------------------------------------------\ncall glmFS_step1:\n'); flush.console() }
    if(verbose){ cat('preprocessing ... '); flush.console() }
    
    # validate and process arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # create output_dir if it does not exist
    if(!dir.exists(output_dir)) dir.create(output_dir, recursive = T, showWarnings = F)
    # input arguments are saved in a specific file in output_dir
    pth <- paste0(output_dir, '/step1-args.RData')
    # if this is the first iteration, save the input arguments in the output_dir
    if(!file.exists(pth)){
        # all args saved except for "num_iteratoins" as it may differ between runs
        save(list = setdiff(ls(), c('num_iterations','verbose','S.TM')), file = pth)
    } else { # otherwise, check the consistency of input arguments with previous call(s)
        wait_for_file_to_stabilize(pth)
        load(pth, previous_args <- new.env())
        for(v in ls(previous_args)){
            if(!all(TRUE == suppressWarnings(all.equal(get(v), previous_args[[v]])))){
                stop(paste0('"', v, '": argument inconsistent with previous iterations'))
            }
        }    
    }
    # validate arg data
    if(!is.data.frame(data)) stop('"data" should be a data.frame or an instance of a class extended from data.frame')
    data <- as.data.frame(data)
    # validate arg response
    stopifnot(is.character(response), length(response)==1, response %in% colnames(data))
    # validate arg event (needed only for survival outcomes)
    if(!identical(NA, event)) stopifnot(is.character(event), length(event)==1, event %in% colnames(data))
    # available data columns
    avDtCls <- setdiff(colnames(data), c(response, event, "response", "event")) # "response" and "event" are **reserved** column names and cannot be used for features
    # validate and prepare vectors of single features
    validPrep_vector_of_single_features <- function(x){
        # empty chr vec if NA
        if(identical(NA, x)) x <- character(0)
        # assert it is subset of available data columns without duplicates
        stopifnot(x %in% avDtCls, !duplicated(x))
        # return itself
        return(x)
    }
    single_features <- validPrep_vector_of_single_features(single_features)
    features_to_discretize <- validPrep_vector_of_single_features(features_to_discretize)    
    features_to_skip_sparsity_prefiltering <- validPrep_vector_of_single_features(features_to_skip_sparsity_prefiltering)
    features_to_skip_univariate_association_prefiltering <- validPrep_vector_of_single_features(features_to_skip_univariate_association_prefiltering)
    features_with_flexible_direction <- validPrep_vector_of_single_features(features_with_flexible_direction)
    # validate and prepare interaction_features
    
    if(identical(NA, interaction_features)){
        interaction_features <- character(0)
    } else {
        if(class(interaction_features) == 'list'){
            stopifnot(2 == length(interaction_features))
            interaction_features[[1]] <- validPrep_vector_of_single_features(interaction_features[[1]])
            interaction_features[[2]] <- validPrep_vector_of_single_features(interaction_features[[2]])
            if(any(duplicated(unlist(interaction_features)))) stop('found duplicated element in "interaction_features"')
            if(!all(unlist(interaction_features) %in% avDtCls)) stop('"interaction_features" must be subset of available data columns')
            interaction_features <- apply(expand.grid(interaction_features), 1, paste, collapse='*')
        } else {
            stopifnot(is.character(interaction_features))
            sp <- strsplit(interaction_features, '\\*')
            if(any(2 != sapply(sp, length))) stop('"interaction_features" must be chr. vector of pairs with "*" or a list of two chr. vectors')
            if(any(2 != sapply(sapply(sp, unique, simplify = F), length))) stop('found self-interaction in "interaction_features"')
            if(any(duplicated(t(sapply(sp, sort))))) stop('found duplicated pair in "interaction_features"')
            if(!all(unlist(sp) %in% avDtCls)) stop('"interaction_features" must be subset of available data columns')
        }
    }
    # throw error if NO input features
    if(0 == (length(single_features)+length(interaction_features))) stop('NO input feature')
    # s"ingle_features" cannot overlap with any element of "interaction_features"
    if(any(single_features %in% unlist(strsplit(interaction_features, '\\*')))){
        stop('"single_features" cannot overlap with any element of "interaction_features"')
    }
    # features to skip any prefiltering cannot overlap with any element of "interaction_features"
    if(any(c(features_to_skip_univariate_association_prefiltering, 
             features_to_skip_sparsity_prefiltering) %in% unlist(strsplit(interaction_features, '\\*')))){
        stop('features to skip any prefiltering cannot overlap with any element of "interaction_features"')
    }
    # validate discretization_method
    if(is.function(discretization_method)){
        discretization_function <- discretization_method
        # test the discretization_function
        n <- 100*(runif(1000)-0.5)
        f <- discretization_method(n)
        if(!is.factor(f)) stop('invalid discretization_method')
        if(length(f) != length(n)) stop('invalid discretization_method')
    } else if(identical('median', discretization_method)) {
        # make the discretization_function
        discretization_function <- function(x){
            stopifnot(is.numeric(x))
            lvs <- c('lower_than_or_equal_to_median', 'higher_than_median')
            return(factor(ifelse(x <= median(x), lvs[1], lvs[2]), levels = lvs))
        }        
    } else {
        stop('invalid discretization_method')
    }    
    # validate "sparsity_criterion"
    if(!all(c(1 == length(sparsity_criterion), grepl('^\\d+_(percent|absolute)$', sparsity_criterion)))){
        stop('invalid sparsity_criterion')
    }
    # validate the remaining single params
    stopifnot(1 == length(univariate_p_cutoff), is.numeric(univariate_p_cutoff), univariate_p_cutoff > 0, univariate_p_cutoff < 1)
    stopifnot(1 == length(subsampling_ratio), is.numeric(subsampling_ratio), subsampling_ratio > 0, subsampling_ratio < 1)
    stopifnot(1 == length(num_iterations), is.numeric(num_iterations), num_iterations > 0)
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    # collect the "input features", validate and prepare data columns >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    features <- c(single_features, interaction_features)
    # *** include individual elements of "interaction_features" in the "features"
    features <- c(features, unique(unlist(strsplit(interaction_features, '\\*'))))
    # sanity check: no overlap between "single_features" and any element of "interaction_features"
    stopifnot(all(!duplicated(features)))
    # set "response" and "event" columns in the data
    data$response <- data[[response]]; data$event <- data[[event]]
    # validate "response" (and "event" in case of survival outcomes) and then set the "family"
    if(!identical(NA, event)){
        # if "event" argument is specified, the model response is survival outcome, set family to "cox"
        family <- 'cox'
        # require response to be "time values"
        if(!all(c(is.numeric(data$response), data$response > 0))) stop('For survival outcomes, "response" column must be numeric with positive values')
        # validate event values
        if(!all(c(is.numeric(data$event), data$event %in% c(0, 1)))) stop('For survival outcomes, "event" column must be numeric with only 0 & 1 values')
        if(sum(data$event == 1) < 5) stop('For survival outcomes, there must be at least 5 events (5 rows with event=1) in "data"')
    } else {
        errMsg <- 'Invalid response values. Expected either "numeric" (for continuous response) or "2-level factor" (for binary response)'
        # is event = NA, the response is either continuous or binary
        if(class(data$response) == 'factor'){
            if(2 != length(levels(data$response))) stop(errMsg)
            # for binary outcome, set family to "binomial"
            family <- 'binomial'
        } else if(class(data$response) == 'numeric'){
            # for continuous outcome, set family to "gaussian"
            family <- 'gaussian'
        } else {
            stop(errMsg)
        }
    }
    # loop over all single elements of features:
    for(f in unique(unlist(strsplit(features, '\\*')))){
        # every element of features should map to either a numeric or factor in data
        # if the column in data is non-numeric & non-factor, convert it to factor ***
        if(!is.numeric(data[[f]]) & !is.factor(data[[f]])){
            data[[f]] <- factor(data[[f]])
            # it is plausible to assume that any feature that is factorized in this step has
            # flexible direction. If not, it should have been set as factor in input data.
            # include this feature in "features_with_flexible_direction"
            features_with_flexible_direction <- union(f, features_with_flexible_direction)
        }
        # any feature in "features_with_flexible_direction" must be either factor or must be in "features_to_discretize"
        if(f %in% features_with_flexible_direction){
            if(!is.factor(data[[f]])){
                if(! f %in% features_to_discretize){
                    stop(paste0('"', f, '" is in "features_with_flexible_direction" so it must be ',
                                'either a "factor" or in "features_to_discretize"'))
                }
            }
        }
        # any factor must have at least two levels
        if(is.factor(data[[f]])){
            # stopifnot(length(levels(data[[f]])) >= 2)
            if(length(levels(data[[f]])) < 2){
                stop(paste0('"', f, '" as a factor must have at least two levels.'))
            }
        }
        # there must NOT be any missing data
        # stopifnot(all(!is.na(data[[f]])))
        if(any(is.na(data[[f]]))){
            stop(paste0('"', f, '" has missing (NA) value(s).'))
        }
        
    }; rm(f)
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    # >>> calculate "min_group_size" for sparsity filtering
    x <- unlist(strsplit(sparsity_criterion, '_'))
    threshold <- as.numeric(x[1])
    method <- x[2]
    if(method == 'percent') min_group_size <- subsampling_ratio * nrow(data) * threshold / 100
    if(method == 'absolute') min_group_size <- threshold

    # >>> discretize "data" for "features_to_discretize" -> "data_disc"
    data_disc <- discretize(data = data, 
                            features = features_to_discretize, 
                            discretization_function = discretization_function)

    # >>> make the table of "features" based on "data_disc" and "input features"
    features <- make_table_of_features(data=data_disc, 
                                       features_vector = features, 
                                       single_features = single_features, 
                                       flexible_features = features_with_flexible_direction)

    # >>> apply sparsity filter to "data_disc" and "features" with "min_group_size" and "features_to_skip_sparsity_prefiltering".
    # the feature IDs detected as sparse at this step (using the whole set of samples) will necessarily 
    # be detected as sparse in the iterations (using subsets of samples). So it's reasonalbe to 
    # remove them upstream of the iterations.
    features <- filter_by_sparsity(data = data_disc, 
                                   features = features, 
                                   skip_ids = features_to_skip_sparsity_prefiltering, 
                                   min_group_size = min_group_size)
    if(0 == nrow(features)) stop('All input features are detected as sparse. Cannot continue to feature selection.')

    # >>> get data table of variables based on "data_disc" and "features" -> "data_vars"
    data_vars <- make_data_for_variables(data = data_disc, features = features)

    # >>> create the "records", save in "output_dir", and check the consistency with previous iterations
    # "records" = list("family", "data", "data_vars", "features")
    # to save in "records", remove "group" column from "features" and unique its rows
    temp_features <- unique(subset(features, select = -group))
    rownames(temp_features) <- NULL
    records <- list(family = family,
                    data = data, 
                    data_vars = data_vars, 
                    features = temp_features)
    rm(temp_features)
    pth <- paste0(output_dir, '/step1-records.rds')
    # if first iteration: save "records" in output_dir, else: "records" must be identical to that of previous iteratoins.
    if(!file.exists(pth)){ saveRDS(records, pth) } else {
        wait_for_file_to_stabilize(pth)
        if(!identical(records, readRDS(pth))){
            stop('the processed data and features are inconsistent with previous iterations')
        }
    }

    F.TM <- Sys.time()
    if(verbose){ cat(format(F.TM - S.TM), '\n'); flush.console() }
    
    # >>> iterative process: 
    # --- 1: subsample
    # --- 2: discretize
    # --- 3: filter by sparsity
    # --- 4: make data table of variables
    # --- 5: filter by univariate associations
    # --- 6: elastic net feature selection
    res <- lapply(1:num_iterations, function(iter){
        S.TM <- Sys.time()
        if(verbose){ cat(paste0('iteration ', iter, '/', num_iterations, ' : ')); flush.console() }

        # >>> (1) subsample "data" -> "data_samp"
        sampleRows <- function(df, ratio) df[sample(nrow(df), floor(ratio*nrow(df))), ]
        if(family == 'cox'){
            # for cox family, subsample by maintaining the ratio of the event in data
            data_samp <- rbind(sampleRows(subset(data, event==0), subsampling_ratio),
                               sampleRows(subset(data, event==1), subsampling_ratio))
        } else if(family == 'binomial'){
            # for binomial family, subsample by maintaining the ratio of the binary outcome in data
            data_samp <- rbind(sampleRows(subset(data, response==levels(data$response)[1]), subsampling_ratio),
                               sampleRows(subset(data, response==levels(data$response)[2]), subsampling_ratio))
        } else if(family == 'gaussian'){
            # for gaussian family, subsample freely
            data_samp <- sampleRows(data, subsampling_ratio)
        } else {
            stop('invalid family')
        }

        # >>> (2) discretize "data_samp" for "features_to_discretize" -> "data_samp_disc"
        data_samp_disc <- discretize(data = data_samp, 
                                     features = features_to_discretize, 
                                     discretization_function = discretization_function)
    
        # >>> (3) apply sparsity filter to "data_samp_disc" and "features" with "min_group_size" and "features_to_skip_sparsity_prefiltering".
        features_flt1 <- filter_by_sparsity(data = data_samp_disc, 
                                            features = features, 
                                            skip_ids = features_to_skip_sparsity_prefiltering, 
                                            min_group_size = min_group_size)
        if(0 == nrow(features_flt1)){
            F.TM <- Sys.time()
            if(verbose){ cat('No feature passed the sparsity filter.', format(F.TM - S.TM), '\n'); flush.console() }
            return(NULL)
        }
        
        # >>> (4) get data table of variables based on "data_samp_disc" and "features_flt1" -> "data_samp_vars"
        data_samp_vars <- make_data_for_variables(data = data_samp_disc, features = features_flt1)

        # >>> (5) apply univariate association filter to "data_samp_vars" for "features_flt1" using "univariate_p_cutoff" taking care of "features_to_skip_univariate_association_prefiltering"
        features_flt2 <- filter_by_univariate_association(family = family,
                                                          data = data_samp_vars, 
                                                          features = features_flt1, 
                                                          skip_ids = features_to_skip_univariate_association_prefiltering, 
                                                          P_cutoff = univariate_p_cutoff)
        # check integrity of "features_flt2" >>>
        # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        if(0 == nrow(features_flt2)){
            F.TM <- Sys.time()
            if(verbose){ cat('No feature passed the filter of univariate association.', format(F.TM - S.TM), '\n'); flush.console() }
            return(NULL)
        }
        if(1 == nrow(features_flt2)){
            F.TM <- Sys.time()
            if(verbose){ cat('Only a single feature passed prefiltering, * selected *.', format(F.TM - S.TM), '\n'); flush.console() }
            return(features_flt2$variable)
        }
        
        # >>> (6) apply Elastic net feature selection to "data_samp_vars" and "features_flt2"
        selected_variables <- elastic_net_variable_selection(family = family, 
                                                             data = data_samp_vars, 
                                                             variables = features_flt2$variable)
        if(0 == length(selected_variables)){
            F.TM <- Sys.time()
            if(verbose){ cat('No feature selected by Elastic net.', format(F.TM - S.TM), '\n'); flush.console() }
            return(NULL)
        }
        F.TM <- Sys.time()
        if(verbose){ cat(length(selected_variables), 'features selected by Elastic net.', format(F.TM - S.TM), '\n'); flush.console() }
        return(selected_variables)
    })

    if(verbose){ cat(paste0('save the results in: ', normalizePath(output_dir),'\n')); flush.console() }
    # >>> save the results in the output_dir
    # check if this is not overwriting an existing file, if so, wait 1 sec and check again...
    while(TRUE){
        # timestamp: microSeconds since the reference date (1/1/1970)
        timestamp <- format(as.numeric(Sys.time())*1e6, scientific=F)
        pth <- paste0(output_dir, 
                      '/step1-results',
                      '-N-iters-', num_iterations, 
                      '-timestamp-', timestamp, 
                      '.rds')
        if(!file.exists(pth)) break
        Sys.sleep(runif(1)) # wait for some random time from 0 to 1
    }
    saveRDS(res, pth)

    F.TM <- Sys.time()
    if(verbose){ cat('elapsed time:', format(F.TM - S.TM), '\n'); flush.console() }
}

#' glmFS Step 2
#'
#' Run the second step of glmFS
#' @param step1_output_dir WIP
#' @param EN_cutoff WIP
#' @param anova_baseline WIP
#' @return WIP
#' @examples 
#' WIP
#' @export
glmFS_step2 <- function(step1_output_dir, 
                         EN_cutoff = 50,
                         anova_baseline = NA, 
                         verbose = FALSE){
    # if(plot_km) require(ggplot2)
    # ---------------------------------------------------------------------------------------------
    records <- readRDS(paste0(step1_output_dir, '/step1-records.rds'))
    family <- records$family
    data <- records$data
    data_vars <- records$data_vars
    features <- records$features
    # ---------------------------------------------------------------------------------------------
    S <- lapply(list.files(step1_output_dir, pattern = 'step1-results', full.names = T), readRDS)
    N <- sum(sapply(S, length))
    if(0 == length(unlist(S))){ if(verbose) cat('No feature selected by EN in step1\n'); return(NULL) }
    S <- sort(100 * table(unlist(S)) / N, decreasing = T)
    features$EN_score <- as.numeric(S[features$variable])
    features <- subset(features, EN_score > EN_cutoff)
    if(nrow(features) == 0){ if(verbose) cat('No feature passed EN_cutoff\n'); return(NULL) }
    features <- features[order(features$EN_score, decreasing=T),]
    # ---------------------------------------------------------------------------------------------
    new_variable <- function() paste0('v', 1 + max(as.numeric(gsub('v', '', grep('^v[0-9]+$', colnames(data_vars), value=T)))))
    get_variable_of_single_feature <- function(name, level){
        # use the original "features" from "records"
        features <- records$features
        df <- subset(features, id == name)
        if(!is.na(level)) df <- subset(df, feature1_level == level)
        if(0 == nrow(df)) return(NA)
        stopifnot(1 == nrow(df))
        return(df$variable)
    }
    # extract the variables of individual features and place them as separate columns in each row
    features$feature1_variable <- mapply(get_variable_of_single_feature, features$feature1, features$feature1_level)
    features$feature2_variable <- mapply(get_variable_of_single_feature, features$feature2, features$feature2_level)
    # remove rows with "interaction_element" from "features" to focus on rows with "single_features" & "interaction_features"
    features <- subset(features, type != 'interaction_element')
    if(nrow(features) == 0){ if(verbose) cat('No feature passed EN_cutoff\n'); return(NULL) }
    # ---------------------------------------------------------------------------------------------
    # validate anova_baseline
    if(identical(NA, anova_baseline)) anova_baseline <- character(0)
    if(!all(anova_baseline %in% colnames(data))){ stop('anova_baseline variables must be in "data"') }
    # make baseline_variables
    baseline_variables <- c()
    for(b in anova_baseline){
        if(! class(data[[b]]) %in% c('factor', 'numeric')) stop('anova_baseline variables must be numeric or factor in "data"')
        lv <- NA; if(class(data[[b]]) == 'factor') lv <- levels(data[[b]])[-1]
        if(length(lv) == 0) stop('a factor as an anova_baseline must have more than two levels')
        for(l in lv){
            v <- get_variable_of_single_feature(b, l)
            if(!is.na(v)){
                baseline_variables <- c(baseline_variables, v)
            } else {
                v <- new_variable()
                baseline_variables <- c(baseline_variables, v)
                if(is.na(l)){
                    data_vars[[v]] <- data[[b]]
                } else {
                    data_vars[[v]] <- as.numeric(data[[b]] == l)
                }
            }
        }
    }
    # ---------------------------------------------------------------------------------------------
    # add p anova and effect size columns to "features"
    features <- cbind(features, do.call(rbind, apply(features, 1, function(r){
        X <- r[['variable']]
        B <- baseline_variables
        V1 <- r[['feature1_variable']]
        V2 <- r[['feature2_variable']]
        isIntr <- grepl('\\*', r[['id']])
        data.frame(P_ANOVA_of_feature = p_anova(family, data_vars, B, X),
                   P_ANOVA_of_feature_controlled_for_first_variable = ifelse(isIntr, p_anova(family, data_vars, c(B, V1), X), NA),
                   P_ANOVA_of_feature_controlled_for_second_variable = ifelse(isIntr, p_anova(family, data_vars, c(B, V2), X), NA),
                   P_ANOVA_of_feature_controlled_for_both_variables = ifelse(isIntr, p_anova(family, data_vars, c(B, V1, V2), X), NA),
                   P_ANOVA_of_first_variable = ifelse(isIntr, p_anova(family, data_vars, B, V1), NA),
                   P_ANOVA_of_second_variable = ifelse(isIntr, p_anova(family, data_vars, B, V2), NA),
                   Effect_size_of_feature = effect_size(family, data_vars, B, X, 'value'),
                   Effect_size_lower.95_of_feature = effect_size(family, data_vars, B, X, 'lower'),
                   Effect_size_upper.95_of_feature = effect_size(family, data_vars, B, X, 'upper'),
                   Effect_size_of_first_variable = ifelse(isIntr, effect_size(family, data_vars, B, V1, 'value'), NA),
                   Effect_size_upper.95_of_first_variable = ifelse(isIntr, effect_size(family, data_vars, B, V1, 'upper'), NA),
                   Effect_size_lower.95_of_first_variable = ifelse(isIntr, effect_size(family, data_vars, B, V1, 'lower'), NA),
                   Effect_size_of_second_variable = ifelse(isIntr, effect_size(family, data_vars, B, V2, 'value'), NA),
                   Effect_size_lower.95_of_second_variable = ifelse(isIntr, effect_size(family, data_vars, B, V2, 'lower'), NA),
                   Effect_size_upper.95_of_second_variable = ifelse(isIntr, effect_size(family, data_vars, B, V2, 'upper'), NA))
    })))
    # ---------------------------------------------------------------------------------------------
    # if(plot_km){
    #     max_time <- max(data_vars$time)
    #     get_km_df <- function(fit, custom_levels){
    #         stopifnot(length(fit$strata) == length(custom_levels))
    #         df <- rbind(
    #             do.call(rbind, lapply(names(fit$strata), function(lv){
    #                 strata <- summary(fit)$strata
    #                 time <- summary(fit)$time[strata == lv]
    #                 surv <- summary(fit)$surv[strata == lv]
    #                 n <- sum(strata == lv)
    #                 data.frame(x = c(0, time, time), 
    #                            xend = c(time, max_time, time), 
    #                            y = c(1, surv, 1, surv[-n]), 
    #                            yend = c(1, surv, surv), 
    #                            strata = lv,
    #                            linewidth = 2)
    #             })), 
    #             subset(data.frame(x = fit$time, 
    #                               xend = fit$time, 
    #                               y = pmax(0, fit$surv - 0.01),
    #                               yend = pmin(1, fit$surv + 0.01),
    #                               strata = rep(names(fit$strata), fit$strata),
    #                               linewidth = 1,
    #                               n.censor = fit$n.censor), n.censor > 0, select=-n.censor))
    #         df$strata <- factor(custom_levels[match(df$strata, names(fit$strata))], levels = rev(custom_levels))
    #         return(df)
    #     }
    #     get_km_plot <- function(this_id){
    #         df <- subset(features, id == this_id)
    #         stopifnot(nrow(df) > 0)
    #         NR <- nrow(df)
    #         
    #         if(grepl('\\*', this_id) & any(is.na(df))) return(NULL)
    #         if(any(is.na(df$feature1_level))) return(NULL)
    #             
    #         
    #         custom_levels <- paste0('Feature ', c('absent', 'present'))
    #         colors <- c('black', 'red')
    #         names(colors) <- custom_levels
    #     
    #         km_df <- do.call(rbind, lapply(df$variable, function(v){
    #             df <- subset(features, variable == v)
    #             stopifnot(1 == nrow(df))
    #             
    #             isIntr <- grepl('\\*', df$id)
    #             fit <- survival::survfit(as.formula(paste0('survival::Surv(time, status) ~ ', v)), data_vars)
    #             km_df <- cbind(get_km_df(fit, custom_levels), 
    #                            title = paste0('Feature: ', gsub('\\*', ' * ', df$id), 
    #                                           '\nLevel: ', ifelse(isIntr, paste(df$feature1_level, '&', 
    #                                                                             df$feature2_level), 
    #                                                               df$feature1_level),
    #                                           '\nN. samples: ', sum(data_vars[[v]]), '/', nrow(data_vars),
    #                                           '\nlogHR: ', signif(df$Effect_size_of_feature, 2),
    #                                           '\nP: ', signif(df$P_ANOVA_of_feature, 2)))
    #             if(isIntr){
    #                 km_df <- rbind(km_df, do.call(rbind, lapply(1:2, function(i){
    #                     v <- df[[paste0('feature',i,'_variable')]]
    #                     n <- ifelse(i == 1, 'first', 'second')
    #                     fit <- survival::survfit(as.formula(paste0('survival::Surv(time, status) ~ ', v)), data_vars)
    #                     cbind(get_km_df(fit, custom_levels), 
    #                           title = paste0('Feature: ', df[[paste0('feature',i)]], 
    #                                          '\nLevel: ', df[[paste0('feature',i,'_level')]],
    #                                          '\nN. samples: ', sum(data_vars[[v]]), '/', nrow(data_vars),
    #                                          '\nlogHR: ', signif(df[[paste0('Effect_size_of_',n,'_variable')]], 2),
    #                                          '\nP: ', signif(df[[paste0('P_ANOVA_of_',n,'_variable')]], 2)))
    #                 })))
    #             }
    #             return(km_df)
    #         }))
    #         
    #         km_df$title <- factor(km_df$title, levels = unique(km_df$title))
    #         ggplot()+theme_bw()+
    #             geom_segment(data = km_df, aes(x=x, xend=xend, y=y, yend=yend, color=strata, linewidth=linewidth))+
    #             scale_linewidth_continuous(range = c(0.3, 0.6), guide = 'none')+
    #             scale_y_continuous(limits = c(0, 1))+
    #             scale_x_continuous(limits = c(0, max_time))+
    #             facet_wrap(~title, nrow = NR)+
    #             scale_color_manual(values = colors, name=NULL)+
    #             ylab('Survival P')+xlab('Time')+
    #             theme(panel.grid = element_blank(), 
    #                   strip.text.x = element_text(hjust=0), 
    #                   strip.background.x = element_blank())
    #     }
    #     km_plots <- sapply(unique(features$id), get_km_plot, simplify = FALSE, USE.NAMES = TRUE)
    # } else {
    #     km_plots <- NULL
    # }
    # ---------------------------------------------------------------------------------------------
    exclude_cols <- c('type', 'variable', 
                      'feature1', 'feature2',
                      'feature1_variable', 'feature2_variable')
    features <- features[, setdiff(colnames(features), exclude_cols)]
    colnames(features)[colnames(features) == 'id'] <- 'feature'
    colnames(features)[colnames(features) == 'feature1_level'] <- 'level_of_first_variable'
    colnames(features)[colnames(features) == 'feature2_level'] <- 'level_of_second_variable'
    rownames(features) <- NULL
    # ---------------------------------------------------------------------------------------------
    return(list(features = features))#, km_plots = km_plots))
}
