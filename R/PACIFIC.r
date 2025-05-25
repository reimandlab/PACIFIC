#' Run step 1 of PACIFIC for survival outcomes
#' 
#' @param data data.frame. The table of input data with samples in rows and features in columns. 
#' @param baseline character vector. The baseline clinical features as columns of "data". Default: c("age", "sex", "stage", "grade"). Put NA for empty set.
#' @param feat1 character vector. The feature set 1 as columns of "data".
#' @param feat2 character vector. The feature set 2 as columns of "data". Default: NA (for empty set)
#' @param discretization_method character. Continuous features are binarized using this method within the feature selection loops. Default: "median"
#' @param sparsity_threshold character. Sparse features are called using this threshold within the feature selection loops. Default: "5_percent"
#' @param univariate_p_cutoff numeric. Features weakly associated with survival are called using this cutoff for P-value of univariate CoxPH models within the feature selection loops. Default: 0.1
#' @param subsampling_ratio numeric. Ratio used for subsampling the "data" for the feature selection loops. Default: 0.8
#' @param num_iterations numeric. How many iterations of the feature selection loop to run. Default: 10
#' @param output_dir character. The directory to store the results.
#' @param verbose character. Whether to print progress messages. Default: FALSE
#' 
#' @export
PACIFIC_survival_step1 <- function(data, 
                                   baseline = NA,
                                   feat1, 
                                   feat2, 
                                   discretization_method = "median", 
                                   sparsity_threshold = "5_percent", 
                                   univariate_p_cutoff = 0.1, 
                                   subsampling_ratio = 0.8, 
                                   num_iterations = 10, 
                                   EN_cutoff = 50,
                                   output_dir, 
                                   verbose = TRUE){
    
    # validate and process arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    # validate arg data
    if(! is.data.frame(data)) stop('"data" must be a data.frame or an instance of a class extended from data.frame')
    if(! 'time' %in% colnames(data)) stop('"data" must contain column "time" for survival time value')
    if(! 'status' %in% colnames(data)) stop('"data" must contain column "status" for survival event status')
    
    # validate arg baseline
    if(identical(NA, baseline)) baseline <- character(0)
    if(! is.character(baseline)) stop('"baseline" must be a character vector')
    if(! all(baseline %in% colnames(data))) stop('all "baseline" features must be in "data" columns')
    
    # validate arg feat1
    if(! is.character(feat1)) stop('"feat1" must be a character vector')
    if(length(feat1) == 0) stop('"feat1" must not be empty')
    if(! all(feat1 %in% colnames(data))) stop('all "feat1" features must be in "data" columns')
    
    # validate arg feat2
    if(! is.character(feat2)) stop('"feat2" must be a character vector')
    if(length(feat2) == 0) stop('"feat2" must not be empty')
    if(! all(feat2 %in% colnames(data))) stop('all "feat2" features must be in "data" columns')
    
    # set single_features and interaction_features
    single_features <- baseline
    interaction_features <- list(feat1, feat2)
    
    # prefiltering (of any type) is only skipped for baseline
    features_to_skip_sparsity_prefiltering <- baseline
    features_to_skip_univariate_association_prefiltering <- baseline
    
    # any numeric feature will go to features_to_discretize & features_with_flexible_direction
    # any other feature must be binary factor
    features_to_discretize <- character(0)
    features_with_flexible_direction <- character(0)
    for(f in c(feat1, feat2)){
        values <- data[[f]]
        if(is.numeric(values)){
            features_to_discretize <- c(features_to_discretize, f)
            features_with_flexible_direction <- c(features_with_flexible_direction, f)
        } else {
            err <- 'features in feat1 & feat2 must be either numeric or factor with two levels'
            if(! is.factor(values)) stop(err)
            if(length(levels(values)) != 2) stop(err)
        }
    }
    
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    PACIFIC_step1(data = data,
                  response = 'time',
                  event = 'status',
                  single_features = single_features,
                  interaction_features = interaction_features,
                  features_to_discretize = features_to_discretize,
                  features_with_flexible_direction = features_with_flexible_direction,
                  features_to_skip_sparsity_prefiltering = features_to_skip_sparsity_prefiltering,
                  features_to_skip_univariate_association_prefiltering = features_to_skip_univariate_association_prefiltering,
                  num_iterations = num_iterations, 
                  output_dir = output_dir, 
                  verbose = verbose)
}



#' Run step 2 of PACIFIC for survival outcomes
#' 
#' @export
PACIFIC_survival_step2 <- function(output_dir, 
                                   EN_cutoff = 50,
                                   anova_baseline = NA, 
                                   verbose = TRUE){
    
    S.TM <- Sys.time()
    if(verbose){ cat('----------------------------------------------------\nPACIFIC step 2:\n'); flush.console() }
    step2 <- PACIFIC_step2(output_dir = output_dir, 
                           anova_baseline = anova_baseline, 
                           EN_cutoff = EN_cutoff)
    # -------------------------------------------------------------------------------------------
    features <- step2$features
    # -------------------------------------------------------------------------------------------
    if(is.null(features)){
        cat('No interaction feature found\n')
        return(NULL)
    }
    # -------------------------------------------------------------------------------------------
    features <- features[grepl('\\*', features$id),]
    # -------------------------------------------------------------------------------------------
    if(nrow(features) == 0){
        cat('No interaction feature found\n')
        return(NULL)
    }
    # -------------------------------------------------------------------------------------------
    records <- step2$records
    data_vars <- records$data_vars
    colnames(data_vars)[colnames(data_vars) == 'response'] <- 'time'
    colnames(data_vars)[colnames(data_vars) == 'event'] <- 'status'
    km_plot_list <- sapply(unique(features$id), get_km_plot, features=features, data_vars=data_vars, simplify = FALSE, USE.NAMES = TRUE)
    # -------------------------------------------------------------------------------------------
    rownames(features) <- NULL
    
    exclude_cols <- c('type', 'variable',
                      'feature1', 'feature2',
                      'feature1_variable', 'feature2_variable')
    features <- features[, setdiff(colnames(features), exclude_cols)]
    
    colnames(features)[colnames(features) == 'id'] <- 'intr'
    colnames(features)[colnames(features) == 'feature1_level'] <- 'feat1_level'
    colnames(features)[colnames(features) == 'feature2_level'] <- 'feat2_level'
    
    colnames(features)[colnames(features) == 'P_ANOVA_of_feature'] <- 'intr_P'
    colnames(features)[colnames(features) == 'P_ANOVA_of_feature_controlled_for_first_variable'] <- 'intr_P_C1'
    colnames(features)[colnames(features) == 'P_ANOVA_of_feature_controlled_for_second_variable'] <- 'intr_P_C2'
    colnames(features)[colnames(features) == 'P_ANOVA_of_feature_controlled_for_both_variables'] <- 'intr_P_C3'
    
    colnames(features)[colnames(features) == 'P_ANOVA_of_first_variable'] <- 'feat1_P'
    
    colnames(features)[colnames(features) == 'P_ANOVA_of_second_variable'] <- 'feat2_P'
    
    colnames(features)[colnames(features) == 'Effect_size_of_feature'] <- 'intr_logHR'
    colnames(features)[colnames(features) == 'Effect_size_lower.95_of_feature'] <- 'intr_logHR_lower95'
    colnames(features)[colnames(features) == 'Effect_size_upper.95_of_feature'] <- 'intr_logHR_upper95'
    
    colnames(features)[colnames(features) == 'Effect_size_of_first_variable'] <- 'feat1_logHR'
    colnames(features)[colnames(features) == 'Effect_size_lower.95_of_first_variable'] <- 'feat1_logHR_lower95'
    colnames(features)[colnames(features) == 'Effect_size_upper.95_of_first_variable'] <- 'feat1_logHR_upper95'
    
    colnames(features)[colnames(features) == 'Effect_size_of_second_variable'] <- 'feat2_logHR'
    colnames(features)[colnames(features) == 'Effect_size_lower.95_of_second_variable'] <- 'feat2_logHR_lower95'
    colnames(features)[colnames(features) == 'Effect_size_upper.95_of_second_variable'] <- 'feat2_logHR_upper95'
    
    # -------------------------------------------------------------------------------------------
    results <- list(top_interactions=features, km_plot_list=km_plot_list)
    
    saveRDS(results, paste0(output_dir, '/aggregated-results.rds'))
    
    F.TM <- Sys.time()
    if(verbose){ cat('elapsed time:', format(F.TM - S.TM), '\n'); flush.console() }
    return(results)
}




# data A data.frame (or an extension of data.framethe, e.g. data.table). The table of input data with rows for samples and columns for features.
# response Name of the column in data which contains the response values.
# event Name of the column in 'data' which contains the event status values for survival outcomes. This column must contain only 0 and 1. Default is 'NA'.
# single_features  (for input features) Vector of "single features". Default is 'NA' (for no single feature).
# interaction_features (for input features) Explicit or list-based specification of "interaction features". Explicit: a vector of the form 'c("A*B", "C*D", ...)'. List-based: list of two vectors of single components for all two-way combinations as interaction features. Default is 'NA' (for no interaction feature). Note that 'single_features' cannot overlap with any component in 'interaction_features'.
# features_to_discretize The subset of input features (being single or component of interaction) to be discretized, i.e. to be converted from numeric to categorical (subject to 'discretization_method'). Default is 'NA' (for no such feature).
# discretization_method The method for discretization. Either "median" (for median dichotomization), or a 'function' designed to get a numeric vector as input and return a factor vector (of the same length) as output. Default is '"median"'.
# features_with_flexible_direction The subset of input features (being single or component of interaction) to be allowed to have flexible direction. Applicable only to either factor variables, or numeric variables in 'features_to_discretize'. This means that the reference level of these variables are to be determined by the pipeline in order to optimize the associations with the 'response' Default is 'NA' (for no such feature).
# features_to_skip_sparsity_prefiltering The subset of input features (being single only) to be skipped for prefiltering by sparsity (within iterations). Default is 'NA' (for no such feature).
# features_to_skip_univariate_association_prefiltering The subset of input features (being single only) to be skipped for prefiltering by weak associations with the 'response' (within iterations) . Default is NA (for no such feature).
# sparsity_criterion The criterion used for detecting sparse features in prefiltering by sparsity (within iterations). Applicable to features which involve factor(s). It should be a string of the form '"X_percent"' or '"X_absolute"' (where 'X' is a number), meaning that the cutoff for the minimum size of the feature is either 'X' percent of the size of the input data (within the iteration) or is just 'X', respectively. Default is '"5_percent"'.
# univariate_p_cutoff Cutoff for the P values of univariate models used for prefiltering by weak association with the 'response'. Default is 0.1.
# subsampling_ratio The ratio used for subsampling 'data' to make the input data for each iteration. Note that the ratio of events (in survival analysis) or each level of response (if it's categorical) in 'data' is maintained in the subsamplings. Default is '0.8'.
# num_iterations Number of iterations to run. Default is '10'.
# output_dir The directory to store the results of iterations.
# verbose Whether to print progress messages. Default is 'FALSE'.
PACIFIC_step1 <- function(data,
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
                          num_iterations,
                          output_dir,
                          verbose = FALSE){
    
    UNIQUE.TAG <- gsub('\\/', '', tempfile(pattern='', tmpdir=''))
    
    S.TM <- Sys.time()
    
    if(verbose){ cat('----------------------------------------------------\nPACIFIC step 1:\n'); flush.console() }
    if(verbose){ cat('preprocessing ... '); flush.console() }
    
    # validate and process arguments >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    # create the output_dir. If it exists, this does nothing
    dir.create(output_dir, recursive = T, showWarnings = F)
    
    # save the list of input arguments in output_dir
    input_args <- sapply(names(formals(sys.function())), get, envir = environment(), simplify = F, USE.NAMES = T)
    saveRDS(input_args, paste0(output_dir, '/step1-arguments', '-tag-', UNIQUE.TAG, '.rds'))
    
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
            if(!all(unlist(sp) %in% avDtCls)) stop('"interaction_features" must be subset of available data columns')
            if(all(2 == sapply(sp, length))){
                if(any(2 != sapply(sapply(sp, unique, simplify = F), length))) stop('found self-interaction in "interaction_features"')
                if(any(duplicated(t(sapply(sp, sort))))) stop('found duplicated pair in "interaction_features"')
            } else if(all(1 == sapply(sp, length))){
                sp <- unlist(sp)
                interaction_features <- apply(combn(sp, 2), 2, paste, collapse='*')
            } else {
                stop('invalid input "interaction_features"')
            }
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
    
    
    # >>> create the "records" and save in output_dir. 
    # "records" = list("family", "data", "data_vars", "features")
    # to save in "records", remove "group" column from "features" and unique its rows
    temp_features <- unique(subset(features, select = -group))
    rownames(temp_features) <- NULL
    records <- list(family = family, data = data, data_vars = data_vars, features = temp_features)
    rm(temp_features)
    saveRDS(records, paste0(output_dir, '/step1-records', '-tag-', UNIQUE.TAG, '.rds'))
    
    
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
    
    
    S.TM.saving <- Sys.time()
    if(verbose){ cat('saving the results ... '); flush.console() }
    saveRDS(res, paste0(output_dir, '/step1-results', '-iter-', num_iterations, '-tag-', UNIQUE.TAG, '.rds'))
    F.TM.saving <- Sys.time()
    if(verbose){ cat(format(F.TM.saving - S.TM.saving), '\n'); flush.console() }
    
    
    F.TM <- Sys.time()
    if(verbose){ cat('total elapsed time:', format(F.TM - S.TM), '\n'); flush.console() }

    return(TRUE)
}

# output_dir WIP
# EN_cutoff WIP
# anova_baseline WIP
PACIFIC_step2 <- function(output_dir, 
                          EN_cutoff = 50,
                          anova_baseline = NA, 
                          verbose = TRUE){
    
    # ---------------------------------------------------------------------------------------------
    read_first_and_check_with_others <- function(fs){
        r <- readRDS(fs[1])
        for(f in fs[-1]){
            if(!identical(r, readRDS(f))){
                stop('Step1 calls are inconsistent.')
            }
        }
        return(r)
    }
    # ---------------------------------------------------------------------------------------------
  
    # if(verbose) cat('Checking the consistency of Step1 calls\n'); flush.console()
    # # --------------------
    # fs <- list.files(output_dir, pattern = 'step1-arguments', full.names = T)
    # args <- readRDS(fs[1])
    # args <- args[! names(args) %in% c('num_iterations', 'verbose')]
    # for(f in fs[-1]){
    #     r <- readRDS(f)
    #     r <- r[! names(r) %in% c('num_iterations', 'verbose')]
    #     if(!identical(r, args)){
    #         stop('Step1 calls are inconsistent.')
    #     }
    # }
    # # --------------------
  
    fs <- list.files(output_dir, pattern = 'step1-records', full.names = T)
    records <- readRDS(fs[1])
    
    # for(f in fs[-1]){
    #     r <- readRDS(f)
    #     if(!identical(r, records)){
    #         stop('Step1 calls are inconsistent.')
    #     }
    # }
    # # --------------------
  
    family <- records$family
    data <- records$data
    data_vars <- records$data_vars
    features <- records$features
    # ---------------------------------------------------------------------------------------------
    if(verbose) cat('Processing the results of Step1 calls\n'); flush.console()
    S <- lapply(list.files(output_dir, pattern = 'step1-results', full.names = T), readRDS)
    N <- sum(sapply(S, length))
    if(0 == length(unlist(S))){ if(verbose) cat('No feature selected by EN in step1\n'); flush.console(); return(NULL) }
    S <- sort(100 * table(unlist(S)) / N, decreasing = T)
    features$EN_score <- as.numeric(S[features$variable])
    features <- subset(features, EN_score > EN_cutoff)
    if(nrow(features) == 0){ if(verbose) cat('No feature passed EN_cutoff\n'); flush.console(); return(NULL) }
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
    if(nrow(features) == 0){ if(verbose) cat('No feature passed EN_cutoff\n'); flush.console(); return(NULL) }
    # ---------------------------------------------------------------------------------------------
    # validate anova_baseline
    if(identical(NA, anova_baseline)) anova_baseline <- character(0)
    if(!all(anova_baseline %in% colnames(data))){ stop('anova_baseline variables must be in "data"') }
    # make baseline_variables
    baseline_variables <- c()
    for(b in anova_baseline){
        if(! class(data[[b]]) %in% c('factor', 'numeric')) stop('anova_baseline variables must be numeric or factor in "data"')
        lv <- NA; if(class(data[[b]]) == 'factor') lv <- levels(data[[b]])[-1]
        if(length(lv) == 0) stop('a factor as an anova_baseline must have more than one level')
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
        suppressWarnings(
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
                       Effect_size_lower.95_of_first_variable = ifelse(isIntr, effect_size(family, data_vars, B, V1, 'lower'), NA),
                       Effect_size_upper.95_of_first_variable = ifelse(isIntr, effect_size(family, data_vars, B, V1, 'upper'), NA),
                       
                       Effect_size_of_second_variable = ifelse(isIntr, effect_size(family, data_vars, B, V2, 'value'), NA),
                       Effect_size_lower.95_of_second_variable = ifelse(isIntr, effect_size(family, data_vars, B, V2, 'lower'), NA),
                       Effect_size_upper.95_of_second_variable = ifelse(isIntr, effect_size(family, data_vars, B, V2, 'upper'), NA))
        )
    })))
    
    return(list(features=features, records=records))
}

# Get the total number of iterations currently stored in a directory
# dir WIP
current_total_iters <- function(dir){
    r <- list.files(dir, pattern = 'step1-results')
    r <- regmatches(r, regexec('-iter-\\s*(.*?)\\s*-tag-', r))
    N <- sum(as.numeric(sapply(r, '[[', 2)))
    return(N)
}
