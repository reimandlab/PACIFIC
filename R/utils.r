do_fit <- function(family, data, variables){
    if(0 == length(variables)) variables <- '1'
    rhs <- paste(variables, collapse = '+')
    if(family == 'cox'){
        frm <- as.formula(paste0('survival::Surv(response, event) ~ ', rhs))
        fit <- suppressWarnings(survival::coxph(formula = frm, data = data))
    } else {
        frm <- as.formula(paste0('response ~ ', rhs))
        fit <- suppressWarnings(glm(formula = frm, data = data, family = family))
    }
    return(fit)
}

p_univariate <- function(family, data, X){
    coefs <- summary(do_fit(family, data, X))$coefficients
    if(! X %in% rownames(coefs)) return(1)
    n <- grep('Pr\\(', colnames(coefs))
    if(1 != length(n)) stop('fit coefficients are not as expected!')
    return(coefs[X, n])
    # this is Wald test for coxPH and logistic regression
    # https://stats.stackexchange.com/questions/59879/logistic-regression-anova-chi-square-test-vs-significance-of-coefficients-ano
    # but it is t-test for linear regression
}

p_anova <- function(family, data, cntrl, X){
    cntrl <- setdiff(cntrl, X)
    fit0 <- do_fit(family, data, cntrl)
    fit1 <- do_fit(family, data, c(cntrl, X))
    a <- anova(fit0, fit1, test = 'Chisq')
    n <- grep('Pr\\(', colnames(a))
    if(1 != length(n)) stop('anova() function is not working as expected!')
    p <- a[2, n]
    if(is.na(p)) p <- 1
    return(c('p.anova' = p))
}

effect_size <- function(family, data, cntrl, X, pick=c(NA, 'value', 'lower', 'upper')){
    pick <- match.arg(pick)
    fit <- do_fit(family, data, c(cntrl, X))
    v <- as.numeric(fit$coefficients[X])
    ci <- as.numeric(suppressMessages(confint(fit)[X,]))
    if(is.na(pick)) return(c('effect.size' = v, 'effect.size.95%.CI.lower' = ci[1], 'effect.size.95%.CI.upper' = ci[2]))
    if(pick == 'value') return(v)
    if(pick == 'lower') return(ci[1])
    if(pick == 'upper') return(ci[2])
}

discretize <- function(data, features, discretization_function){
    for(f in features){
        if(!is.numeric(data[[f]])) stop(paste0('"',f,'": cannot discretize non-numeric'))
        data[[f]] <- discretization_function(data[[f]])
    }; rm(f)
    return(data)
}

make_table_of_features <- function(data, features_vector, single_features, flexible_features){
    # this table contains information of feature-IDs: types, levels of elements, base groups, variables.
    features <- do.call(rbind, lapply(features_vector, function(id){
        get_levels <- function(f){
            # NA level for NA feature
            if(is.na(f)) return(NA)
            v <- data[[f]]
            # sanity check: v must be either numeric or factor with at least two levels
            stopifnot(is.numeric(v) | (is.factor(v) & length(levels(v)) >=2))
            # NA level for numeric feature
            if(is.numeric(v)) return(NA)
            # otherwise, input feature is factor, return levels
            return(levels(v))
        }
        get_bases <- function(f){
            # get levels of this feature
            L <- get_levels(f)
            # if levels=NA, return NA (meaning that this feature is NA or numeric)
            if(identical(NA, L)) return(NA)
            # otherwise, input feature is factor, return bases based on flexibility
            # if flexible, all levels can be potential bases
            if(f %in% flexible_features) return(L)
            # otherwise, the first level is the only base
            return(L[1])
        }
        sp <- unlist(strsplit(id, '\\*'))
        df <- expand.grid(id=id,
                          type = ifelse(2==length(sp), 
                                        'interaction_feature', 
                                        ifelse(id %in% single_features, 
                                               'single_feature', 
                                               'interaction_element')), 
                          feature1 = sp[1], feature1_level = get_levels(sp[1]),
                          feature2 = sp[2], feature2_level = get_levels(sp[2]),
                          stringsAsFactors = F)
        if(length(sp) == 1) sp <- c(sp, NA)
        ex <- expand.grid(lapply(sp, get_bases))
        do.call(rbind, lapply(1:nrow(ex), function(i){
            L1 <- ex[i,1]; if(!is.na(L1)) df <- subset(df, feature1_level != L1)
            L2 <- ex[i,2]; if(!is.na(L2)) df <- subset(df, feature2_level != L2)
            cbind(df, group = i)
        }))
    }))
    # assign variable names (avoid duplicated variables by taking unique ids of combined row values)
    uid <- apply(features[, !grepl('group', colnames(features))], 1, paste, collapse=',')
    features$variable <- paste0('v', as.numeric(factor(uid, levels = unique(uid))))
    # return the table of "features"
    return(features)
}

make_data_for_variables <- function(data, features){
    # make a new data table with the same {response, event} columns
    new_data <- data[, intersect(colnames(data), c('response','event')), drop=F]
    # then start adding new "variable" columns to new_data
    # these "variables" map to numeric columns based on "features" table:
    # binary when the feature contains only factor(s) {fac, fac*fac}
    # continuous when the feature contains at least one numeric element {num, num*num, num*fac, fac*num}
    for(i in 1:nrow(features)){
        r <- features[i,]
        get_v <- function(name, level){
            if(is.na(name)) return(1)
            v <- data[[name]]
            if(!is.na(level)) v <- as.numeric(v == level)
            return(v)
        }
        v1 <- get_v(r$feature1, r$feature1_level)
        v2 <- get_v(r$feature2, r$feature2_level)
        new_data[[r$variable]] <- v1 * v2
    }
    return(new_data)
}

filter_by_sparsity <- function(data, features, skip_ids, min_group_size){
    # loop over unique ids in "features" and keep the ids if either of the following is true:
    # 1 - pass the sparsity filter
    # 2 - included in "skip_ids"
    keep_ids <- names(which(sapply(unique(features$id), function(f){
        # keep this id if in "skip_ids"
        if(f %in% skip_ids) return(TRUE)
        # extract the factors from elements of this feature
        fct <- names(which(sapply(data[, unlist(strsplit(f, '\\*')), drop=F], is.factor)))
        # if this id does not contain any factor feature (i.e. num, or num*num), just keep it
        if(length(fct) == 0) return(TRUE)
        # otherwise (i.e. if there are extracted factors), keep this id only if the smallest group size > cutoff
        return(min(table(data[, fct, drop=F])) > min_group_size)
    })))
    # if an interaction feature passes the sparsity filter, its corresponding single features must have passed too
    # the other side of this should be resolved in a second pass:
    # if an interaction feature is filtered by sparsity, its corresponding single features must be removed too unless they 
    # contribute to other interaction features that have passed the filter. 
    # for this, remove ids that remain "alone" when splitting (by "*") the ids of "keep_ids" after excluding "single_features".
    single_features <- subset(features, type=='single_feature')$id
    rm_ids <- names(which(1 == table(unlist(strsplit(setdiff(keep_ids, single_features), '\\*')))))
    keep_ids <- setdiff(keep_ids, rm_ids)
    # subset to ids passing the filter, and return the filtered features.
    features <- subset(features, id %in% keep_ids)
    return(features)
}

filter_by_univariate_association <- function(family, data, features, skip_ids, P_cutoff){
    # calculate univariate P values for all variables in "features"
    features$univariate_P <- sapply(features$variable, function(v){
        p <- p_univariate(family, data, v)
        if(is.na(p)) p <- 1
        stopifnot(is.finite(p))    # ensure appropriate value of P
        stopifnot(p >= 0 & p <= 1) # ensure appropriate value of P
        return(p)
    })        
    pass_filter <- function(this_id){
        # process the variables in "features" subset with <this_id>
        df <- subset(features, id == this_id)
        # choose the group with "the most optimum" P values
        x <- do.call(rbind, lapply(unique(df$group), function(g){
            df <- subset(df, group == g)
            data.frame(group = g, N_pass = sum(df$univariate_P < P_cutoff), best_univariate_P = min(df$univariate_P))
        }))
        x <- x[x$N_pass == max(x$N_pass), ]
        if(length(unique(x$group)) > 1) x <- x[x$best_univariate_P == min(x$best_univariate_P), ]
        chosen_group <- unique(x$group)
        # stopifnot(1 == length(chosen_group))
        if(length(chosen_group) > 1){
            # in the case of tie, always choose the first group
            chosen_group <- chosen_group[1]
        }
        # subset to the chosen group
        df <- subset(df, group == chosen_group)
        # return the whole group if <this_id> is in "skip_ids"
        if(this_id %in% skip_ids) return(df)
        # otherwise, return the subset of chosen group where univariate P < "P_cutoff"
        return(subset(df, univariate_P < P_cutoff))
    }
    # first, apply filter to "single_features" & "interaction_features", captured by rows with type != interaction_element
    subset_1 <- do.call(rbind, lapply(subset(features, type != 'interaction_element')$id, pass_filter))    
    # then add the individual variables that contribute to the interaction features that passed the filter
    extract_single_feature <- function(name, level){
        if(is.na(name)) return(NULL)
        if(is.na(level)) return(subset(features, id == name & is.na(feature1_level)))
        return(subset(features, id == name & feature1_level == level))
    }
    subset_2 <- subset(subset_1, grepl('\\*', id))
    if(nrow(subset_2) > 0){
        subset_2 <- do.call(rbind, apply(subset_2, 1, function(r){
            rbind(extract_single_feature(r[['feature1']], r[['feature1_level']]), 
                  extract_single_feature(r[['feature2']], r[['feature2_level']]))
        }))
    }
    # rbind subset_1 & subset_2 to get the updated (filtered) features table. Do unique, some rows might be duplicated.
    features <- unique(rbind(subset_1, subset_2))
    return(features)
}

elastic_net_variable_selection <- function(family, data, variables){
    if(family == 'cox'){
        Y <- as.matrix(data[, c('response', 'event')])
        colnames(Y) <- c('time', 'status')
    } else {
        Y <- as.matrix(data[, 'response', drop=F])
    }
    X <- as.matrix(data[, variables, drop=F])
    fit <- suppressWarnings(glmnet::cv.glmnet(x=X, y=Y, family=family, alpha=0.5))
    coefs <- as.matrix(glmnet::coef.glmnet(fit, s = 'lambda.min'))
    selected_variables <- rownames(coefs)[which(coefs[,1] != 0)]
    return(selected_variables)
}

wait_for_file_to_stabilize <- function(p){
    while(TRUE){ 
        s <- file.info(p)$size
        Sys.sleep(3)
        if(s == file.info(p)$size) break 
    } 
}

