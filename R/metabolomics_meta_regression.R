#' Metabolomics meta-regression
#' 
#' Perform meta-regression for repeated measurements.
#' Return merged timewise and training differential analysis summary statistics, 
#' where results for features with multiple measurements are replaced with the
#' meta-regression results when appropriate. This method was used to generate the 
#' [MotrpacRatTraining6moData::METAB_DA_METAREG] differential analysis results,
#' which are the version of results used in the manuscript analyses. 
#'
#' @param tissue `r tissue()`
#' @param timewise_input r data frame, custom input. To see the expected format, look at 
#'   a table returned by [load_metab_da()] with \code{type="timewise"}. 
#' @param training_input r data frame, custom input. To see the expected format, look at 
#'   a table returned by [load_metab_da()] with \code{type="training"}. 
#' @param het_p_threshold numeric, meta-regression cases with a heterogeneity p-value 
#'   below this are considered to have high heterogeneity. Default: 0.001
#'
#' @return named list where \code{meta_reg_timewise_dea} is a data frame with 
#'   the adjusted timewise results, \code{training_meta_regression} is a data frame with the adjusted training
#'   results, \code{meta_regression_results} is a named list with meta-regression 
#'   results per redundant metabolite, \code{meta_regression_models} is a table of the number of each 
#'   type of model used for meta-regression, and \code{metabolite_categories} is a named list 
#'   of the RefMet IDs of metabolites corresponding to each category described in the details. 
#' 
#' @export
#' 
#' @seealso [load_metab_da()], [metabolite_meta_regression()], [MotrpacRatTraining6moData::METAB_DA_METAREG]
#' @details
#' 
#' We try multiple models per repeated analyte:  
#' \itemize{
#'   \item{Model 1:}{Two random effects factors if \code{platform} and \code{is_targeted} are not redundant. Default optimization.}
#'   \item{Model 2:}{Two random effects factors if \code{platform} and \code{is_targeted} are not redundant with alternative optimization.}
#'   \item{Model 3:}{\code{platform} and \code{is_targeted} are redundant. Use a single RE factor with default optimization. 
#'     Also, use this if \code{QMp} is NA, which is an indication of over-parameterization of the model.}
#'   \item{Model 4:}{\code{platform} and \code{is_targeted} are redundant. Use a single RE factor with alternative optimization.}
#'   \item{Model 5:}{If all previous analyses failed, use a simple fixed-effects approach.}
#' }
#' 
#' After performing meta-regression, we define four categories of metabolites:  
#' 
#' 1. measured once  
#' 2. measured multiple times, high heterogeneity, has a targeted platform  
#' 3. measured multiple times, high heterogeneity, no targeted platform  
#' 4. measured multiple times, low heterogeneity  
#' 
#' For categories 1 and 3 we keep the results as is. 
#' For category 2 we take the targeted data only. 
#' Finally, for category 4 we take the meta-regression results.
#'
#' @examples
#' # Perform meta-regression for gastrocnemius
#' res = metab_meta_regression("SKM-GN")
#' names(res)
metab_meta_regression = function(tissue, timewise_input = NULL, training_input = NULL, het_p_threshold = 0.001){
  
  if(!requireNamespace("metafor", quietly = TRUE)){
    stop(
      "Package 'metafor' must be installed to perform meta-regression.",
      call. = FALSE
    )
  }
  if (!requireNamespace("nloptr", quietly = TRUE)) {
    stop(
      "Package 'nloptr' must be installed to perform meta-regression'.",
      call. = FALSE
    )
  }
  
  ##############################################################################
  # load the original differential analysis summary statistics
  ##############################################################################
  
  if(tissue == "BLOOD"){
    stop("Metabolomics was not performed in whole blood. Did you mean tissue='PLASMA'?")
  }
  
  req_input_cols = c("logFC","logFC_se","metabolite_refmet","site","sex","comparison_group","dataset","feature_ID","p_value")
  # req_input_cols = c("feature_ID","assay","tissue","dataset","site","is_targeted","sex",
  #                    "comparison_group","metabolite_refmet",
  #                    "effect_size","effect_size_se",
  #                    "comparison_average_intensity","reference_average_intensity",
  #                    "p_value","cv")
  
  # load redundant timewise differential analysis results 
  if(is.null(timewise_input)){
    timewise = load_metab_da(tissue, type="timewise")
  }else{
    timewise = timewise_input
  }
  
  # load redundant training differential analysis results 
  if(is.null(training_input)){
    training = load_metab_da(tissue, type="training")
    setnames(training, c("tissue","tissue_abbreviation"), c("tissue_code","tissue"))
  }else{
    training = training_input
  }
  req_training_cols = c("metabolite_refmet","p_value","mz","rt","neutral_mass","p_value_male","p_value_female","fscore_male","fscore_female","site","dataset")
  if(!all(req_training_cols %in% colnames(training))){
    stop(sprintf("The following columns are expected in the training input:\n %s", paste0(req_input_cols, collapse=", ")))
  }
  
  # check/change col names
  if("effect_size" %in% colnames(timewise)){
    setnames(timewise, c("effect_size"), c("logFC"))
  }
  if("effect_size_se" %in% colnames(timewise)){
    setnames(timewise, c("effect_size_se"), c("logFC_se"))
  }
  if(!all(req_input_cols %in% colnames(timewise))){
    stop(sprintf("The following columns are expected in the timewise input:\n %s", paste0(req_input_cols, collapse=", ")))
  }

  # get redundant metabolites only 
  tb = table(timewise$metabolite_refmet,timewise$dataset)
  non_unique_metabs = rownames(tb)[rowSums(tb>0)>1]
  # check if we have ref standards
  if(sum(grepl("standard",non_unique_metabs,ignore.case = T))>0){
    warning(paste(tissue,"differential analysis results have ref standard results. Excluding them from meta-analysis."))
    non_unique_metabs = non_unique_metabs[!grepl("standard",non_unique_metabs,ignore.case = T)]
  }
  timewise = timewise[timewise$metabolite_refmet %in% non_unique_metabs,]

  ##############################################################################
  # perform meta-regression
  ##############################################################################
  
  message(sprintf("Performing meta-regression for %s redundant metabolites in %s...", length(non_unique_metabs), tissue))
  # for each metabolite, perform meta-regression 
  meta_res = lapply(non_unique_metabs,
                    metabolite_meta_regression, 
                    timewise = timewise,
                    training = training)
  names(meta_res) = non_unique_metabs
  message("Done.\n")
  
  ##############################################################################
  # provide summary of results 
  ##############################################################################
  
  #message(sprintf("Computed %s meta-regression results for %s.", length(meta_res), tissue))
  
  # plot of heterogeneity p-value versus meta-analysis p-value 
  het_p = unlist(sapply(meta_res,function(x)x$het_p_value))
  max_cv = unlist(sapply(meta_res,function(x)max(x$rawdata$control_cv,na.rm=T)))
  metaanal_p = unlist(sapply(meta_res,function(x)x$p_value))
  cols = rep("black",length(metaanal_p))
  cols[metaanal_p < 0.001] = "red"
  plot(x=-log10(het_p),y=-log10(metaanal_p),pch=20,xlab = "-log10 het_p",ylab = "meta-analysis p",col=cols)
  graphics::legend(x="topright",c("meta p<0.001","meta p>0.001"),fill=c("red","black"))
  
  # We now go over all meta-regression analysis objects and count the number of models that were used.
  # For example, a metabolite that is measured by two untargeted platforms will only have a single random effects component based on the platform used. 
  # Moreover, in some rare cases the random effects optimization process fails and we resort to running the fixed effects analysis.
  
  # Check which RE models were used directly from the call:
  metareg_calls = unlist(sapply(meta_res,function(y)as.character(y$call)[5]))
  # check number of platforms per meta-analysis
  metareg_nplatform = unlist(sapply(meta_res,function(y)y$nplatform))
  metareg_run_table = table(metareg_calls,metareg_nplatform)
  rownames(metareg_run_table)[3] = "FE"
  
  message(sprintf("Number of models that were fit for %s:", tissue))
  print(metareg_run_table)
  
  # look at distributions of heterogeneity p-values 
  
  # We show the distribution and count the number of cases with a p-value below the specified threshold.
  hist(het_p,main="Meta-regression heterogeneity p-values",xlab="")
  message(sprintf("Total number of cases with high heterogeneity: %s", sum(het_p < het_p_threshold)))
  message(sprintf("Total number of cases with low heterogeneity: %s", sum(het_p > het_p_threshold)))
  
  ##############################################################################
  # determine 4 categories of metabolites
  ##############################################################################
  
  # We now define four classes of metabolites:
  # (1) measured once
  # (2) measured multiple times, high heterogeneity, has a targeted platform
  # (3) measured multiple times, high heterogeneity, no targeted platform
  # (4) measured multiple times, low heterogeneity
  # 
  # For classes 1 and 3 we keep the results as is. 
  # For class 2 we take the targeted data only. 
  # Finally, for class 4 we take the meta-regression results.
  
  basic_cols_to_keep = c("feature_ID","assay","tissue","dataset","site","is_targeted",
                         "metabolite_refmet","cv","covariates","control_cv")
  additional_training_cols = c(
    "fscore_male","fscore_female",
    "p_value_male","p_value_female",
    "mz","rt","neutral_mass"
  )
  
  # class 1
  unique_metabs_inds = !(training$metabolite_refmet %in% names(meta_res))
  unique_metabs = unique(training$metabolite_refmet[unique_metabs_inds])
  # examine the heterogeneity
  het_ps = unlist(sapply(meta_res,function(x)x$het_p_value))
  high_het_metabs = names(meta_res)[het_ps < het_p_threshold]
  # get class 2
  high_het_metabs_targeted = c()
  for(m in high_het_metabs){
    if(any(training$is_targeted[training$metabolite_refmet == m],na.rm=T)){
      high_het_metabs_targeted = c(high_het_metabs_targeted,m)
    }
  }
  # get class 3
  high_het_metabs_untargeted_only = setdiff(high_het_metabs,high_het_metabs_targeted)
  # class 4, meta-analysis results
  meta_anal_metabs = setdiff(names(meta_res),high_het_metabs)
  metab_classes = list(
    unique_metabs = unique_metabs,
    high_het_metabs_targeted = high_het_metabs_targeted,
    high_het_metabs_untargeted_only = high_het_metabs_untargeted_only,
    meta_anal_metabs = meta_anal_metabs
  )

  message(sprintf("\nSummary of the number of %s metabolites in each category:", tissue))
  class_size_sum_table = t(t(sapply(metab_classes,length)))
  colnames(class_size_sum_table) = tissue
  print(class_size_sum_table)
  
  ##############################################################################
  # make the timewise differential analysis results tables 
  ##############################################################################
  
  unique_metabs_inds = !(timewise$metabolite_refmet %in% names(meta_res))
  # change the indices of some high het cases to TRUE here
  # this will make the algorithm add them
  tissue_metab_classes = metab_classes
  if(length(tissue_metab_classes[["high_het_metabs_targeted"]]) > 0){
    unique_metabs_inds[
      timewise$metabolite_refmet %in% tissue_metab_classes[["high_het_metabs_targeted"]] &
        timewise$is_targeted
    ] = TRUE
  }
  if(length(tissue_metab_classes[["high_het_metabs_untargeted_only"]]) > 0){
    unique_metabs_inds[
      timewise$metabolite_refmet %in% tissue_metab_classes[["high_het_metabs_untargeted_only"]]
    ] = TRUE
  }
  # finally, for the meta_reg obj keep the low het metabs only
  meta_res = meta_res[tissue_metab_classes[["meta_anal_metabs"]]]
  meta_res_coeff_tables = lapply(meta_res,function(x)stats::coefficients(summary(x)))
  
  curr_timewise_df = c()
  for(sex in unique(timewise$sex)){
    if(length(meta_res)==0){next}
    for(comp_g in unique(timewise$comparison_group)){
      curr_rowname = paste0("analysis_group",sex,",",comp_g)
      # transform results to a summary data frame
      basic_df = t(sapply(meta_res,function(x,y)unlist(x[y]),y=basic_cols_to_keep))
      basic_df = as.data.frame(basic_df,stringsAsFactors = F)
      basic_df$sex = sex
      basic_df$comparison_group = comp_g
      basic_df$covariates = NA
      basic_df$zscore = sapply(meta_res_coeff_tables,function(x,y)x[y,"zval"],y=curr_rowname)
      basic_df$p_value = sapply(meta_res_coeff_tables,function(x,y)x[y,"pval"],y=curr_rowname)
      basic_df$logFC = sapply(meta_res_coeff_tables,function(x,y)x[y,"estimate"],y=curr_rowname)
      basic_df$logFC_se = sapply(meta_res_coeff_tables,function(x,y)x[y,"se"],y=curr_rowname)
      basic_df$meta_reg_het_p = sapply(meta_res,function(x)x$QEp)
      basic_df$meta_reg_pvalue = sapply(meta_res,function(x)x$p_value)
      curr_timewise_df = rbind(curr_timewise_df,basic_df)
    }
  }
  
  # take the shared cols, merge meta-analysis results with the DEA results of the unique metabolites
  meta_reg_timewise_dea = merge_two_dea_dfs(curr_timewise_df,timewise[unique_metabs_inds,])

  ##############################################################################
  # make the training differential analysis results tables 
  ##############################################################################

  training_dea_without_meta = training
  unique_metabs_inds = !(training$metabolite_refmet %in% names(meta_res))
  # change the indices of some high het cases to TRUE here
  # this will make the algorithm add them
  tissue_metab_classes = metab_classes
  if(length(tissue_metab_classes[["high_het_metabs_targeted"]]) > 0){
    unique_metabs_inds[
      training$metabolite_refmet %in% tissue_metab_classes[["high_het_metabs_targeted"]] &
        training$is_targeted
    ] = T
  }
  if(length(tissue_metab_classes[["high_het_metabs_untargeted_only"]]) > 0){
    unique_metabs_inds[
      training$metabolite_refmet %in% tissue_metab_classes[["high_het_metabs_untargeted_only"]]
    ] = T
  }
  # finally, for the meta_reg obj keep the low het metabs only
  meta_res = meta_res[tissue_metab_classes[["meta_anal_metabs"]]]
  
  basic_df = c()
  if(length(meta_res)>0){
    basic_df = t(sapply(meta_res,function(x,y)unlist(x[y]),
                        y=union(basic_cols_to_keep,additional_training_cols)))
    basic_df = as.data.frame(basic_df,stringsAsFactors = F)
    basic_df$covariates = NA
    basic_df$meta_reg_het_p = sapply(meta_res,function(x)x$QEp)
    basic_df$p_value = sapply(meta_res,function(x)x$p_value)
    basic_df$original_ftest_ps = sapply(meta_res,function(x)x$original_ftest_ps)
  }
  
  meta_reg_training_dea = merge_two_dea_dfs(basic_df,training[unique_metabs_inds,])
    
  ##############################################################################
  # return results  
  ##############################################################################
  
  return(list(timewise_meta_regression = meta_reg_timewise_dea,
              training_meta_regression = meta_reg_training_dea,
              meta_regression_results = meta_res,
              meta_regression_models = metareg_run_table,
              metabolite_categories = metab_classes))
}


#' Concatenate data frames
#'
#' Concatenate (rbind) two data frames and expand the result to include all columns from both tables. 
#' 
#' @param x1 data frame
#' @param x2 data frame
#' 
#' @export 
#'
#' @examples 
#' merge_two_dea_dfs(
#'   data.frame(a=1:3,b=5:7),
#'   data.frame(a=1:3,b=5:7)
#' )
#' merge_two_dea_dfs(
#'   data.frame(a=1:3,b=5:7),
#'   data.frame(a=1:3,b=5:7,c=10:12)
#' )
#' merge_two_dea_dfs(
#'   data.frame(a=1:3,b=5:7,d=10:12),
#'   data.frame(a=1:3,b=5:7,c=10:12)
#' )
merge_two_dea_dfs = function(x1,x2){
  if(length(x1)==0 || nrow(x1)==0){return(x2)}
  if(length(x2)==0 || nrow(x2)==0){return(x1)}
  newrn = c(rownames(x1), rownames(x2))
  concat = data.table::rbindlist(list(x1,x2), fill=TRUE)
  concat = as.data.frame(concat)
  if(length(newrn) == length(unique(newrn))){
    rownames(concat) = newrn
  }
  return(concat)
}


#' Meta-regression for a metabolite 
#' 
#' Worker function wrapped by [metab_meta_regression()]. 
#' 
#' @param metabolite_refmet character, RefMet name of metabolite
#' @param timewise data frame with all timewise summary statistics
#' @param training optional data frame with all training summary statistics 
#'
#' @return a named list with meta-regression results, metabolite metadata, and 
#'  training summary statistics if \code{training} is not NULL. 
#' 
#' @details 
#' 
#' We try multiple models per analyte:  
#' \itemize{
#'   \item{Model 1:}{Two random effects factors if \code{platform} and \code{is_targeted} are not redundant. Default optimization.}
#'   \item{Model 2:}{Two random effects factors if \code{platform} and \code{is_targeted} are not redundant with alternative optimization.}
#'   \item{Model 3:}{\code{platform} and \code{is_targeted} are redundant. Use a single RE factor with default optimization. 
#'     Also, use this if \code{QMp} is NA, which is an indication of over-parameterization of the model.}
#'   \item{Model 4:}{\code{platform} and \code{is_targeted} are redundant. Use a single RE factor with alternative optimization.}
#'   \item{Model 5:}{If all previous analyses failed, use a simple fixed-effects approach.}
#' }
#' 
#' See <https://www.metafor-project.org/doku.php/tips:models_with_or_without_intercept> for some explanations.
#' See <https://www.publichealth.columbia.edu/research/population-health-methods/meta-regression> for general background,
#' including discussion on FE models for replicated experiments (as in our case). 
#'
#' "Use of a fixed effect meta-analysis model assumes all studies are estimating the same (common) treatment effect.":
#' - Riley et al. <https://www-bmj-com.stanford.idm.oclc.org/content/342/bmj.d549>
#' 
#' Hypothetically, if all studies had an infinite sample size,
#' there would be no differences due to chance and the differences in study estimates would completely disappear.
#'
#' Unlike the other assays, having greater sampling variance in
#' one of the sexes is not a problem for the model and we therefore do not need to split the sexes. 
#'
#' Most longitudinal meta-analyses include random effect terms per study. For examples see
#' <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0164898>, 
#' <https://metafor-project.org/doku.php/tips:two_stage_analysis>. 
#' 
metabolite_meta_regression = function(metabolite_refmet, timewise, training = NULL){
  
  if(!requireNamespace("metafor", quietly = TRUE)){
    stop(
      "Package 'metafor' must be installed to run 'metabolite_meta_regression()'.",
      call. = FALSE
    )
  }
  if (!requireNamespace("nloptr", quietly = TRUE)) {
    stop(
      "Package 'nloptr' must be installed to run 'metabolite_meta_regression()'.",
      call. = FALSE
    )
  }
  
  name = metabolite_refmet
  x = timewise

  # take the metabolite's dataset, add required variables
  x_subset = x[x$metabolite_refmet==name,]
  x_subset$analysis_group = factor(paste(x_subset$sex,x_subset$comparison_group,sep=","))
  x_subset$platform = paste(x_subset$site,x_subset$dataset,sep="_")
  x_subset$istargeted = x_subset$site %in% c("duke","emory","mayo")
  two_models_tb = table(x_subset$platform,x_subset$istargeted)
  use_two_re_factors = length(two_models_tb) > 2 &&
    (any(rowSums(two_models_tb>0) > 1) || any(colSums(two_models_tb>0) > 1))
  
  # stop if we cannot meta-analyze
  if(is.null(dim(table(x_subset$site,x_subset$dataset)))){
    stop(paste("Metabolite",name,
               "has either data from a single site or a single platform."))
  }
  
  res = NULL
  # try different solutions
  # model 1: two random effects factors if platform and is_targeted are not redundant
  # with default optimization
  if(use_two_re_factors){
    try({
      res = metafor::rma.mv(yi = logFC,
                   V = logFC_se^2,
                   mods = ~ 0 + analysis_group,
                   random = list(~ analysis_group|platform,
                                 ~analysis_group|istargeted),
                   data = x_subset,
                   slab = paste(site,sex,comparison_group))
      res$p_value = res$QMp
    }, silent = TRUE)
  }
  # model 2: two random effects factors if platform and is_targeted are not redundant
  # with alternative optimization
  if((is.null(res) || is.na(res$QMp)) && use_two_re_factors){
    try({
      res = metafor::rma.mv(yi = logFC,
                   V = logFC_se^2,
                   mods = ~ 0 + analysis_group,
                   random = list(~ analysis_group|platform,
                                 ~analysis_group|istargeted),
                   data = x_subset,
                   control=list(optimizer="nloptr", algorithm="NLOPT_LN_SBPLX", ftol_rel=1e-6),
                   slab = paste(site,sex,comparison_group))
      res$p_value = res$QMp
    }, silent = TRUE)
  }
  
  # model 3: is_targeted are redundant: use a single RE factor
  # with default optimization
  # also, use this if QMp is NA, which is an indication of over-parameterization
  # of the model
  if((is.null(res) && !use_two_re_factors) ||
     (!is.null(res) && is.na(res$QMp))){
    try({
      res = metafor::rma.mv(yi = logFC,
                   V = logFC_se^2,
                   mods = ~ 0 + analysis_group,
                   random = ~ analysis_group|platform,
                   data = x_subset,
                   slab = paste(site,sex,comparison_group))
      res$p_value = res$QMp
    }, silent = TRUE)
  }
  # model 4: is_targeted are redundant: use a single RE factor
  # with alternative optimization
  # also, use this if QMp is NA, which is an indication of over-parameterization
  # of the model
  if((is.null(res) && !use_two_re_factors) ||
     (!is.null(res) && is.na(res$QMp))){
    try({
      res = metafor::rma.mv(yi = logFC,
                   V = logFC_se^2,
                   mods = ~ 0 + analysis_group,
                   random = ~ analysis_group|platform,
                   data = x_subset,
                   control=list(optimizer="nloptr", algorithm="NLOPT_LN_SBPLX", ftol_rel=1e-6),
                   slab = paste(site,sex,comparison_group))
      res$p_value = res$QMp
    }, silent = TRUE)
  }
  # model 5: if all prev analyses failed, use a simple FE approach
  if(is.null(res) || is.na(res$QMp)){
    res = metafor::rma(yi = logFC,
              se = logFC_se,method = "FE",
              mods = ~ 0 + analysis_group, data = x_subset,
              slab = paste(site,sex,comparison_group))
    non_site_inds = which(grepl("analysis_group",rownames(res$b)))
    res$p_value = metafor::anova.rma(res,btt = non_site_inds)$QMp
  }
  
  # add extra information to the res object
  res$assay = "metab"
  res$tissue = x_subset$tissue[1]
  res$dataset = "meta-reg"
  # give QEp a friendly name
  res$het_p_value = res$QEp
  
  res$site = paste(unique(x_subset$platform),collapse=",")
  res$is_targeted = any(x_subset$is_targeted)
  res$sex = x_subset$sex[1]
  res$comparison_group = x_subset$comparison_group[1]
  res$metabolite = x_subset$metabolite_refmet[1]
  res$feature_ID = res$metabolite
  res$metabolite_refmet = x_subset$metabolite_refmet[1]
  res$nplatform = length(unique(x_subset$platform))
  
  res$cv = mean(x_subset$cv)
  res$comparison_average_intensity = NA
  res$reference_average_intensity = NA
  res$covariates = NA
  
  res$rawdata = x_subset
  
  training_dea = training
  if(!is.null(training_dea)){
    currtab = training_dea[which(training_dea$metabolite_refmet==name),]
    training_dea_ps = currtab[,"p_value"]
    res$original_ftest_ps = paste(training_dea_ps,collapse=",")
    mzs = currtab[,"mz"]
    res$mz = paste(mzs,collapse=",")
    rts = currtab[,"rt"]
    res$rt = paste(rts,collapse=",")
    mass = currtab[,"neutral_mass"]
    res$neutral_mass = paste(mass,collapse=",")
    ps_males = currtab[,"p_value_male"]
    res$p_value_male = paste(ps_males,collapse=",")
    ps_females = currtab[,"p_value_female"]
    res$p_value_female = paste(ps_females,collapse=",")
    f_males = currtab[,"fscore_male"]
    res$fscore_male = paste(f_males,collapse=",")
    f_females = currtab[,"fscore_female"]
    res$fscore_female = paste(f_females,collapse=",")
    sites = currtab[,"site"]
    datasets = currtab[,"dataset"]
    res$site = paste(paste(sites,datasets,sep="_"),collapse=",")
  }
  return(res)
}


#' Print forest plot
#' 
#' Make a forest plot to present multiple measurements and consensus effect
#' size from meta-analysis or meta-regression. Note that using [metab_meta_regression()] 
#' results as input generates a single forest plot while using [metab_meta_analysis()]
#' generates one forest plot per sex and training time point. 
#'
#' @param results named list returned by [metab_meta_analysis()] or [metab_meta_regression()] 
#' @param metabolite_refmet character, RefMet name of metabolite 
#'
#' @export
#' 
#' @seealso [metab_meta_analysis()], [metab_meta_regression()]
#'
#' @examples
#' # Get meta-analysis results in gastrocnemius
#' res = metab_meta_analysis("SKM-GN")
#' # Pick a feature
#' res$merged_res[res$merged_res$min_nominal_p > 0.01 & 
#'   res$merged_res$p_value < 0.001 & 
#'   res$merged_res$I2 < 30,]
#' forest_plot(res, metabolite_refmet="Hydroxyproline")
#' 
#' # Look at meta-regression results for a feature in the plasma
#' res = metab_meta_regression("PLASMA")
#' forest_plot(res, metabolite_refmet="CAR(18:2)")
#' 
#' res = metab_meta_regression("BAT")
#' forest_plot(res, metabolite_refmet="Aspartic acid")
forest_plot = function(results, metabolite_refmet){
  
  if(!requireNamespace("metafor", quietly = TRUE)){
    stop(
      "Package 'metafor' must be installed to run 'forest_plot()'.",
      call. = FALSE
    )
  }
  
  if("meta_analysis_res" %in% names(results)){
    type = "metaanalysis"
  }else{
    type = "metaregression"
  }
  
  if(type=="metaanalysis"){
    names = paste0(paste0(metabolite_refmet, ",", paste0(c(1,2,4,8),"w")), ",", rep(c("male","female"),each=4))
    for(name in names){
      
      splits = unname(unlist(strsplit(name, ",")))
      metabolite_refmet = splits[1]
      timepoint = splits[2]
      sex = splits[3]
      
      if(!name %in% names(results$meta_analysis_res)){
        warning(sprintf("Meta-analysis results for %s are not available in these results.", name))
        next
      }
      metap = unique(results$merged_res$p_value[results$merged_res$metabolite_refmet == metabolite_refmet &
                                                  results$merged_res$sex == sex &
                                                  results$merged_res$comparison_group == timepoint])
      metap = round(metap, digits=5)
      minp = unique(results$merged_res$min_nominal_p[results$merged_res$metabolite_refmet == metabolite_refmet &
                                                       results$merged_res$sex == sex &
                                                       results$merged_res$comparison_group == timepoint])
      minp = round(minp, digits=5)
      metafor::forest(results$meta_analysis_res[[name]],
                      main = sprintf("\n\n%s\nmeta-analysis p = %s\nmin p = %s", name, metap, minp),
                      col = "blue") 
    }
  }else if(type=="metaregression"){
    tissue = results$training_meta_regression$tissue[1]
    obj = results$meta_regression_results
    minrawp = sapply(obj, function(x){
      return(min(x$rawdata[,"p_value"]))
    })
    rawdsize = sapply(obj,function(x)nrow(x$rawdata))
    metaps = sapply(obj,function(x)x$p_value)
    meta_hetp = sapply(obj,function(x)x$QEp)
    currnames = strsplit(attr(obj[[metabolite_refmet]]$yi,"slab"),split = " ")
    currnames = t(sapply(currnames,function(x)x))
    curr_cvs = unique(obj[[metabolite_refmet]]$rawdata[,c("site","dataset","control_cv")])
    obj[[metabolite_refmet]]$rawdata$control_cv = round(obj[[metabolite_refmet]]$rawdata$control_cv,2)
    curr_labels  = paste0(obj[[metabolite_refmet]]$slab," (",obj[[metabolite_refmet]]$rawdata$control_cv,")")
    metafor::forest(obj[[metabolite_refmet]], 
                    slab = curr_labels,
                    addfit = F,
                    main=paste0( "\n\n\n\n",tissue,": ",metabolite_refmet,"\n",
                                 "meta-reg p=",format(round(obj[[metabolite_refmet]]$p_value,digits = 20),digits=2),"\n",
                                 "min p=", format(round(minrawp[metabolite_refmet],digits = 20),digits=2),", ",
                                 "meta-reg het p=",format(round(obj[[metabolite_refmet]]$QEp,digits = 20),digits=2)),
                    top = 8,
                    col = "blue",
                    cex = 0.9,
                    order = order(currnames[,2],currnames[,1],currnames[,3]))
  }else{
    stop("How did we get a type other than 'metaanalysis', 'metaregression'?")
  }
}
