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
#'   a table returned by [load_metabolomics_da()] with \code{type="timewise"}. 
#' @param training_input r data frame, custom input. To see the expected format, look at 
#'   a table returned by [load_metabolomics_da()] with \code{type="training"}. 
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
#' @seealso [load_metabolomics_da()], [metabolite_meta_regression()], [MotrpacRatTraining6moData::METAB_DA_METAREG]
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
    timewise = load_metabolomics_da(tissue, type="timewise")
  }else{
    timewise = timewise_input
  }
  
  # load redundant training differential analysis results 
  if(is.null(training_input)){
    training = load_metabolomics_da(tissue, type="training")
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
  
  basic_cols_to_keep = c("feature_ID","assay","tissue",'tissue_abbreviation',
                         "dataset","site","is_targeted","metabolite_refmet","cv","covariates","control_cv")
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



# all_meta_reg_results = meta_reg_training_dea[meta_reg_training_dea$dataset == "meta-reg",]
# 
# # compare p-values with and without the meta-regression
# get_char_minp<-function(x){
#   ps = as.numeric(strsplit(x,split=",")[[1]])
#   return(min(ps))
# }
# all_meta_reg_ps = pmax(1e-20,all_meta_reg_results$p_value)
# raw_min_p = unname(sapply(all_meta_reg_results$original_ftest_ps,get_char_minp))
# plot(
#   x = -log10(all_meta_reg_ps),
#   y = -log10(raw_min_p),
#   cex = 0.8,
#   col="gray",pch=20,
#   xlab = "-log 10 meta-regression p-value",ylab = "-log 10 min p-value in raw data"
# )
# abline(0,1)
# 
# print("Summary of meta-regression metabolites with training dea p-val < 0.001:")
# print(table(all_meta_reg_ps < 0.001))
# print("Compare the results above to the raw minimal p-value using the same 0.001 threshold:")
# table(all_meta_reg_ps < 0.001,raw_min_p < 0.001 )
# 
# # examine some extreme cases
# extreme_case = which(all_meta_reg_results$p_value < 1e-50)
# print("Here are some meta-regression results with extremely low p-value (< 1e-50:")
# all_meta_reg_results[extreme_case,]
# 
# 
# 
# ## Print the results to the bucket:
# 
# 
# ####################################################
# # Print the meta-regression results into the mawg bucket
# ####################################################
# tissues = unique(meta_reg_timewise_dea$tissue)
# mawg_bucket = "gs://mawg-data/pass1b-06/metabolomics/dea/data-freeze/meta-regression/"
# 
# suffix = "_metab-meta-reg_timewise-dea_20211014.txt"
# for(tissue in tissues){
#   x = meta_reg_timewise_dea[meta_reg_timewise_dea$tissue == tissue,]
#   if(is.null(x) || nrow(x)==0){next}
#   if(grepl("t54-hypo",tissue)){x$tissue = "t54-hypothalamus"}
#   local_fname = paste0(local_data_dir, "pass1b-06_",tissue,suffix)
#   write.table(x,file=local_fname,sep="\t",quote=F,
#               row.names = F,col.names = T)
#   cmd = paste(gsutil_cmd,"cp",local_fname,mawg_bucket)
#   system(cmd)
# }
# 
# suffix = "_metab-meta-reg_training-dea_20211014.txt"
# for(tissue in tissues){
#   x = meta_reg_training_dea[meta_reg_training_dea$tissue == tissue,]
#   if(is.null(x) || nrow(x)==0){next}
#   if(grepl("t54-hypo",tissue)){x$tissue = "t54-hypothalamus"}
#   local_fname = paste0(local_data_dir, "pass1b-06_",tissue,suffix)
#   write.table(x,file=local_fname,sep="\t",quote=F,
#               row.names = F,col.names = T)
#   cmd = paste(gsutil_cmd,"cp",local_fname,mawg_bucket)
#   system(cmd)
# }
# 
# 
# 
# ## Print the forest plots into pdfs
# 
# 
# min_rawp_sig<-function(x){return(min(x$rawdata[,"p_value"]))}
# report_hetp_thr = het_p_threshold
# report_pthr = 0.001
# curr_stats = c()
# 
# local_fname = "pass1b-06_metab_meta-reg-significant-high-het.pdf"
# pdf(local_fname)
# for(tissue in names(tissue2meta_res_objects)){
#   obj = tissue2meta_res_objects[[tissue]]
#   minrawp = sapply(obj,min_rawp_sig)
#   rawdsize = sapply(obj,function(x)nrow(x$rawdata))
#   metaps = sapply(obj,function(x)x$p_value)
#   meta_hetp = sapply(obj,function(x)x$QEp)
#   high_het_cases =  metaps < report_pthr & meta_hetp < report_hetp_thr
#   high_het_cases = names(which(high_het_cases))
#   low_het_cases =  names(which(meta_hetp > report_hetp_thr & metaps < report_pthr))
#   insignificant_res = metaps > report_pthr & minrawp > report_pthr
#   curr_stats = rbind(curr_stats,
#                      c(length(low_het_cases),length(high_het_cases),
#                        sum(is.na(insignificant_res)) + sum(insignificant_res,na.rm=T)))
#   rownames(curr_stats)[nrow(curr_stats)] = tissue
#   if(length(high_het_cases)==0){next}
#   for(currexample in high_het_cases){
#     currnames = strsplit(attr(obj[[currexample]]$yi,"slab"),split = " ")
#     currnames = t(sapply(currnames,function(x)x))
#     curr_cvs = unique(obj[[currexample]]$rawdata[,c("site","dataset","control_cv")])
#     obj[[currexample]]$rawdata$control_cv = round(obj[[currexample]]$rawdata$control_cv,2)
#     curr_labels  = paste0(obj[[currexample]]$slab," (",obj[[currexample]]$rawdata$control_cv,")")
#     forest(obj[[currexample]], slab = curr_labels,addfit = F,
#            main=paste0( "\n\n\n\n",tissue,": ",currexample,"\n",
#                         "meta-reg p=",format(round(obj[[currexample]]$p_value,digits = 20),digits=2),"\n",
#                         "min p=", format(round(minrawp[currexample],digits = 20),digits=2),", ",
#                         "meta-reg het p=",format(round(obj[[currexample]]$QEp,digits = 20),digits=2)
#                         #"I2 = ", round(metai[currexample],digits = 2)
#            ),top = 8,col = "blue",cex = 0.9,
#            order = order(currnames[,2],currnames[,1],currnames[,3])
#     )
#   }
# }
# dev.off()
# cmd = paste(gsutil_cmd,"cp",local_fname,mawg_bucket)
# system(cmd)
# 
# #curr_stats
# 
# local_fname = "pass1b-06_metab_meta-reg-significant-low-het-lowmetap.pdf"
# pdf(local_fname)
# for(tissue in names(tissue2meta_res_objects)){
#   obj = tissue2meta_res_objects[[tissue]]
#   minrawp = sapply(obj,min_rawp_sig)
#   rawdsize = sapply(obj,function(x)nrow(x$rawdata))
#   metaps = sapply(obj,function(x)x$p_value)
#   metai = sapply(obj,function(x)x$I2)
#   meta_hetp = sapply(obj,function(x)x$QEp)
#   low_het_cases =  names(which(meta_hetp > report_hetp_thr & metaps < report_pthr))
#   for(currexample in low_het_cases){
#     currnames = strsplit(attr(obj[[currexample]]$yi,"slab"),split = " ")
#     currnames = t(sapply(currnames,function(x)x))
#     curr_cvs = unique(obj[[currexample]]$rawdata[,c("site","dataset","control_cv")])
#     obj[[currexample]]$rawdata$control_cv = round(obj[[currexample]]$rawdata$control_cv,2)
#     curr_labels  = paste0(obj[[currexample]]$slab," (",obj[[currexample]]$rawdata$control_cv,")")
#     forest(obj[[currexample]], slab = curr_labels,addfit = F,
#            main=paste0( "\n\n\n\n",tissue,": ",currexample,"\n",
#                         "meta-reg p=",format(round(obj[[currexample]]$p_value,digits = 20),digits=2),"\n",
#                         "min p=", format(round(minrawp[currexample],digits = 20),digits=2),", ",
#                         "meta-reg het p=",format(round(obj[[currexample]]$QEp,digits = 20),digits=2)
#            ),top = 8,col = "blue",cex = 0.9,
#            order = order(currnames[,2],currnames[,1],currnames[,3])
#     )
#   }
# }
# dev.off()
# cmd = paste(gsutil_cmd,"cp",local_fname,mawg_bucket)
# system(cmd)
# 
# local_fname = "pass1b-06_metab_meta-reg-insignificant-minrawp-significant.pdf"
# pdf(local_fname)
# for(tissue in names(tissue2meta_res_objects)){
#   obj = tissue2meta_res_objects[[tissue]]
#   minrawp = sapply(obj,min_rawp_sig)
#   rawdsize = sapply(obj,function(x)nrow(x$rawdata))
#   metaps = sapply(obj,function(x)x$p_value)
#   metai = sapply(obj,function(x)x$I2)
#   meta_hetp = sapply(obj,function(x)x$QEp)
#   lost_results =  minrawp < report_pthr & metaps > report_pthr
#   lost_results = names(which(lost_results))
#   for(currexample in lost_results){
#     currnames = strsplit(attr(obj[[currexample]]$yi,"slab"),split = " ")
#     currnames = t(sapply(currnames,function(x)x))
#     curr_cvs = unique(obj[[currexample]]$rawdata[,c("site","dataset","control_cv")])
#     obj[[currexample]]$rawdata$control_cv = round(obj[[currexample]]$rawdata$control_cv,2)
#     curr_labels  = paste0(obj[[currexample]]$slab," (",obj[[currexample]]$rawdata$control_cv,")")
#     forest(obj[[currexample]], slab = curr_labels,addfit = F,
#            main=paste0( "\n\n\n\n",tissue,": ",currexample,"\n",
#                         "meta-reg p=",format(round(obj[[currexample]]$p_value,digits = 20),digits=2),"\n",
#                         "min p=", format(round(minrawp[currexample],digits = 20),digits=2),", ",
#                         "meta-reg het p=",format(round(obj[[currexample]]$QEp,digits = 20),digits=2)
#            ),top = 8,col = "blue",cex = 0.9,
#            order = order(currnames[,2],currnames[,1],currnames[,3])
#     )
#   }
# }
# dev.off()
# cmd = paste(gsutil_cmd,"cp",local_fname,mawg_bucket)
# system(cmd)
# 
# local_fname = "pass1b-06_metab_meta-reg-gained_results.pdf"
# pdf(local_fname)
# for(tissue in names(tissue2meta_res_objects)){
#   obj = tissue2meta_res_objects[[tissue]]
#   minrawp = sapply(obj,min_rawp_sig)
#   rawdsize = sapply(obj,function(x)nrow(x$rawdata))
#   metaps = sapply(obj,function(x)x$p_value)
#   metai = sapply(obj,function(x)x$I2)
#   meta_hetp = sapply(obj,function(x)x$QEp)
#   gained_results =  minrawp > report_pthr & metaps < report_pthr
#   gained_results = names(which(gained_results))
#   for(currexample in gained_results){
#     currnames = strsplit(attr(obj[[currexample]]$yi,"slab"),split = " ")
#     currnames = t(sapply(currnames,function(x)x))
#     curr_cvs = unique(obj[[currexample]]$rawdata[,c("site","dataset","control_cv")])
#     obj[[currexample]]$rawdata$control_cv = round(obj[[currexample]]$rawdata$control_cv,2)
#     curr_labels  = paste0(obj[[currexample]]$slab," (",obj[[currexample]]$rawdata$control_cv,")")
#     forest(obj[[currexample]], slab = curr_labels,addfit = F,
#            main=paste0( "\n\n\n\n",tissue,": ",currexample,"\n",
#                         "meta-reg p=",format(round(obj[[currexample]]$p_value,digits = 20),digits=2),"\n",
#                         "min p=", format(round(minrawp[currexample],digits = 20),digits=2),", ",
#                         "meta-reg het p=",format(round(obj[[currexample]]$QEp,digits = 20),digits=2)
#            ),top = 8,col = "blue",cex = 0.9,
#            order = order(currnames[,2],currnames[,1],currnames[,3])
#     )
#   }
# }
# dev.off()
# cmd = paste(gsutil_cmd,"cp",local_fname,mawg_bucket)
# system(cmd)
# 
# 
# ## Show a few forest plots
# 
# 
# # add this code here once because we sometimes want to rerun the flow without running
# # the code above that updates the buckets.
# min_rawp_sig<-function(x){return(min(x$rawdata[,"p_value"]))}
# 
# # An example for the presentation
# currexample = "CAR(18:2)"
# tissue = "t31-plasma"
# obj = tissue2meta_res_objects[[tissue]]
# minrawp = sapply(obj,min_rawp_sig)
# rawdsize = sapply(obj,function(x)nrow(x$rawdata))
# metaps = sapply(obj,function(x)x$p_value)
# meta_hetp = sapply(obj,function(x)x$QEp)
# currnames = strsplit(attr(obj[[currexample]]$yi,"slab"),split = " ")
# currnames = t(sapply(currnames,function(x)x))
# curr_cvs = unique(obj[[currexample]]$rawdata[,c("site","dataset","control_cv")])
# obj[[currexample]]$rawdata$control_cv = round(obj[[currexample]]$rawdata$control_cv,2)
# curr_labels  = paste0(obj[[currexample]]$slab," (",obj[[currexample]]$rawdata$control_cv,")")
# forest(obj[[currexample]], slab = curr_labels,addfit = F,
#        main=paste0( "\n\n\n\n",tissue,": ",currexample,"\n",
#                     "meta-reg p=",format(round(obj[[currexample]]$p_value,digits = 20),digits=2),"\n",
#                     "min p=", format(round(minrawp[currexample],digits = 20),digits=2),", ",
#                     "meta-reg het p=",format(round(obj[[currexample]]$QEp,digits = 20),digits=2)
#        ),top = 8,col = "blue",cex = 0.9,
#        order = order(currnames[,2],currnames[,1],currnames[,3])
# )
# 
# currexample = "Aspartic acid"
# tissue = "t69-brown-adipose"
# obj = tissue2meta_res_objects[[tissue]]
# minrawp = sapply(obj,min_rawp_sig)
# rawdsize = sapply(obj,function(x)nrow(x$rawdata))
# metaps = sapply(obj,function(x)x$p_value)
# meta_hetp = sapply(obj,function(x)x$QEp)
# currnames = strsplit(attr(obj[[currexample]]$yi,"slab"),split = " ")
# currnames = t(sapply(currnames,function(x)x))
# curr_cvs = unique(obj[[currexample]]$rawdata[,c("site","dataset","control_cv")])
# obj[[currexample]]$rawdata$control_cv = round(obj[[currexample]]$rawdata$control_cv,2)
# curr_labels  = paste0(obj[[currexample]]$slab," (",obj[[currexample]]$rawdata$control_cv,")")
# forest(obj[[currexample]], slab = curr_labels,addfit = F,
#        main=paste0( "\n\n\n\n",tissue,": ",currexample,"\n",
#                     "meta-reg p=",format(round(obj[[currexample]]$p_value,digits = 20),digits=2),"\n",
#                     "min p=", format(round(minrawp[currexample],digits = 20),digits=2),", ",
#                     "meta-reg het p=",format(round(obj[[currexample]]$QEp,digits = 20),digits=2)
#        ),top = 8,col = "blue",cex = 0.9,
#        order = order(currnames[,2],currnames[,1],currnames[,3])
# )
# 
