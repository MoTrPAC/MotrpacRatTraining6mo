#' Metabolomics timewise meta-analysis
#' 
#' Perform meta-analysis for repeated measurements and return non-redundant results.
#' If there were no repeated measurements for a feature, the input summary statistics
#' are returned. **These results were not used in the manuscript.**
#'
#' @param tissue `r tissue()`
#' @param input r data frame, custom input. To see the expected format, look at 
#'   a table returned by [load_metabolomics_da()] with \code{type="timewise"}. 
#'
#' @return a named list, where \code{res$meta_analysis_res} is a list with one [metafor::rma.uni()] object per
#'   meta-analyzed feature, and \code{res$merged_res} is a data frame non-redundant measurements for all input features
#'   with the following columns: 
#' \describe{
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{dataset}}{`r dataset_metab()`}
#'   \item{\code{assay}}{`r assay()` ("METAB")}
#'   \item{\code{sex}}{`r sex()`}
#'   \item{\code{is_targeted}}{logical, is this a targeted platform?}
#'   \item{\code{site}}{character, Chemical Analysis Site (CAS) name}
#'   \item{\code{metabolite_refmet}}{character, RefMet name of metabolite}
#'   \item{\code{logFC}}{`r logFC()`}
#'   \item{\code{logFC_se}}{`r logFC_se()`}
#'   \item{\code{zscore}}{numeric, z-score (either original t-score or z-score from meta-analysis)}
#'   \item{\code{cv}}{double, feature coefficient of variation in the dataset}
#'   \item{\code{comparison_group}}{`r comparison_group()`}
#'   \item{\code{p_value}}{`r p_value_da()`}
#'   \item{\code{comparison_average_intensity}}{`r comparison_average_intensity()`}
#'   \item{\code{reference_average_intensity}}{`r reference_average_intensity()`}
#'   \item{\code{I2}}{double, heterogeneity coefficient. 
#'     Higher values mean the measurements between sites/platforms are more discordant.}
#'   \item{\code{QEp}}{double, heterogeneity p-value. A smaller p-value means
#'     the measurements between sites/platforms are more discordant.} 
#'   \item{\code{min_nominal_p}}{double, minimum nominal p-value for this sex
#'     and time point across potentially multiple measurements in the input} 
#' }
#' 
#' @export
#' 
#' @examples
#' # Perform meta-analysis on gastrocnemius timewise summary statistics
#' res = metab_meta_analysis("SKM-GN")
metab_meta_analysis = function(tissue, input=NULL){
  
  if(!requireNamespace("metafor", quietly = TRUE)){
    stop(
      "Package 'metafor' must be installed to perform meta-analysis.",
      call. = FALSE
    )
  }
  
  req_input_cols = c("feature_ID","assay","tissue","dataset","site","is_targeted","sex",
                     "comparison_group","metabolite_refmet",
                     "effect_size","effect_size_se",
                     "comparison_average_intensity","reference_average_intensity",
                     "p_value","cv")
  
  if(is.null(input)){
    if(tissue == "BLOOD"){
      stop("Metabolomics was not performed in whole blood. Did you mean tissue='PLASMA'?")
    }
    # load non-redundant timewise differential analysis results 
    timewise = load_metabolomics_da(tissue, type="timewise")
  }else{
    timewise = input
  }
  # check/change col names
  if("logFC" %in% colnames(timewise)){
    setnames(timewise, c("logFC"), c("effect_size"))
  }
  if("logFC_se" %in% colnames(timewise)){
    setnames(timewise, c("logFC_se"), c("effect_size_se"))
  }
  if("tscore" %in% colnames(timewise)){
    setnames(timewise, c("tscore"), c("zscore"))
  }

  if(!all(req_input_cols %in% colnames(timewise))){
    stop(sprintf("The following columns are expected in the input:\n %s", paste0(req_input_cols, collapse=", ")))
  }
  
  # run the timewise analyses
  cols_to_keep = c(req_input_cols, "I2","QEp","zscore")
  
  x = timewise
  # merge the meta-analysis and the unique results
  analysis_names = apply(x[,c("metabolite_refmet","comparison_group","sex")],
                         1,paste,collapse=",")
  unique_mets = names(which(table(analysis_names)==1))
  inds = analysis_names %in% unique_mets
  x_unique = x[inds,]
  rownames(x_unique) = analysis_names[inds]
  if(length(unique_mets)==nrow(x)){
    message("No redundant measurements. Returning original summary statistics.")
    return(list(meta_analysis_res = list(),
                merged_res = x_unique))
  }
  
  # keep track of min p-values per feature
  minp = data.table::as.data.table(x)[,list(min_nominal_p = min(p_value, na.rm=TRUE)), 
                                      by=.(metabolite_refmet, sex, comparison_group)]
  
  # perform meta-analysis 
  subset_res = get_repeated_met_data(x)
  metaanalysis_mets = unique(subset_res[[2]])
  names(metaanalysis_mets) = metaanalysis_mets
  meta_res = lapply(metaanalysis_mets,
                    get_meta_analysis_results,
                    x_subset = subset_res[[1]],
                    subset_names = subset_res[[2]])
  meta_res_list = meta_res
  message(paste("Computed",length(meta_res),"meta-analysis results."))

  # transform results to a summary data frame
  meta_res = t(sapply(meta_res,function(x,y)unlist(x[y]),y=cols_to_keep))
  meta_res = as.data.frame(meta_res)
  meta_res$I2 = as.numeric(as.character(meta_res$I2))
  meta_res$effect_size = as.numeric(as.character(meta_res$effect_size))
  meta_res$effect_size_se = as.numeric(as.character(meta_res$effect_size_se))
  meta_res$p_value = as.numeric(as.character(meta_res$p_value))
  meta_res$zscore = as.numeric(as.character(meta_res$zscore))
  meta_res$comparison_average_intensity =
    as.numeric(as.character(meta_res$comparison_average_intensity))
  meta_res$reference_average_intensity =
    as.numeric(as.character(meta_res$reference_average_intensity))
  meta_res$QEp = as.numeric(as.character(meta_res$QEp))
  meta_res$cv = as.numeric(as.character(meta_res$cv))
  meta_res$is_targeted = as.logical(as.character(meta_res$is_targeted))
  names(meta_res)[names(meta_res)=="pval"] = "p_value"
  for(j in 1:length(meta_res)){
    if(is.factor(meta_res[[j]])){
      meta_res[[j]] = as.character(meta_res[[j]])
    }
  }
  
  # take the shared cols between original non-redundant results and meta-analysis merged results
  shared_cols = intersect(names(x_unique),names(meta_res))
  merged_results = rbind(meta_res[,shared_cols],x_unique[,shared_cols])
  # add meta-res cols
  merged_results$I2 = NA
  merged_results$QEp = NA
  merged_results[rownames(meta_res),"I2"] = meta_res$I2
  merged_results[rownames(meta_res),"QEp"] = meta_res$QEp
  
  # change colnames back
  setnames(merged_results, c("effect_size","effect_size_se"), c("logFC","logFC_se"))
  
  # merge back to save min_nominal_p
  merged_results2 = merge(merged_results, minp, by=c("sex","comparison_group","metabolite_refmet"))
  stopifnot(nrow(merged_results2) == nrow(merged_results))
  
  return(list(meta_analysis_res = meta_res_list,
              merged_res = merged_results2))
}

get_repeated_met_data = function(x){
  
  if(!requireNamespace("metafor", quietly = TRUE)){
    stop(
      "Package 'metafor' must be installed to run 'get_repeated_met_data()'.",
      call. = FALSE
    )
  }
  
  analysis_names = apply(x[,c("metabolite_refmet","comparison_group","sex")],
                         1,paste,collapse=",")
  unique_mets = names(which(table(analysis_names)==1))
  inds = !(analysis_names %in% unique_mets)
  x_subset = x[inds,]
  subset_names = analysis_names[inds]
  return(list(x_subset,subset_names))
}


get_meta_analysis_results = function(name, x_subset, subset_names){
  x_subset = x_subset[subset_names==name,]
  labs = paste(x_subset$site, x_subset$dataset, sep="_")
  res = metafor::rma.uni(yi = x_subset$effect_size,
                sei = x_subset$effect_size_se,
                method = "FE",
                slab = labs)
  # add extra information to the res object
  res$assay = "metab"
  res$tissue = x_subset$tissue[1]
  res$dataset = "meta-analysis"
  res$site = paste(labs,collapse=",")
  res$is_targeted = any(x_subset$is_targeted)
  res$sex = x_subset$sex[1]
  res$comparison_group = x_subset$comparison_group[1]
  res$metabolite = x_subset$metabolite_refmet[1]
  res$feature_ID = res$metabolite
  res$metabolite_refmet = x_subset$metabolite_refmet[1]
  res$cv = mean(x_subset$cv)
  res$comparison_average_intensity = mean(x_subset$comparison_average_intensity)
  res$reference_average_intensity = mean(x_subset$reference_average_intensity)
  res$covariates = NA
  res$is_named = x_subset$is_named[1]
  res$effect_size = res$beta
  res$effect_size_se = res$se
  res$p_value = res$pval
  res$zscore = res$zval
  
  return(res)
}


#' Print forest plot
#' 
#' Make a forest plot to present multiple measurements and consensus effect
#' size from meta-analysis. 
#'
#' @param results named list returned by [metab_meta_analysis()] 
#' @param metabolite_refmet character, RefMet name of metabolite 
#' @param sex `r sex()`
#' @param timepoint `r comparison_group()`
#' @param name character, string in the format \code{metabolite_refmet,sex,timepoint}. 
#'   Use *either* this argument *or* \code{metabolite_refmet}, \code{sex}, and \code{timepoint}
#'   to specify a result. 
#'
#' @export
#' 
#' @seealso [metab_meta_analysis()] 
#'
#' @examples
#' # Get meta-analysis results for gastrocnemius
#' res = metab_meta_analysis("SKM-GN")
#' # Pick a feature
#' res$merged_res[res$merged_res$min_nominal_p > 0.01 & 
#'   res$merged_res$p_value < 0.001 & 
#'   res$merged_res$I2 < 30,]
#' forest_plot(res, metabolite_refmet="Hydroxyproline", sex="female", timepoint="8w")
#' forest_plot(res, name="Hydroxyproline,8w,female")
forest_plot = function(results, metabolite_refmet=NULL, sex=NULL, timepoint=NULL, name=NULL){
  
  if(!requireNamespace("metafor", quietly = TRUE)){
    stop(
      "Package 'metafor' must be installed to run 'forest_plot()'.",
      call. = FALSE
    )
  }
  
  if(is.null(name) & any(is.null(c(metabolite_refmet, sex, timepoint)))){
    stop("If 'name' is NULL, all other arguments must be specified.")
  }
  if(!is.null(name) & any(!is.null(c(metabolite_refmet, sex, timepoint)))){
    warning("'name' is specified along with other character arguments. Using 'name' to extract the result.")
  }
  
  if(is.null(name)){
    if(!tolower(sex) %in% c("male","female")){
      stop("'sex' must be one of 'male', 'female'.")
    }
    if(!tolower(timepoint) %in% paste0(c(1,2,4,8),"w")){
      stop(sprintf("'timepoint' must be one of %s.", paste0(paste0(c(1,2,4,8),"w"), collapse=", ")))
    }
    name = sprintf("%s,%s,%s", metabolite_refmet, tolower(timepoint), tolower(sex))
  }else{
    splits = unname(unlist(strsplit(name, ',')))
    metabolite_refmet = splits[1]
    timepoint = splits[2]
    sex = splits[3]
  }

  if(!name %in% names(results$meta_analysis_res)){
    stop(sprintf("Meta-analysis results for %s are not available in these results.", name))
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
