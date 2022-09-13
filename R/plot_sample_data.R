# TODO
#' Combine sample-level data
#' 
#' Combine data from the specified tissues and omes(s)/assay(s). 
#' If no tissues or omes are specified, all data is returned. 
#' 
#' @param tissues optional character vector of tissue abbreviations, one of [MotrpacRatTraining6moData::TISSUE_ABBREV].
#' @param assays optional character vector of assay abbreviations, one of [MotrpacRatTraining6moData::ASSAY_ABBREV]
#' @param training_regulated_only bool, whether to return data only for training-regulated features
#' 
#' @export
#' 
#' @return data frame with features in rows and Participant IDs (PIDs) in columns
#' 
#' @examples
#' #TODO
combine_sample_level_data = function(tissues=NULL, 
                                     assays=NULL,
                                     training_regulated_only=TRUE){
  return()
}


#' Combine differential analysis results 
#' 
#' Combine differential analysis results from the specified tissues and omes(s)/assay(s). 
#' If no tissues or omes are specified, all differential analysis results are returned. 
#'
#' @param tissues optional character vector of tissue abbreviations, 
#'   one of [MotrpacRatTraining6moData::TISSUE_ABBREV].
#' @param assays optional character vector of assay abbreviations, 
#'   one of [MotrpacRatTraining6moData::ASSAY_ABBREV]
#' @param metareg bool, whether to use the meta-regression results for METAB.
#'   Use the upstream differential analysis results in \code{FALSE}.
#'   \code{TRUE} by default. 
#' @param include_epigen bool, whether to include the full ATAC or METHYL 
#'   differential analysis results from Google Cloud Storage. 
#'   Only relevant if \code{assays} includes "ATAC" or "METHYL".
#'   \code{FALSE} by default. 
#' @param scratchdir character, local directory in which to download data from the web. 
#'   Current working directory by default. Only relevant if \code{assays} includes "ATAC" or "METHYL".
#' 
#' @export
#' 
#' @return data frame. Depending on the specified assays, some of these columns may not be included:
#' \describe{
#'   \item{\code{feature}}{`r feature()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{sex}}{`r sex()`}
#'   \item{\code{comparison_group}}{`r comparison_group()`}
#'   \item{\code{p_value}}{`r p_value_da()`}
#'   \item{\code{adj_p_value}}{`r adj_p_value_da()`}
#'   \item{\code{logFC}}{`r logFC()`}
#'   \item{\code{logFC_se}}{`r logFC_se()`}
#'   \item{\code{tscore}}{`r tscore()`}
#'   \item{\code{covariates}}{`r covariates()`}
#'   \item{\code{numNAs}}{`r numNAs()`}
#'   \item{\code{comparison_average_intensity}}{`r comparison_average_intensity()`}
#'   \item{\code{reference_average_intensity}}{`r reference_average_intensity()`}
#'   \item{\code{selection_fdr}}{`r selection_fdr()`}
#'   \item{\code{dataset}}{character, immune panel, metabolomics platform, or ATAC-seq dataset name}
#'   \item{\code{site}}{character, Chemical Analysis Site (CAS) name. METAB only}
#'   \item{\code{is_targeted}}{logical, is this a targeted platform? METAB only}
#'   \item{\code{metabolite_refmet}}{character, RefMet name of metabolite. METAB only}
#'   \item{\code{cv}}{double, feature coefficient of variation in the dataset. METAB only}
#'   \item{\code{metabolite}}{character, name of metabolite as appears in the CAS's data. METAB only}
#'   \item{\code{control_cv}}{double, feature coefficient of variation in the dataset. METAB only}
#'   \item{\code{mz}}{double, mass over charge. METAB only}
#'   \item{\code{rt}}{double, retention time. METAB only}
#'   \item{\code{neutral_mass}}{double, neutral mass. METAB only}
#'   \item{\code{meta_reg_het_p}}{`r meta_reg_het_p()` METAB only}
#'   \item{\code{meta_reg_pvalue}}{`r meta_reg_pvalue()` METAB only}
#'   \item{\code{shrunk_logFC}}{double, log fold-change with shrinkage applied}
#'   \item{\code{shrunk_logFC_se}}{double, standard error of the shrunken log fold-change}
#'   \item{\code{zscore}}{`r zscore()`}
#'   \item{\code{removed_samples}}{`r removed_samples()`}
#'   \item{\code{comparison_average_intensity_se}}{`r comparison_average_intensity_se()`}
#'   \item{\code{reference_average_intensity_se}}{`r comparison_average_intensity_se()`} 
#'   \item{\code{Chr}}{integer, chromosome. METHYL only}
#'   \item{\code{Locus}}{character, base pair range of feature. METHYL only}
#'   \item{\code{EntrezID}}{character, Entrez ID of closest gene. METHYL only}
#'   \item{\code{Symbol}}{character, gene symbol of closest gene. METHYL only}
#' }
#' 
#' @examples
#' # Return all non-epigenetic differential analysis results, 
#' # including meta-regression results for metabolomics
#' combine_da_results()
#' 
#' # Return all global proteomics differential analysis results
#' combine_da_results(assays="PROT")
#' 
#' \dontrun{
#' # Return METHYL and ATAC differential analysis results for gastrocnemius 
#' combine_da_results(tissues="SKM-GN", 
#'                    assays=c("ATAC","METHYL"),
#'                    include_epigen=TRUE)
#' }
combine_da_results = function(tissues = MotrpacRatTraining6moData::TISSUE_ABBREV, 
                              assays = MotrpacRatTraining6moData::ASSAY_ABBREV, 
                              metareg = TRUE,
                              include_epigen = FALSE,
                              scratchdir = "."){
  
  if( ("ATAC" %in% assays | "METHYL" %in% assays) & !include_epigen){
    warning("'include_epigen' is FALSE. Excluding ATAC and METHYL results.")
    assays = assays[!assays %in% c("ATAC","METHYL")]
  }
  
  available_data = list_available_data("MotrpacRatTraining6moData")
  # add available epigen data
  possible_tissues = c(
    'BAT',
    'HEART',
    'HIPPOC',
    'KIDNEY',
    'LIVER',
    'LUNG',
    'SKMGN',
    'WATSC'
  )
  available_data = c(available_data, paste0(rep(c("ATAC", "METHYL"), each=8), "_", rep(possible_tissues,2), "_DA"))
  
  reslist = list()
  i = 1
  for(a in assays){
    for(t in tissues){
      if(a == "METAB" & metareg){
        object_name = sprintf("%s_%s_DA_METAREG", a, gsub("-","",t))
      }else{
        object_name = sprintf("%s_%s_DA", a, gsub("-","",t))
      }
      # check if object exists
      if(object_name %in% available_data){
        message(object_name)
        if (a %in% c("ATAC","METHYL")){
          # load from URL
          data = get_rdata_from_url(tissue=t, assay=a, suffix="DA", scratchdir=scratchdir)
        }else{
          data = fetch_object(object_name)
        }
        if ("removed_samples" %in% colnames(data)){
          data$removed_samples = as.character(data$removed_samples)
        }
        reslist[[i]] = data
        i = i+1
      }
    }
  }
  if(length(reslist) == 0){
    warning("No results returned.")
    return()
  }
  data = data.table::rbindlist(reslist, fill=TRUE)
  return(data)
}


# TODO
plot_feature_normalized_data = function(){}

# TODO
plot_feature_logfc = function(){}

# TODO
# formerly plot_group_mean_trajectories
plot_cluster_trajectory = function(){}
