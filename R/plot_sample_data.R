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
#' @return data frame... TODO
#' 
#' @examples
#' print("TODO")
combine_sample_level_data = function(tissues=NULL, 
                                     assays=NULL,
                                     training_regulated_only=TRUE){
  
}


# TODO
#' Combine differential analysis results 
#' 
#' Combine differential analysis results from the specified tissues and omes(s)/assay(s). 
#' If no tissues or omes are specified, 
#' all differential analysis results are returned. 
#'
#' @param tissues optional character vector of tissue abbreviations, one of [MotrpacRatTraining6moData::TISSUE_ABBREV].
#' @param assays optional character vector of assay abbreviations, one of [MotrpacRatTraining6moData::ASSAY_ABBREV]
#' @param training_regulated_only bool, whether to return data only for training-regulated features
#' 
#' @export
#' 
#' @return data frame... TODO
#' 
#' @details For **all** epigenetic differential analysis results, see the 
#'   documentation for [MotrpacRatTraining6moData::ATAC_DA] and 
#'   [MotrpacRatTraining6moData::METHYL_DA], e.g., \code{?ATAC_DA}. 
#' 
combine_da_results = function(){}

# TODO
plot_feature_normalized_data = function(){}

# TODO
plot_feature_logfc = function(){}

# TODO
# formerly plot_group_mean_trajectories
plot_cluster_trajectory = function(){}
