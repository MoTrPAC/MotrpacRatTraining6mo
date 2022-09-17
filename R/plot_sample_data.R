# TODO
plot_feature_normalized_data = function(){}

# TODO
plot_feature_logfc = function(){}


#' Plot feature trajectories 
#' 
#' Plot group means of a set of features from normalized sample-level data.
#' 
#' @param features character vector of features to plot in the format 
#'   '[MotrpacRatTraining6moData::ASSAY_ABBREV];[MotrpacRatTraining6moData::TISSUE_ABBREV];feature_ID'
#' @param exclude_outliers bool, whether to remove sample outliers specified by [MotrpacRatTraining6moData::OUTLIERS].
#'   \code{TRUE} by default. 
#' @param center bool, whether to center the trajectories. \code{TRUE} by default.
#' @param scale bool, whether to scale the trajectories. \code{TRUE} by default.
#' @param title optional character, plot title
#' @param return_data bool, whether to return the normalized sample-level data 
#'   corresponding to \code{features} instead of a plot. \code{FALSE} by default. 
#' @param scratchdir character, local directory in which to download data from the web. 
#'   Current working directory by default. Only relevant if \code{features} includes ATAC or METHYL features. 
#' 
#' @return a [ggplot2::ggplot()] object if \code{return_data=FALSE} or a data frame otherwise
#' 
#' @export 
#' 
#' @examples 
#' # Pick largest cluster in gastrocnemius 
#' clust = extract_tissue_sets("SKM-GN", k=1, add_week8=FALSE)
#' # Extract features 
#' names(clust)
#' features = clust[["1w_F1_M1->2w_F1_M1->4w_F1_M1->8w_F1_M1"]]
#' plot_feature_trajectories(features)
plot_feature_trajectories = function(features, exclude_outliers=TRUE, center=TRUE, scale=TRUE, title=NULL, return_data=FALSE, scratchdir="."){
  
  if(!requireNamespace("viridis", quietly = TRUE) | !requireNamespace("viridisLite", quietly = TRUE)){
    stop(
      "Packages 'viridis' and 'viridisLite' must be installed to run 'plot_feature_trajectories()'.",
      call. = FALSE
    )
  }
  
  message("Identifying data sets...")
  features = unique(features)
  feature_dt = data.table::as.data.table(check_cluster_res_format(data.frame(feature=features, cluster="dummy")))
  datasets = unique(feature_dt[,.(ome, tissue)])
  omes = unique(datasets[,ome])
  
  # choose which epigenetic dataset to use 
  training_regulated_only = FALSE
  if("ATAC" %in% omes | "METHYL" %in% omes){
    epigen_input = feature_dt[ome %in% c("ATAC","METHYL"), feature]
    epigen_training_reg = MotrpacRatTraining6moData::TRAINING_REGULATED_FEATURES$feature[MotrpacRatTraining6moData::TRAINING_REGULATED_FEATURES$assay %in% c("METHYL", "ATAC")]
    if(all(epigen_input %in% epigen_training_reg)){
      training_regulated_only = TRUE
    }
  }
  
  message("Compiling sample-level data...")
  res = list()
  for(i in 1:nrow(datasets)){
    .ome = datasets[i, ome]
    .tissue = datasets[i, tissue]
    data = NULL
    # get sample-level data 
    if (.ome %in% c("ATAC","METHYL")){
      data = load_sample_data(.tissue, 
                              .ome, 
                              normalized=TRUE, 
                              training_regulated_only=training_regulated_only, 
                              exclude_outliers=exclude_outliers, 
                              scratchdir=scratchdir,
                              warnings=TRUE)
    }else{
      data = load_sample_data(.tissue, 
                              .ome, 
                              normalized=TRUE, 
                              training_regulated_only=FALSE, 
                              exclude_outliers=exclude_outliers, 
                              scratchdir=scratchdir,
                              warnings=TRUE)
    }
    if(is.null(data)) next
    # convert colnames to PID
    viallabel_cols = colnames(data)[grepl("^9", colnames(data))]
    if(length(viallabel_cols)>0){
      pids = viallabel_to_pid(viallabel_cols)
      stopifnot(length(pids) == length(viallabel_cols))
      stopifnot(length(pids) == length(unique(pids)))
      # rename columns 
      new_colnames = as.character(unname(pids[viallabel_cols]))
      colnames(data)[grepl("^[0-9]", colnames(data))] = new_colnames
    }
    # add new feature column 
    data = data.table::as.data.table(data)
    data[,new_feature := sprintf("%s;%s;%s", assay, tissue, feature_ID)]
    # select features 
    curr_feat = feature_dt[ome==.ome & tissue==.tissue, feature]
    data = data[feature %in% curr_feat | new_feature %in% curr_feat]
    data[,feature := ifelse(feature %in% curr_feat, feature, new_feature)]
    data[,new_feature := NULL]
    if(nrow(data)>0){
      # add to result
      res[[sprintf("%s_%s",.ome,.tissue)]] = data
    }else{
      warning(sprintf("No unfiltered features for %s %s.", .ome, .tissue))
    }
  }
  if(length(res)==0){
    warning(sprintf("No normalized data returned for datasets %s.",
                    paste(paste(datasets[,ome], datasets[,tissue], sep=";"), collapse=", ")))
    return()
  }
  
  sample_level_data = data.table::rbindlist(res, fill=TRUE)
  
  # check if features are present 
  if(!all(features %in% sample_level_data[,feature])){
    # what's missing?
    missing = features[!features %in% sample_level_data[,feature]]
    warning(sprintf("%s out of %s features were not found in the normalized sample-level data:\n%s", 
                    length(missing), 
                    length(features), 
                    paste(missing, collapse="\n")))
  }
  
  # handle duplicate row names 
  sample_level_data[,feature := paste0("feature", 1:nrow(sample_level_data))]
  
  # melt
  sample_cols = colnames(sample_level_data)[grepl("^[0-9]", colnames(sample_level_data))]
  melted_subset = data.table::melt(sample_level_data, 
                                   id.vars=c('feature'), 
                                   measure.vars=sample_cols, 
                                   variable.name='sample')
  melted_subset = melted_subset[!is.na(value)]
  melted_subset[,sample := as.integer(as.character(sample))]
  
  meta = data.table::as.data.table(MotrpacRatTraining6moData::PHENO)
  
  # merge by pid 
  meta = unique(meta[,.(group,sex,pid)])
  subset_meta = merge(melted_subset, meta, by.x='sample', by.y='pid')
  
  bygroup = subset_meta[,list(expr = mean(value, na.rm=T)),
                        by=.(sex, group, feature)]
  tmm_wide = data.frame(data.table::dcast(bygroup, feature ~ sex + group, value.var='expr'))
  rownames(tmm_wide) = tmm_wide$feature
  tmm_wide$feature = NULL
  
  if(center | scale){
    tmm_wide = as.data.frame(t(scale(t(tmm_wide), center=center, scale=scale)))
  }
  ylab = "Normalized value"
  
  tmm_wide = data.table::data.table(cbind(data.table(feature=rownames(tmm_wide)), tmm_wide))
  
  tmm_melt = data.table::melt(tmm_wide, id.vars="feature")
  tmm_melt[,sex := gsub('_.*','',variable)]
  tmm_melt[,group := gsub('.*_','',variable)]
  means = tmm_melt[,list(value = mean(value, na.rm=T)), by=c('sex','group')]
  
  # calculate distance to average
  m2 = data.table::data.table(merge(tmm_melt, means, by=c('sex','group'), suffixes = c('_feature','_mean')))
  distances = m2[,list(ss = sum((value_feature - value_mean)^2)), by=feature]
  tmm_melt = merge(tmm_melt, distances, by='feature')
  
  if(return_data){
    return(list(norm_melt=tmm_melt,
                centroids=means))
  }
  
  g = ggplot2::ggplot() +
    ggplot2::geom_line(data=tmm_melt, alpha=0.5, ggplot2::aes(x=group, y=value, colour=log(ss), group=feature)) +
    ggplot2::geom_line(data=means, ggplot2::aes(x=group, y=value, group=sex)) +
    ggplot2::facet_wrap(~sex) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_discrete(limits=c('control','1w','2w','fill','4w',rep('fill',3),'8w'),
                              breaks=c('control','1w','2w','4w','8w'),
                              labels=c('0','1','2','4','8')) +
    ggplot2::labs(x='Time trained (weeks)', y=ylab, title=title) +
    viridis::scale_color_viridis(direction = 1, option = "magma", guide='none') +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank())
  
  return(g)
}
