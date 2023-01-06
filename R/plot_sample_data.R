#' Plot sample-level data for a feature
#' 
#' Plot normalized sample-level data for a single feature. 
#' Points are mean values across samples in each group, 
#' and error bars indicate standard deviations. 
#'
#' @param assay NULL or `r assay()`
#' @param tissue NULL or `r tissue()`
#' @param feature_ID NULL or `r feature_ID()`
#' @param feature NULL or `r feature()`. If NULL, \code{assay}, \code{tissue}, and 
#'   \code{feature_ID} must all be specified. 
#' @param title character, plot title. By default, the plot ID is \code{feature}. 
#'   If \code{add_gene_symbol = TRUE}, the gene symbol is also added to the plot title.
#' @param add_gene_symbol bool, whether to add corresponding gene symbol to 
#'   plot title. Default: FALSE
#' @param facet_by_sex bool, whether to facet the plot by sex. If \code{TRUE},
#'   lines are colored by tissue. If \code{FALSE}, lines are colored by sex. Default: FALSE
#' @param scale_x_by_time bool, whether to scale the x-axis by time. If \code{FALSE},
#'   space the time points (0w, 1w, 2w, 4w, 8w) evenly. Default: TRUE
#' @param return_data bool, whether to return data instead of plot. Default: FALSE
#' @param exclude_outliers bool, whether to exclude data from sample outliers. Default: TRUE
#'   (see [MotrpacRatTraining6moData::OUTLIERS])
#' @param add_adj_p bool, whether to include the training adjusted p-value (AKA selection FDR)
#'   in the plot subtitle. Default: TRUE
#' @param ... additional arguments passed to [load_sample_data()]
#'
#' @return a [ggplot2::ggplot()] object or a data frame if \code{return_data = TRUE}
#'    or NULL if the data cannot be found
#' @export
#'
#' @examples
#' # Multiple ways of plotting the same data are shown for each example below
#' 
#' # Plot a differential feature 
#' plot_feature_normalized_data(feature = "ACETYL;HEART;NP_001003673.1_K477k",
#'                              add_gene_symbol = TRUE)
#' plot_feature_normalized_data(assay = "ACETYL",
#'                              tissue = "HEART",
#'                              feature_ID = "NP_001003673.1_K477k",
#'                              add_gene_symbol = TRUE,
#'                              scale_x_by_time = FALSE)
#' plot_feature_normalized_data(assay = "ACETYL",
#'                              tissue = "HEART",
#'                              feature_ID = "NP_001003673.1_K477k",
#'                              add_gene_symbol = TRUE,
#'                              facet_by_sex = TRUE)
#' 
#' # Plot a redundant differential feature
#' plot_feature_normalized_data(feature = "IMMUNO;PLASMA;BDNF",
#'                              add_gene_symbol = FALSE)
#' plot_feature_normalized_data(assay = "IMMUNO",
#'                              tissue = "PLASMA",
#'                              feature_ID = "BDNF",
#'                              add_gene_symbol = FALSE)
#' plot_feature_normalized_data(assay = "IMMUNO",
#'                              tissue = "PLASMA",
#'                              feature_ID = "BDNF",
#'                              add_gene_symbol = FALSE,
#'                              facet_by_sex = TRUE)
#'                              
#' # Plot one measurement of a redundant feature
#' plot_feature_normalized_data(feature = "IMMUNO;PLASMA;rat-myokine:BDNF",
#'                              add_gene_symbol = FALSE)
#' plot_feature_normalized_data(assay = "IMMUNO",
#'                              tissue = "PLASMA",
#'                              feature_ID = "rat-myokine:BDNF",
#'                              add_gene_symbol = FALSE)
#' plot_feature_normalized_data(assay = "IMMUNO",
#'                              tissue = "PLASMA",
#'                              feature_ID = "rat-myokine:BDNF",
#'                              add_gene_symbol = FALSE,
#'                              facet_by_sex = TRUE)
#'                              
#' # Plot a non-differential feature
#' plot_feature_normalized_data(feature = "PROT;SKM-GN;YP_665629.1",
#'                              add_gene_symbol = TRUE)
#' plot_feature_normalized_data(assay = "PROT",
#'                              tissue = "SKM-GN",
#'                              feature_ID = "YP_665629.1",
#'                              add_gene_symbol = TRUE)
#' plot_feature_normalized_data(assay = "PROT",
#'                              tissue = "SKM-GN",
#'                              feature_ID = "YP_665629.1",
#'                              add_gene_symbol = TRUE,
#'                              facet_by_sex = TRUE)
#'                              
#' # Plot a merged feature from meta-regression
#' plot_feature_normalized_data(assay = "METAB",
#'                              tissue = "PLASMA",
#'                              feature_ID = "Glucose",
#'                              facet_by_sex = TRUE)
#' plot_feature_normalized_data(assay = "METAB",
#'                              tissue = "PLASMA",
#'                              feature_ID = "glucose",
#'                              scale_x_by_time = FALSE)
#' plot_feature_normalized_data(assay = "METAB",
#'                              tissue = "PLASMA",
#'                              feature_ID = "glucose",
#'                              scale_x_by_time = FALSE,
#'                              add_adj_p = TRUE)
#'                              
plot_feature_normalized_data = function(assay = NULL,
                                        tissue = NULL, 
                                        feature_ID = NULL,
                                        feature = NULL, 
                                        title = NULL, 
                                        add_gene_symbol = FALSE, 
                                        facet_by_sex = FALSE, 
                                        scale_x_by_time = TRUE, 
                                        return_data = FALSE,
                                        exclude_outliers = TRUE,
                                        add_adj_p = FALSE,
                                        ...){
  
  curr_feature = feature 
  if(is.null(curr_feature)){
    if(any(is.null(c(assay, tissue, feature_ID)))){
      stop("If 'feature' is not specified, 'assay', 'tissue', and 'feature_ID' must all be specified.")
    }
    curr_feature = sprintf("%s;%s;%s", assay, tissue, feature_ID)
  }
  if(is.null(tissue)){
    splits = unname(unlist(strsplit(curr_feature, ";")))
    assay = splits[1]
    tissue = splits[2]
    feature_ID = splits[3]
  }
  # rename to avoid conflict with data.table columns
  FEATURE_ID = feature_ID
  ASSAY = assay 
  TISSUE = tissue 
  FEATURE = curr_feature 
  redundant_feature = FEATURE
  
  # is this a differential feature?
  differential = TRUE
  training_reg = data.table::as.data.table(TRAINING_REGULATED_FEATURES)
  if(!FEATURE %in% training_reg[,feature]){
    differential = FALSE
    # check if it is a repeated feature - TRAINING_REGULATED_FEATURES uses non-redundant feature
    if(ASSAY %in% c("METAB", "IMMUNO")){
      if(FEATURE %in% MotrpacRatTraining6moData::REPEATED_FEATURES$feature){
        FEATURE = MotrpacRatTraining6moData::REPEATED_FEATURES$new_feature[MotrpacRatTraining6moData::REPEATED_FEATURES$feature == FEATURE]
        differential = TRUE
        # now handle FEATURE_ID. get redundant one
        FEATURE_ID = MotrpacRatTraining6moData::REPEATED_FEATURES$feature_ID[MotrpacRatTraining6moData::REPEATED_FEATURES$new_feature == FEATURE[1]]
      }else{
        if(ASSAY == "METAB"){
          # check for refmet
          new_feature_id = unique(MotrpacRatTraining6moData::METAB_FEATURE_ID_MAP$feature_ID_metareg[MotrpacRatTraining6moData::METAB_FEATURE_ID_MAP$metabolite_name == FEATURE_ID])
          new_feature_id = unique(new_feature_id[new_feature_id %in% training_reg[assay == ASSAY & tissue == TISSUE,feature_ID]])
          FEATURE = unique(training_reg[feature_ID == new_feature_id & tissue == TISSUE & assay == ASSAY, feature])
          if(length(FEATURE) > 0){
            differential = TRUE
          }
        }
      }
    }
  }
  
  if(differential){
    if(exclude_outliers){
      all_sample_level_data = data.table::as.data.table(MotrpacRatTraining6moData::TRAINING_REGULATED_NORM_DATA_NO_OUTLIERS)
    }else{
      all_sample_level_data = data.table::as.data.table(MotrpacRatTraining6moData::TRAINING_REGULATED_NORM_DATA)
    }
    # subset data
    sample_level_data = all_sample_level_data[feature %in% FEATURE]
  }else{
    message(sprintf("'%s' is not a training-regulated feature. Looking in all sample-level data.", FEATURE))
    all_sample_level_data =  data.table::as.data.table(load_sample_data(TISSUE, ASSAY, exclude_outliers = exclude_outliers, ...))
    if(!FEATURE_ID %in% all_sample_level_data[,feature_ID]){
      warning(sprintf("'%s' not found in the %s %s sample-level data.", FEATURE_ID, ASSAY, TISSUE))
      return()
    }
    # subset data
    sample_level_data = all_sample_level_data[feature_ID == FEATURE_ID]
    sample_level_data[is.na(feature), feature := FEATURE]
  }
  
  multiple_measurements = FALSE
  if(nrow(sample_level_data) > 1){
    warning(sprintf("Multiple features correspond to '%s'. Plotting them together.", redundant_feature))
    # make feature non-redundant
    sample_level_data[,feature := dataset]
    multiple_measurements = TRUE
  }
  
  if(add_gene_symbol){
    if(ASSAY %in% c("METHYL","ATAC") & !differential){
      feature_to_gene = data.table::data.table(MotrpacRatTraining6moData::FEATURE_TO_GENE)
    }else{
      feature_to_gene = data.table::data.table(MotrpacRatTraining6moData::FEATURE_TO_GENE_FILT)
    }
    gene_symbol = feature_to_gene[feature_ID == FEATURE_ID, gene_symbol][1]
  }
  
  if(is.null(title)){
    if(add_gene_symbol){
      title = sprintf("%s (%s)", redundant_feature, gene_symbol)
    }else{
      title = redundant_feature
    }
  }else{
    if(add_gene_symbol){
      title = sprintf("%s (%s)", title, gene_symbol)
    }
  }
  
  # melt
  value_cols = colnames(sample_level_data)[grepl("^[0-9]", colnames(sample_level_data))]
  melted_subset = data.table::melt(sample_level_data, 
                                   id.vars=c('feature'), 
                                   measure.vars=value_cols,
                                   variable.name='sample')
  melted_subset = melted_subset[!is.na(value)]
  melted_subset[,sample := as.character(sample)]
  
  meta = unique(data.table::as.data.table(MotrpacRatTraining6moData::PHENO[,c('group','sex','pid','viallabel')]))
  meta[,pid := as.character(pid)]
  
  # use pid or viallabel to merge?
  if(all(melted_subset[,sample] %in% meta[,viallabel])){
    col = 'viallabel'
  }else if(all(melted_subset[,sample] %in% meta[,pid])){
    col = 'pid'
    meta[,viallabel := NULL]
    meta = unique(meta)
  }else{
    stop(sprintf("Sample names in sample-level data do not correspond to vial labels or PIDs: %s...",
         paste(utils::head(melted_subset[,sample]), collapse=", ")))
  }
  
  # merge
  subset_meta = merge(melted_subset, meta, by.x='sample', by.y=col)
  
  bygroup = subset_meta[,list(expr = mean(value, na.rm=T),
                              sd = sd(value, na.rm = T)),
                        by=.(sex, group, feature)]
  
  if(return_data){
    return(as.data.frame(bygroup))
  }
  
  bygroup[,plot_group := sprintf("%s_%s", feature, sex)]
  if(multiple_measurements){
    if(!facet_by_sex){
      g = ggplot2::ggplot(bygroup, ggplot2::aes(x=group, y=expr, group=plot_group, colour=sex)) +
        ggplot2::geom_line(position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::geom_point(aes(shape=feature), size=3, position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = expr-sd, ymax = expr+sd), width=0.2, position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::theme_classic() +
        ggplot2::scale_colour_manual(values=c(female=MotrpacRatTraining6moData::SEX_COLORS[['F']],
                                              male=MotrpacRatTraining6moData::SEX_COLORS[['M']])) +
        ggplot2::labs(x='Time trained (weeks)', y='Normalized value', title=title) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, size=11),
                       legend.title = ggplot2::element_blank(),
                       plot.subtitle = ggplot2::element_text(hjust=0.5),
                       legend.position = "top",
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       legend.margin = ggplot2::margin(t=-5, b=-5, unit="pt"),
                       legend.spacing.y = ggplot2::unit(0, "pt"))  
    }else{
      g = ggplot2::ggplot(bygroup, ggplot2::aes(x=group, y=expr, group=plot_group)) +
        ggplot2::geom_line(colour=MotrpacRatTraining6moData::TISSUE_COLORS[[TISSUE]],
                           position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::geom_point(colour=MotrpacRatTraining6moData::TISSUE_COLORS[[TISSUE]], 
                            aes(shape=feature), 
                            size=3, 
                            position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::geom_errorbar(aes(ymin = expr-sd, ymax = expr+sd), 
                               width=0.2, 
                               colour=MotrpacRatTraining6moData::TISSUE_COLORS[[TISSUE]], 
                               position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::theme_classic() +
        ggplot2::labs(x='Time trained (weeks)', y='Normalized value', title=title) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, size=11),
                       legend.title = ggplot2::element_blank(),
                       plot.subtitle = ggplot2::element_text(hjust=0.5), 
                       legend.position = "top",
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       legend.margin = ggplot2::margin(t=-5, b=-5, unit="pt")) +
        ggplot2::facet_wrap(~sex)
    }
  }else{
    if(!facet_by_sex){
      g = ggplot2::ggplot(bygroup, ggplot2::aes(x=group, y=expr, group=plot_group, colour=sex)) +
        ggplot2::geom_line(position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::geom_point(position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = expr-sd, ymax = expr+sd), width=0.2, position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::theme_classic() +
        ggplot2::scale_colour_manual(values=c(female=MotrpacRatTraining6moData::SEX_COLORS[['F']],
                                              male=MotrpacRatTraining6moData::SEX_COLORS[['M']])) +
        ggplot2::labs(x='Time trained (weeks)', y='Normalized value', title=title) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, size=11),
                       legend.title = ggplot2::element_blank(),
                       plot.subtitle = ggplot2::element_text(hjust=0.5), 
                       legend.position = "top",
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       legend.margin = ggplot2::margin(t=-5, b=-5, unit="pt"),
                       legend.spacing.y = ggplot2::unit(0, "pt"))  
    }else{
      g = ggplot2::ggplot(bygroup, ggplot2::aes(x=group, y=expr, group=plot_group)) +
        ggplot2::geom_line(colour=MotrpacRatTraining6moData::TISSUE_COLORS[[TISSUE]]) +
        ggplot2::geom_point(colour=MotrpacRatTraining6moData::TISSUE_COLORS[[TISSUE]]) +
        ggplot2::geom_errorbar(aes(ymin = expr-sd, ymax = expr+sd), width=0.2, colour=MotrpacRatTraining6moData::TISSUE_COLORS[[TISSUE]]) +
        ggplot2::theme_classic() +
        ggplot2::labs(x='Time trained (weeks)', y='Normalized value', title=title) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, size=11),
                       legend.title = ggplot2::element_blank(),
                       plot.subtitle = ggplot2::element_text(hjust=0.5), 
                       legend.position = "none",
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::facet_wrap(~sex)
    }
  }
    
  if(scale_x_by_time){
    g = g +
      scale_x_discrete(limits=c('control','1w','2w','fill','4w',rep('fill',3), '8w'),
                       labels=c('0','1','2','4','8'),
                       breaks=c('control','1w','2w','4w','8w'))
  }else{
    g = g +
      scale_x_discrete(limits=c('control','1w','2w','4w','8w'),
                       labels=c('0','1','2','4','8'),
                       breaks=c('control','1w','2w','4w','8w'))
  }
  
  if(add_adj_p){
    message("Adding differential analysis p-value...")
    da = plot_feature_logfc(assay = ASSAY,
                            tissue = TISSUE,
                            feature_ID = feature_ID,
                            add_adj_p = TRUE, 
                            return_data = TRUE)
    if(!is.null(da)){
      adj_p_value = min(unique(da$selection_fdr), na.rm=TRUE)
      subtitle = sprintf("adj. p-value: %s", round(adj_p_value, 3))
      g = g + labs(subtitle = subtitle)
    }
  }
  
  return(g)
}

#' Plot differential analysis results for a feature
#' 
#' Plot timewise differential analysis results for a single feature. 
#' Points are log fold-changes, and error bars indicate standard errors. 
#'
#' @param assay NULL or `r assay()`
#' @param tissue NULL or `r tissue()`
#' @param feature_ID NULL or `r feature_ID()`
#' @param feature NULL or `r feature()`. If NULL, \code{assay}, \code{tissue}, and 
#'   \code{feature_ID} must all be specified. 
#' @param title character, plot title. By default, the plot ID is \code{feature}. 
#'   If \code{add_gene_symbol = TRUE}, the gene symbol is also added to the plot title.
#' @param add_gene_symbol bool, whether to add corresponding gene symbol to 
#'   plot title. Default: FALSE 
#' @param facet_by_sex bool, whether to facet the plot by sex. If \code{TRUE},
#'   lines are colored by tissue. If \code{FALSE}, lines are colored by sex. Default: FALSE
#' @param scale_x_by_time bool, whether to scale the x-axis by time. If \code{FALSE},
#'   space the time points (0w, 1w, 2w, 4w, 8w) evenly. Default: TRUE
#' @param add_adj_p bool, whether to include the training adjusted p-value (AKA selection FDR)
#'   in the plot subtitle. Default: TRUE
#' @param metareg bool, whether to use the meta-regression version of differential
#'   analysis results for metabolomics data. If \code{FALSE}, use the redundant,
#'   non-meta-analyzed results. Default: TRUE 
#' @param return_data bool, whether to return data instead of plot. Default: FALSE
#' @param ... additional arguments passed to [get_file_from_url()]
#' 
#' @export
#' 
#' @return a [ggplot2::ggplot()] object or a data frame if \code{return_data = TRUE}
#'    or NULL if the data cannot be found
#' 
#' @examples
#' # 3 ways of plotting the same data are shown in each example below
#' 
#' # Plot a differential feature 
#' plot_feature_logfc(feature = "ACETYL;HEART;NP_001003673.1_K477k",
#'                    add_gene_symbol = TRUE)
#' plot_feature_logfc(assay = "ACETYL",
#'                    tissue = "HEART",
#'                    feature_ID = "NP_001003673.1_K477k",
#'                    add_gene_symbol = TRUE,
#'                    scale_x_by_time = FALSE)
#' plot_feature_logfc(assay = "ACETYL",
#'                    tissue = "HEART",
#'                    feature_ID = "NP_001003673.1_K477k",
#'                    add_gene_symbol = TRUE,
#'                    facet_by_sex = TRUE)
#' 
#' # Plot a redundant differential feature
#' plot_feature_logfc(feature = "IMMUNO;PLASMA;BDNF",
#'                    add_gene_symbol = FALSE)
#' plot_feature_logfc(assay = "IMMUNO",
#'                    tissue = "PLASMA",
#'                    feature_ID = "BDNF",
#'                    add_gene_symbol = FALSE,
#'                    scale_x_by_time = FALSE)
#' plot_feature_logfc(assay = "IMMUNO",
#'                    tissue = "PLASMA",
#'                    feature_ID = "BDNF",
#'                    add_gene_symbol = FALSE,
#'                    facet_by_sex = TRUE)
#'                              
#' # Plot one measurement of a redundant feature
#' plot_feature_logfc(feature = "IMMUNO;PLASMA;rat-myokine:BDNF",
#'                    add_gene_symbol = FALSE)
#' plot_feature_logfc(assay = "IMMUNO",
#'                    tissue = "PLASMA",
#'                    feature_ID = "rat-myokine:BDNF",
#'                    add_gene_symbol = FALSE,
#'                    scale_x_by_time = FALSE)
#' plot_feature_logfc(assay = "IMMUNO",
#'                    tissue = "PLASMA",
#'                    feature_ID = "rat-myokine:BDNF",
#'                    add_gene_symbol = FALSE,
#'                    facet_by_sex = TRUE)
#'                              
#' # Plot a non-differential feature
#' plot_feature_logfc(feature = "PROT;SKM-GN;YP_665629.1",
#'                    add_gene_symbol = TRUE)
#' plot_feature_logfc(assay = "PROT",
#'                    tissue = "SKM-GN",
#'                    feature_ID = "YP_665629.1",
#'                    add_gene_symbol = TRUE,
#'                    scale_x_by_time = FALSE)
#' plot_feature_logfc(assay = "PROT",
#'                    tissue = "SKM-GN",
#'                    feature_ID = "YP_665629.1",
#'                    add_gene_symbol = TRUE,
#'                    facet_by_sex = TRUE)
#'                    
#' # Plot a merged feature from meta-regression
#' plot_feature_logfc(assay = "METAB",
#'                    tissue = "PLASMA",
#'                    feature_ID = "Glucose",
#'                    facet_by_sex = TRUE)
#' plot_feature_logfc(assay = "METAB",
#'                    tissue = "PLASMA",
#'                    feature_ID = "Glucose",
#'                    scale_x_by_time = FALSE,
#'                    metareg = FALSE)
#' plot_feature_logfc(assay = "METAB",
#'                    tissue = "PLASMA",
#'                    feature_ID = "glucose",
#'                    scale_x_by_time = FALSE,
#'                    metareg = FALSE)
#' plot_feature_logfc(assay = "METAB",
#'                    tissue = "PLASMA",
#'                    feature_ID = "glucose",
#'                    facet_by_sex = TRUE,
#'                    metareg = TRUE)
#' plot_feature_logfc(assay = "METAB",
#'                    tissue = "PLASMA",
#'                    feature_ID = "metab-u-ionpneg:glucose",
#'                    scale_x_by_time = FALSE,
#'                    metareg = FALSE)                  
#'    
plot_feature_logfc = function(assay = NULL,
                              tissue = NULL, 
                              feature_ID = NULL,
                              feature = NULL, 
                              title = NULL, 
                              add_gene_symbol = FALSE,
                              facet_by_sex = FALSE, 
                              scale_x_by_time = TRUE, 
                              add_adj_p = TRUE,
                              metareg = TRUE,
                              return_data = FALSE,
                              ...){
  
  curr_feature = feature 
  if(is.null(curr_feature)){
    if(any(is.null(c(assay, tissue, feature_ID)))){
      stop("If 'feature' is not specified, 'assay', 'tissue', and 'feature_ID' must all be specified.")
    }
    curr_feature = sprintf("%s;%s;%s", assay, tissue, feature_ID)
  }
  if(is.null(tissue)){
    splits = unname(unlist(strsplit(curr_feature, ";")))
    assay = splits[1]
    tissue = splits[2]
    feature_ID = splits[3]
  }
  # rename to avoid conflict with data.table columns
  FEATURE_ID = feature_ID
  ASSAY = assay 
  TISSUE = tissue 
  FEATURE = curr_feature 
  redundant_feature = FEATURE
  
  if(ASSAY %in% c("ATAC","METHYL")){
    include_epigen = TRUE
  }else{
    include_epigen = FALSE
  }
  
  # get timewise differential analysis results
  timewise = data.table::data.table(combine_da_results(tissues = TISSUE, assays = ASSAY, metareg = metareg))
  # make feature column
  timewise[is.na(feature), feature := sprintf("%s;%s;%s", assay, tissue, feature_ID)]
  if(!FEATURE %in% timewise[,feature]){
    # check redundant features - timewise results use non-redundant version
    if(FEATURE %in% MotrpacRatTraining6moData::REPEATED_FEATURES$feature){
      FEATURE = MotrpacRatTraining6moData::REPEATED_FEATURES$new_feature[MotrpacRatTraining6moData::REPEATED_FEATURES$feature == FEATURE]
      # now handle FEATURE_ID. get redundant one
      FEATURE_ID = MotrpacRatTraining6moData::REPEATED_FEATURES$feature_ID[MotrpacRatTraining6moData::REPEATED_FEATURES$new_feature == FEATURE[1]]
    # check if this is a metabolite name, in which case there may be multiple results under a refmet ID from meta-reg results 
    }else if(ASSAY == "METAB"){
      if(FEATURE_ID %in% timewise[,feature_ID]){
        FEATURE = unique(timewise[feature_ID == FEATURE_ID & tissue == TISSUE, feature])
      }else if(grepl(":", FEATURE_ID)){
        # or maybe they gave a feature ID and site? e.g. metab-u-ionpneg:glucose
        ds = gsub(":.*","",FEATURE_ID)
        new_feature_id = gsub(".*:","",FEATURE_ID)
        if(nrow(timewise[feature_ID == new_feature_id & dataset == ds & tissue == TISSUE]) > 0){
          timewise = timewise[feature_ID == new_feature_id & dataset == ds & tissue == TISSUE]
          timewise[,feature := sprintf("%s;%s;%s", ASSAY, TISSUE, FEATURE_ID)]
          FEATURE = unique(timewise[,feature])
        }
      }else if(metareg){
        # they could have supplied a metabolite ID that is only available in results as refmet
        if(FEATURE_ID %in% MotrpacRatTraining6moData::METAB_FEATURE_ID_MAP$metabolite_name){
          new_feature_id = unique(MotrpacRatTraining6moData::METAB_FEATURE_ID_MAP$feature_ID_metareg[MotrpacRatTraining6moData::METAB_FEATURE_ID_MAP$metabolite_name == FEATURE_ID])
          new_feature_id = unique(new_feature_id[new_feature_id %in% timewise[,feature_ID]])
          FEATURE = unique(timewise[feature_ID == new_feature_id & tissue == TISSUE, feature])
        }
      }
    }
  }
  curr_timewise_dea = timewise[feature %in% FEATURE]
  
  if(nrow(curr_timewise_dea) == 0){
    warning(sprintf("No DEA results for '%s'.",curr_feature))
    return()
  }
  
  if(add_gene_symbol){
    if(ASSAY %in% c("METHYL","ATAC")){
      feature_to_gene = data.table::data.table(MotrpacRatTraining6moData::FEATURE_TO_GENE)
    }else{
      feature_to_gene = data.table::data.table(MotrpacRatTraining6moData::FEATURE_TO_GENE_FILT)
    }
    gene_symbol = feature_to_gene[feature_ID == FEATURE_ID, gene_symbol][1]
  }
  
  multiple_measurements = FALSE
  ADJ_P = unique(curr_timewise_dea[,selection_fdr])
  if(length(ADJ_P)>1){
    warning(sprintf("Multiple measurements for feature '%s'. Taking the smallest training-dea FDR for the plot label.",curr_feature))
    ADJ_P = min(ADJ_P, na.rm=TRUE)
    curr_timewise_dea[,feature := dataset]
    multiple_measurements = TRUE
  }

  # make title
  if(is.null(title)){
    if(add_gene_symbol){
      title = sprintf("%s (%s)", redundant_feature, gene_symbol)
    }else{
      title = redundant_feature
    }
  }else{
    if(add_gene_symbol){
      title = sprintf("%s (%s)", title, gene_symbol)
    }
  }

  # add 0
  dummy = unique(curr_timewise_dea[,.(tissue, assay, feature)])
  dummy[,logFC := 0]
  dummy[,logFC_se := 0]
  dummy[,comparison_group := 'control']
  dlist = list()
  i = 1
  for(SEX in na.omit(unique(curr_timewise_dea[,sex]))){
    for(f in unique(curr_timewise_dea[,feature])){
      d = copy(dummy)
      d[,sex := SEX]
      d[,feature := f]
      d[,dataset := f]
      dlist[[i]] = d
      i = i+1
    }
  }
  res = rbindlist(c(list(curr_timewise_dea), dlist),fill=TRUE)
  
  if(return_data){
    res[,selection_fdr := ADJ_P]
    return(as.data.frame(res))
  }

  if(multiple_measurements){
    if(facet_by_sex){
      g = ggplot2::ggplot(res, ggplot2::aes(y=logFC,x=comparison_group,group=paste0(tissue,feature),color=tissue))+
        ggplot2::geom_point(size=3, aes(shape=dataset), position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::geom_line(position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin=logFC-logFC_se, ymax=logFC+logFC_se),
                               width=0.2, 
                               position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::theme_classic() +
        ggplot2::geom_hline(yintercept = 0,linetype="dotted") +
        ggplot2::facet_wrap(~sex) +
        ggplot2::labs(title=title, x="Time trained (weeks)", y="Log fold-change") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5),
                       plot.subtitle = ggplot2::element_text(hjust=0.5),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       legend.position = "top",
                       legend.margin = ggplot2::margin(t=-5, b=-5, unit="pt"),
                       legend.spacing.y = ggplot2::unit(0, "pt"),
                       legend.title = ggplot2::element_blank()) +
        ggplot2::scale_colour_manual(values=MotrpacRatTraining6moData::TISSUE_COLORS[names(MotrpacRatTraining6moData::TISSUE_COLORS) %in% res[,tissue]], name="Tissue")
    }else{
      g = ggplot2::ggplot(res, ggplot2::aes(y=logFC,x=comparison_group,group=paste0(tissue,sex,feature),color=sex))+
        ggplot2::geom_point(size=3, aes(shape=dataset), position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::geom_line(position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::geom_errorbar(aes(ymin=logFC-logFC_se, ymax=logFC+logFC_se), width=0.2, position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::theme_classic() +
        ggplot2::geom_hline(yintercept = 0,linetype="dotted") +
        ggplot2::labs(title=title, x="Time trained (weeks)", y="Log fold-change") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5),
              plot.subtitle = ggplot2::element_text(hjust=0.5),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              legend.position = "top",
              legend.margin = ggplot2::margin(t=-5, b=-5, unit="pt"),
              legend.spacing.y = ggplot2::unit(0, "pt"),
              legend.title = ggplot2::element_blank()) +
        ggplot2::scale_colour_manual(values=MotrpacRatTraining6moData::SEX_COLORS[names(MotrpacRatTraining6moData::SEX_COLORS) %in% res[,sex]], name="Sex")
    }
  }else{
    if(facet_by_sex){
      g = ggplot2::ggplot(res, ggplot2::aes(y=logFC,x=comparison_group,group=paste0(tissue,feature),color=tissue))+
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::geom_errorbar(ggplot2::aes(ymin=logFC-logFC_se, ymax=logFC+logFC_se),width=0.2) +
        ggplot2::theme_classic() +
        ggplot2::geom_hline(yintercept = 0,linetype="dotted") +
        ggplot2::facet_wrap(~sex) +
        ggplot2::labs(title=title, x="Time trained (weeks)", y="Log fold-change") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5),
                       plot.subtitle = ggplot2::element_text(hjust=0.5),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       legend.position = "none") +
        ggplot2::scale_colour_manual(values=MotrpacRatTraining6moData::TISSUE_COLORS[names(MotrpacRatTraining6moData::TISSUE_COLORS) %in% res[,tissue]], name="Tissue")
    }else{
      g = ggplot2::ggplot(res, ggplot2::aes(y=logFC,x=comparison_group,group=paste0(tissue,sex,feature),color=sex))+
        ggplot2::geom_point(position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::geom_line(position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::geom_errorbar(aes(ymin=logFC-logFC_se, ymax=logFC+logFC_se),width=0.2,position=ggplot2::position_dodge(width=0.3)) +
        ggplot2::theme_classic() +
        ggplot2::geom_hline(yintercept = 0,linetype="dotted") +
        ggplot2::labs(title=title, x="Time trained (weeks)", y="Log fold-change") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5),
                       plot.subtitle = ggplot2::element_text(hjust=0.5),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       legend.position = "top",
                       legend.margin = ggplot2::margin(t=-5, b=-5, unit="pt"),
                       legend.spacing.y = ggplot2::unit(0, "pt"),
                       legend.title = ggplot2::element_blank()) +
        ggplot2::scale_colour_manual(values=MotrpacRatTraining6moData::SEX_COLORS[names(MotrpacRatTraining6moData::SEX_COLORS) %in% res[,sex]])
    }
  }

  if(add_adj_p){
    subtitle = sprintf("adj. p-value: %s", round(ADJ_P, 3))
    g = g + labs(subtitle = subtitle)
  }

  if(scale_x_by_time){
    g = g +
      scale_x_discrete(limits=c('control','1w','2w','fill','4w',rep('fill',3), '8w'),
                       labels=c('0','1','2','4','8'),
                       breaks=c('control','1w','2w','4w','8w'))
  }else{
    g = g +
      scale_x_discrete(limits=c('control','1w','2w','4w','8w'),
                       labels=c('0','1','2','4','8'),
                       breaks=c('control','1w','2w','4w','8w'))
  }

  return(g)
}


#' Plot feature trajectories 
#' 
#' Plot group means of a set of features from normalized sample-level data.
#' 
#' @param features character vector of features to plot in the format 
#'   '[MotrpacRatTraining6moData::ASSAY_ABBREV];[MotrpacRatTraining6moData::TISSUE_ABBREV];feature_ID'
#' @param exclude_outliers bool, whether to remove sample outliers specified by [MotrpacRatTraining6moData::OUTLIERS].
#'   \code{TRUE} by default. 
#' @param training_regulated_only bool, whether all input features are training-regulated at 5% FDR.
#'   \code{FALSE} by default. If \code{TRUE}, data is loaded from 
#'   [MotrpacRatTraining6moData::TRAINING_REGULATED_NORM_DATA] instead of with [load_sample_data()],
#'   which dramatically improves performance. 
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
#' @details Note that while features in the format 
#'   '[MotrpacRatTraining6moData::ASSAY_ABBREV];[MotrpacRatTraining6moData::TISSUE_ABBREV];feature_ID'
#'   are only given for training-regulated features in the data objects provided by 
#'   \code{MotrpacRatTraining6moData}, one could specify non-training-regulated features by 
#'   concatenating [MotrpacRatTraining6moData::ASSAY_ABBREV], [MotrpacRatTraining6moData::TISSUE_ABBREV],
#'   and \code{feature_ID} for features of interest (semi-colon-separated). 
#' 
#' @examples 
#' # Pick largest cluster in gastrocnemius 
#' clust = extract_tissue_sets("SKM-GN", k=1, add_week8=FALSE)
#' # Extract features 
#' names(clust)
#' features = clust[["1w_F1_M1->2w_F1_M1->4w_F1_M1->8w_F1_M1"]]
#' plot_feature_trajectories(features)
#' 
#' # Since we're only considering training-regulated features in this example,
#' # set training_regulated_only to TRUE to make it slightly faster 
#' plot_feature_trajectories(features, training_regulated_only=TRUE)
#' 
#' # Plot a mix of training-regulated and non-training-regulated features
#' # Note this takes longer because the original datasets are downloaded 
#' features = c(features, "TRNSCRPT;SKM-GN;ENSRNOG00000000008")
#' plot_feature_trajectories(features)
plot_feature_trajectories = function(features, 
                                     training_regulated_only=FALSE, 
                                     exclude_outliers=TRUE, 
                                     center=TRUE, 
                                     scale=TRUE, 
                                     title=NULL, 
                                     return_data=FALSE, 
                                     scratchdir="."){
  
  if(!requireNamespace("viridis", quietly = TRUE) | !requireNamespace("viridisLite", quietly = TRUE)){
    stop(
      "Packages 'viridis' and 'viridisLite' must be installed to run 'plot_feature_trajectories()'.",
      call. = FALSE
    )
  }
  
  features = unique(features)
  
  # this is a default argument, but it's likely that all features are training-regulated 
  # if this is the case, set training_regulated_only to TRUE to make things faster
  if(!training_regulated_only){
    training_reg = unique(MotrpacRatTraining6moData::TRAINING_REGULATED_FEATURES$feature)
    if(length(features) <= length(training_reg)){
      if(all(features %in% training_reg)){
        training_regulated_only = TRUE
      }
    }
  }
  
  if(training_regulated_only){
    if(exclude_outliers){
      data = data.table::as.data.table(fetch_object("TRAINING_REGULATED_NORM_DATA_NO_OUTLIERS"))
    }else{
      data = data.table::as.data.table(fetch_object("TRAINING_REGULATED_NORM_DATA"))
    }
    data = fix_cols(data)
    # select features 
    data = data[feature %in% features]
    if(nrow(data)==0){
      stop("No features in the input match training-regulated features in TRAINING_REGULATED_NORM_DATA.")
    }
    sample_level_data = data
  }else{
    message("Identifying data sets...")
    feature_dt = data.table::as.data.table(check_cluster_res_format(data.frame(feature=features, cluster="dummy")))
    datasets = unique(feature_dt[,.(ome, tissue)])
    omes = unique(datasets[,ome])
    
    # choose which epigenetic dataset to use 
    epigen_training_regulated_only = FALSE
    if("ATAC" %in% omes | "METHYL" %in% omes){
      epigen_input = feature_dt[ome %in% c("ATAC","METHYL"), feature]
      epigen_training_reg = MotrpacRatTraining6moData::TRAINING_REGULATED_FEATURES$feature[MotrpacRatTraining6moData::TRAINING_REGULATED_FEATURES$assay %in% c("METHYL", "ATAC")]
      if(all(epigen_input %in% epigen_training_reg)){
        epigen_training_regulated_only = TRUE
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
                                training_regulated_only=epigen_training_regulated_only, 
                                exclude_outliers=exclude_outliers, 
                                scratchdir=scratchdir,
                                warnings=TRUE)
      }else{
        data = load_sample_data(.tissue, 
                                .ome, 
                                normalized=TRUE, 
                                training_regulated_only=training_regulated_only, 
                                exclude_outliers=exclude_outliers, 
                                scratchdir=scratchdir,
                                warnings=TRUE)
      }
      if(is.null(data)) next
      curr_feat = feature_dt[ome==.ome & tissue==.tissue, feature]
      data = fix_cols(data)
      # add new feature column 
      data[,new_feature := sprintf("%s;%s;%s", assay, tissue, feature_ID)]
      # select features 
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
  }
  
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

fix_cols = function(data){
  data = as.data.frame(data, check.names=FALSE)
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
  data = data.table::as.data.table(data)
  return(data)
}
