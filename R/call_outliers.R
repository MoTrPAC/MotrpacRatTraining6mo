#' Plot 2D scatter plot of principal components 
#' 
#' Plot pairs of principal components and highlight indicated outliers.
#' Each point is a sample. 
#' 
#' @param pcA character, PC for x-axis, e.g. "PC1"; also a column name in \code{pcax}
#' @param pcB character, PC for y-axis, e.g. "PC2"; also a column name in \code{pcax}
#' @param pcax e.g., \code{prcomp(data)$x}
#' @param outliers vector of viallabels corresponding to PC outliers 
#' @param pca result returned by \code{prcomp()}
#' @param title optional character, plot title 
#' 
#' @return ggplot object
#' 
#' @seealso [call_pca_outliers()]
#' 
#' @export
#'
plot_pcs = function(pcA, pcB, pcax, outliers, pca, title=NULL){
  
  if(length(outliers) > 0 & !requireNamespace("ggrepel", quietly = TRUE)){
    stop(
      "Package 'ggrepel' must be installed to run 'plot_pcs()' with outliers.",
      call. = FALSE
    )
  }
  
  pcax = as.data.frame(pcax)
  pcax$viallabel = rownames(pcax)
  
  # get VE 
  explained_var = summary(pca)[["importance"]][2,1:ncol(summary(pca)[["importance"]])]
  xlab = sprintf("%s (%s%%)", pcA, round(explained_var[[pcA]]*100, 1))
  ylab = sprintf("%s (%s%%)", pcB, round(explained_var[[pcB]]*100, 1))
  
  g = ggplot2::ggplot(pcax, aes(x=get(pcA), y=get(pcB))) +
    ggplot2::geom_point() +
    ggplot2::theme_classic() +
    ggplot2::labs(x=xlab, y=ylab, title=title)
  
  if(length(outliers) > 0){
    g = g + ggplot2::geom_point(data=pcax[rownames(pcax) %in% outliers,], size=2, colour='red') +
      ggrepel::geom_text_repel(data=pcax[rownames(pcax) %in% outliers,], ggplot2::aes(label=viallabel)) +
      ggplot2::labs(title=title)
  }
  
  return(g)
}


#' Call PCA outliers
#' 
#' Identify samples that fall outside of the specified range for principal
#' components that explain some minimum variance. 
#' 
#' @param norm feature by sample data frame of normalized data 
#' @param min_pc_ve numeric, minimum percent variance explained by a PC to check it for outliers 
#' @param scale bool, whether to scale input data before PCA. \code{TRUE} by default. 
#' @param plot bool, whether to print PC plots before and after removing outliers. \code{TRUE} by default. 
#' @param verbose bool, whether to print descriptive strings. \code{TRUE} by default.
#' @param iqr_coef numeric, flag PC outliers if they are outside of IQR * \code{iqr_coef}
#' @param M integer, select M most variable features
#' @param title character, substring to include in PC plot titles 
#' 
#' @return named list of four items: 
#' \describe{
#'   \item{\code{pca_outliers}}{character vector of viallabels identified as outliers}
#'   \item{\code{prcomp_obj}}{result returned by \code{prcomp()} from PCA of normalized data without outliers removed}
#'   \item{\code{num_pcs}}{integer, number of PCs checked for outliers}
#'   \item{\code{pc_outliers_report}}{matrix of results with one row per outlier}
#' }
#'
#' @examples
#' bat_rna_data = transcript_prep_data("BAT", covariates = NULL, outliers = NULL)
#' bat_rna_outliers = call_pca_outliers(bat_rna_data$norm_data, 
#'                                      min_pc_ve=0.05, 
#'                                      iqr_coef=5, 
#'                                      M=1000, 
#'                                      title="Brown Adipose")
#'
#' @export
#' @importFrom grDevices boxplot.stats
#'
call_pca_outliers = function(norm, min_pc_ve, scale=TRUE, plot=TRUE, verbose=TRUE, iqr_coef=3, M=Inf, title=NULL){
  
  # check input 
  norm = as.data.frame(norm)
  stopifnot(is.numeric(norm[,1]))
  
  # keep M features with highest CVs
  if(M < nrow(norm)){
    cv = apply(norm, 1, function(x) (stats::sd(x)/mean(x))*100)
    cv = cv[order(cv, decreasing=TRUE)]
    norm = norm[names(cv)[1:M],]
  }
  
  # remove features with 0 variance 
  tnorm = as.data.frame(t(norm))
  novar = names(which(apply(tnorm, 2, stats::var)==0))
  tnorm[,novar] = NULL
  pca = stats::prcomp(x = tnorm,center=TRUE,scale.=scale)
  cum_var = summary(pca)[["importance"]][3,1:ncol(summary(pca)[["importance"]])]
  # save for later
  pca_before = pca
  
  # save as many PCs that explain at least X variance
  explain_vars = summary(pca)[["importance"]][2,1:ncol(summary(pca)[["importance"]])]
  if(verbose){
    print(summary(pca)[["importance"]][2,1:ncol(summary(pca)[["importance"]])][1:10])
  }
  explain_var = explain_vars[explain_vars > min_pc_ve]
  num_pcs = length(explain_var)
  if(num_pcs == 0){
    num_pcs = 1
    explain_var[explain_vars[1]]
  }
  
  if(verbose){
    message(sprintf("The first %s PCs were selected to identify outliers.", num_pcs))
  }
  pcax = pca$x[,1:max(2, num_pcs)]
  
  # ID outliers 
  # Univariate: use IQRs
  pca_outliers_report = c()
  pca_outliers = c()
  for(j in 1:num_pcs){
    outlier_values = grDevices::boxplot.stats(pcax[,j],coef=iqr_coef)$out # flag samples beyond IQR*iqr_coef
    for(outlier in names(outlier_values)){
      pca_outliers_report = rbind(pca_outliers_report,
                                  c(paste("PC",j,sep=""),outlier,
                                    unname(format(outlier_values[outlier],digits=5)))
      )
      if(!is.element(outlier,names(pca_outliers))){
        pca_outliers[outlier] = outlier_values[outlier]
      }
    }
  }
  if(verbose){
    if(length(pca_outliers) == 0){
      message("No PC outliers.")
    }else{
      message("PC outliers:")
      colnames(pca_outliers_report) = c('PC','viallabel','PC value')
      print(pca_outliers_report)
    }
  }
  
  if(plot){
    # plot before
    if(length(pca_outliers) > 0){
      if(verbose) message("Plotting PCs with outliers flagged...")
    }else{
      if(verbose) message("Plotting PCs. No outliers flagged.")
    }
    for(i in 2:max(2,num_pcs)){
      print(plot_pcs("PC1", paste0("PC", i), pcax, names(pca_outliers), pca, title=sprintf('%s before outlier removal',title)))
    }
    
    # plot after 
    if(length(pca_outliers) > 0){
      if(verbose) message("Plotting PCs with outliers removed...")
      # plot after 
      filt_norm = norm[,!colnames(norm) %in% names(pca_outliers)]
      # remove features with 0 variance 
      tnorm <- as.data.frame(t(filt_norm))
      novar <- names(which(apply(tnorm, 2, stats::var)==0))
      tnorm[,novar] = NULL
      pca = stats::prcomp(x = tnorm,center=TRUE,scale.=TRUE)
      pcax = pca$x[,1:max(2, num_pcs)]
      for(i in 2:max(2, num_pcs)){
        print(plot_pcs("PC1", paste0("PC", i), pcax, c(), pca, title=sprintf('%s after outlier removal',title)))
      }
    }
  }
  
  # return outlier list 
  return(list(pca_outliers = names(pca_outliers),
              prcomp_obj = pca_before,
              num_pcs = num_pcs,
              pc_outliers_report = pca_outliers_report))
}


#' Call RNA-seq sample outliers
#' 
#' Identify samples that are outside of 5 times the interquartile range of principal
#' components that explain at least 5% of variance in each tissue. Use only the 1000 most variable genes. 
#' This specifies RNA-seq outliers excluded from differential analysis by MoTrPAC.  
#' 
#' @param tissues character vector of tissue abbreviations for which to call RNA-seq outliers. 
#'     See [MotrpacRatTraining6moData::TISSUE_ABBREV] for accepted values. 
#' 
#' @return NULL if there are no outliers, or a data frame with three columns and one row per outlier:
#' \describe{
#'   \item{\code{viallabel}}{character, `r viallabel()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{reason}}{character, PC(s) in which the sample was flagged}
#' }
#'
#' @seealso [call_pca_outliers()] for workhorse function and [plot_pcs()] for plotting function
#'
#' @export
#'
#' @examples
#' transcript_call_outliers("SKM-VL")
#' \dontrun{
#' transcript_call_outliers(c("SKM-GN","BLOOD"))
#' }
transcript_call_outliers = function(tissues){
  pca_outliers_list = list()
  for(.tissue in tissues){
    
    data = transcript_prep_data(.tissue, covariates = NULL, outliers = NULL)
    meta = data.table::data.table(data$metadata)
    counts = data$filt_counts
    tmm = data$norm_data
    
    message(sprintf("%s:",.tissue))
    tmm_pca_1k = call_pca_outliers(tmm, 0.05, iqr_coef=5, M=1000, title=.tissue)
    
    if(length(tmm_pca_1k$pca_outliers)>0){
      rep = data.table::as.data.table(tmm_pca_1k$pc_outliers_report)
      rep[,tissue := .tissue]
      pca_outliers_list[[.tissue]] = rep
    }
    
  }
  pca_outliers = rbindlist(pca_outliers_list)
  if(nrow(pca_outliers)==0){
    return(NULL)
  }
  out = pca_outliers[,list(tissue=.tissue,
                           reason=paste0(PC, collapse=',')),
                     by = viallabel]
  out = unique(out)
  return(as.data.frame(out))
}


#' Call ATAC-seq sample outliers
#' 
#' Identify samples that are outside of 3 times the interquartile range of principal
#' components that explain at least 7.5% of variance in each tissue. Use only the 10,000 most variable peaks. 
#' Outlier calling is performed separately in each sex.  
#' This specifies ATAC-seq outliers excluded from differential analysis by MoTrPAC.  
#' 
#' @param tissues character vector of tissue abbreviations for which to call ATAC-seq outliers, 
#'   one of "BAT", "HEART", "HIPPOC", "KIDNEY", "LIVER", "LUNG", "SKM-GN", "WAT-SC" 
#' @param scratchdir character, local directory in which to download data from Google Cloud Storage. 
#'   Current working directory by default.
#' 
#' @return NULL if there are no outliers, or a data frame with three columns and one row per outlier:
#' \describe{
#'   \item{\code{viallabel}}{character, `r viallabel()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{reason}}{character, PC(s) in which the sample was flagged}
#' }
#'
#' @seealso [call_pca_outliers()] for workhorse function and [plot_pcs()] for plotting function
#'
#' @export 
#'
#' @examples
#' \dontrun{
#' atac_call_outliers("LIVER","/tmp")
#' atac_call_outliers(c("LIVER","BAT"),"/tmp")
#' }
atac_call_outliers = function(tissues, scratchdir = "."){
  pca_outliers_list = list()
  for(.tissue in tissues){
    
    data = atac_prep_data(.tissue, covariates = NULL, outliers = NULL, return_normalized_data = TRUE, scratchdir = scratchdir)
    meta = data.table::data.table(data$metadata)
    tmm = data$norm_data
    
    for(s in unique(meta[,sex])){
      dataset = sprintf("%s %s", .tissue, s)
      message(dataset)
      sub_quant = tmm[,meta[sex == s,viallabel]]
      pca_out = call_pca_outliers(sub_quant, 0.075, iqr_coef=3, title=dataset, M=10000)
      if(length(pca_out$pca_outliers)>0){
        rep = data.table::as.data.table(pca_out$pc_outliers_report)
        rep[,tissue := .tissue]
        pca_outliers_list[[dataset]] = rep
      }
    }
  }
  
  pca_outliers = rbindlist(pca_outliers_list)
  if(nrow(pca_outliers)==0){
    return(NULL)
  }
  out = pca_outliers[,list(reason=paste0(PC, collapse=',')),
                     by = .(tissue, viallabel)]
  out = unique(out)
  return(as.data.frame(out))
}
