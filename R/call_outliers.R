#' Plot 2D scatter plot of principal components 
#' 
#' Internal function to plot pairs of principal components and highlight indicated outliers.
#' Each point is a sample. 
#' 
#' @param pcaA string of PC for x-axis, e.g. "PC1"; also a column name in \code{pcax}
#' @param pcaB string of PC for y-axis, e.g. "PC2"; also a column name in \code{pcax}
#' @param pcax e.g., \code{prcomp(data)$x}
#' @param outliers list of viallabels corresponding to PC outliers 
#' @param pca result returned by \code{prcomp()}
#' @param title optional plot title 
#' 
#' @return ggplot object
#' 
#' @seealso [call_pca_outliers()]
#' 
#' @examples
#' TODO
#' 
#' @import ggplot2 
#'
plot_pcs = function(pcA, pcB, pcax, outliers, pca, title=NULL){
  
  pcax = as.data.frame(pcax)
  pcax$viallabel = rownames(pcax)
  
  # get VE 
  explained_var = summary(pca)[["importance"]][2,1:ncol(summary(pca)[["importance"]])]
  xlab = sprintf("%s (%s%%)", pcA, round(explained_var[[pcA]]*100, 1))
  ylab = sprintf("%s (%s%%)", pcB, round(explained_var[[pcB]]*100, 1))
  
  g = ggplot(pcax, aes(x=get(pcA), y=get(pcB))) +
    geom_point() +
    theme_classic() +
    labs(x=xlab, y=ylab, title=title)
  
  if(length(outliers) > 0){
    g = g + geom_point(data=pcax[rownames(pcax) %in% outliers,], size=2, colour='red') +
      geom_text_repel(data=pcax[rownames(pcax) %in% outliers,], aes(label=viallabel)) +
      labs(title=title)
  }
  
  return(g)
}


#' Use PCA method to call outliers
#' 
#' TODO: Description
#' 
#' @param norm filtered, normalized values 
#' @param min_pc_ve minimum percent variance explained by a PC to check it for outliers 
#' @param plot bool, whether or not to print plots
#' @param verbose bool, whether or not to print descriptive strings
#' @param iqr_coef numeric, flag PC outliers if they are outside of IQR * \code{iqr_coef}
#' @param N integer, select N most variable features
#' @param TITLE string, substring to include in PC plot titles 
#' 
#' @return TODO 
#'
#' @examples
#' TODO
#'
#' @export
#' @import data.table
#' @importFrom grDevices boxplot.stats
#'
call_pca_outliers = function(norm, min_pc_ve, plot, verbose, iqr_coef=3, N=Inf, TITLE=NULL){
  
  # keep N features with highest CVs
  if(N < nrow(norm)){
    cv = apply(norm, 1, function(x) (sd(x)/mean(x))*100)
    cv = cv[order(cv, decreasing=T)]
    norm = norm[names(cv)[1:N],]
  }
  
  # remove features with 0 variance 
  tnorm = as.data.frame(t(norm))
  novar = names(which(apply(tnorm, 2, stats::var)==0))
  tnorm[,novar] = NULL
  pca = prcomp(x = tnorm,center=T,scale.=T)
  cum_var = summary(pca)[["importance"]][3,1:ncol(summary(pca)[["importance"]])]
  # save for later
  pca_before = pca
  
  # save as many PCs that explain at least X variance
  explain_vars = summary(pca)[["importance"]][2,1:ncol(summary(pca)[["importance"]])]
  if(verbose){print(summary(pca)[["importance"]][2,1:ncol(summary(pca)[["importance"]])][1:10])}
  explain_var = explain_vars[explain_vars > min_pc_ve]
  num_pcs = length(explain_var)
  if(num_pcs == 0){
    num_pcs = 1
    explain_var[explain_vars[1]]
  }
  
  if(verbose){cat(sprintf("The first %s PCs were selected to identify outliers.\n", num_pcs))}
  pcax = pca$x[,1:max(2, num_pcs)]
  
  # ID outliers 
  # Univariate: use IQRs
  pca_outliers_report = c()
  pca_outliers = c()
  for(j in 1:num_pcs){
    outlier_values = boxplot.stats(pcax[,j],coef=iqr_coef)$out # flag samples beyond IQR*iqr_coef
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
      cat("No PC outliers.\n")
    }else{
      cat("PC outliers:\n")
      colnames(pca_outliers_report) = c('PC','viallabel','PC value')
      print(pca_outliers_report)
    }
  }
  
  if(plot){
    # plot before
    cat("PC plots with any outliers flagged:\n")
    for(i in 2:max(2,num_pcs)){
      print(plot_pcs("PC1", paste0("PC", i), pcax, names(pca_outliers), pca, title=sprintf('%s before outlier removal',TITLE)))
    }
    if(length(pca_outliers) > 0){
      cat("PC plots with outliers removed:\n")
      # plot after 
      filt_norm = norm[,!colnames(norm) %in% names(pca_outliers)]
      # remove features with 0 variance 
      tnorm <- as.data.frame(t(filt_norm))
      novar <- names(which(apply(tnorm, 2, stats::var)==0))
      tnorm[,novar] = NULL
      pca = prcomp(x = tnorm,center=T,scale.=T)
      pcax = pca$x[,1:max(2, num_pcs)]
      for(i in 2:max(2, num_pcs)){
        print(plot_pcs("PC1", paste0("PC", i), pcax, c(), pca, title=sprintf('%s after outlier removal',TITLE)))
      }
    }
  }
  
  # return outlier list 
  return(list(pca_outliers = names(pca_outliers),
              prcomp_obj = pca_before,
              num_pcs = num_pcs,
              pc_outliers_report = pca_outliers_report))
}


#' Call RNA-seq outliers
#' 
#' Description
#' 
#' @param tissues list of tissue abbreviations for which to call RNA-seq outliers. 
#'     See [TISSUE_ABBREV()] for accepted values. 
#' 
#' @return something
#'
#' @seealso [call_pca_outliers()] for workhorse function, [plot_pcs()] for plotting function, 
#'     [TISSUE_ABBREV()] for list of accepted tissue abbrevations 
#'
#' @examples
#' transcript_call_outliers("SKM-GN")
#' transcript_call_outliers(c("SKM-GN","BLOOD"))
#'
#' @export
#' @import MotrpacBicQC
#' @import data.table
#'
transcript_call_outliers = function(tissues){
  pca_outliers_list = list()
  for(tissue in tissues){
    
    data = preprocess_pass1b_rnaseq_gcp(tissue_code, 'all', gsutil_path = gsutil)
    if(is.null(data)){next}
    meta = data$meta
    counts = data$counts
    tmm = data$norm
    meta[,sex_group := paste0(sex, ';', group)]
    
    tmm_pca_1k = call_pca_outliers(tmm, 0.05, plot=T, verbose=T, iqr_coef=5, N=1000, TITLE=tissue_code)
    
    if(length(tmm_pca_1k$pca_outliers)>0){
      rep = as.data.table(tmm_pca_1k$pc_outliers_report)
      rep[,tissue := tissue_code]
      pca_outliers_list[[tissue_code]] = rep
    }
    
  }
  pca_outliers = rbindlist(pca_outliers_list)
  pca_outliers
  out = pca_outliers[,list(dataset=tissue_code,
                           reason=paste0(PC, collapse=',')),
                     by=viallabel]
  out=unique(out)
  return(out)
}
