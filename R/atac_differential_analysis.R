#' ATAC-seq training differential analysis 
#' 
#' Use limma to perform an F-test to test the effect of training
#' across time points. Analysis is performed separately for males and females. 
#'
#' @param tissue `r tissue()`
#' @param covariates character vector of covariates that correspond to column names of [MotrpacRatTraining6moData::ATAC_META].
#'   Defaults to covariates that were used for the manuscript. 
#' @param outliers vector of viallabels to exclude during differential analysis. Defaults
#'   to \code{[MotrpacRatTraining6moData::OUTLIERS]$viallabel}
#' @param scratchdir character, local directory in which to download data from 
#'   Google Cloud Storage. Current working directory by default. 
#' @param rdata_outfile NULL or path in which to save eBayes objects in an RData file 
#' @param overwrite boolean, whether to overwrite the file if \code{rdata_outfile} exists
#' @param verbose boolean, whether to print messages
#' @param n_features integer, number of peaks to perform differential analysis on. Useful for troubleshooting. 
#'
#' @return data frame with one row per feature:
#' \describe{
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{removed_samples_male}}{character, comma-separated list of male outliers (vial labels) removed from differential analysis}
#'   \item{\code{removed_samples_female}}{character, comma-separated list of female outliers (vial labels) removed from differential analysis}
#'   \item{\code{fscore_male}}{double, F statistic for males}
#'   \item{\code{fscore_female}}{double, F statistic for females}
#'   \item{\code{p_value_male}}{double, nominal F-test p-value for males}
#'   \item{\code{p_value_female}}{double, nominal F-test p-value for females}
#'   \item{\code{full_model_male}}{character, full model used in F-test for males}
#'   \item{\code{full_model_female}}{character, full model used in F-test for females}
#'   \item{\code{reduced_model_male}}{character, reduced model used in F-test for males}
#'   \item{\code{reduced_model_female}}{character, reduced model used in F-test for females}
#'   \item{\code{p_value}}{double, combined male and female nominal p-value using the sum of logs}
#' }
#' 
#' @export
#' @importFrom metap sumlog
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom data.table data.table as.data.table rbindlist
#' @import MotrpacRatTraining6moData
#'
#' @examples
#' 
#' \dontrun{
#' # Perform differential analysis for chromatin accessibility peaks in 
#' # brown adipose tissue with default parameters, 
#' # i.e., outliers and covariates used for the manuscript
#' da = atac_training_da("BAT")
#' 
#' # Same as above but save the [limma::eBayes] MArrayLM objects in an RData file 
#' da = atac_training_da("BAT", rdata_outfile = "~/test/BAT_ATAC_training-da.RData", overwrite = TRUE)
#' }
atac_training_da = function(tissue, 
                            covariates = c("Sample_batch", "peak_enrich.frac_reads_in_peaks.macs2.frip"), 
                            outliers = na.omit(MotrpacRatTraining6moData::OUTLIERS$viallabel[MotrpacRatTraining6moData::OUTLIERS$assay == "ATAC"]),
                            scratchdir = ".",
                            rdata_outfile = NULL,
                            overwrite = FALSE,
                            verbose = TRUE,
                            n_features = Inf){
  # data.table workaround
  .tissue = tissue
  
  check_da_args(.tissue, rdata_outfile, overwrite)
  
  if(verbose) message("Loading data...")
  data = atac_prep_data(tissue, 
                        covariates = covariates,
                        filter_counts = TRUE,
                        return_normalized_data = FALSE, 
                        scratchdir = scratchdir, 
                        outliers = outliers)
  
  meta_df = as.data.frame(data$metadata)
  if("sample_batch" %in% tolower(covariates)){
    which = covariates[grepl("sample_batch", covariates, ignore.case=T)]
    meta_df[,which] = factor(meta_df[,which])
  }
  meta_df$group = factor(meta_df$group, levels=c('control','1w','2w','4w','8w')) # IMPORTANT
  filt_counts = data$counts
  filt_counts = filt_counts[1:n_features,]
  
  # split by sex 
  sex_res = list()
  ebayes_list = list()
  for(SEX in c('male','female')){
    
    if(verbose) message(sprintf("Performing F-tests for %s %ss...", tissue, SEX))
    
    # subset counts and meta
    curr_meta = meta_df[meta_df$sex==SEX,]
    curr_counts = filt_counts[,curr_meta$viallabel]
    curr_outliers = filter_outliers(TISSUE=tissue, SEX=SEX, outliers=data$outliers)
    
    full = paste0(c("~ 1", "group", covariates), collapse=" + ")
    reduced = paste0(c("~ 1", covariates), collapse=" + ")
    
    # normalize and get voom weights 
    design = model.matrix(eval(parse(text=full)), data = curr_meta)
    # check if full rank
    if(!is.fullrank(design)){
      if("sample_batch" %in% tolower(covariates)){
        which = covariates[grepl("sample_batch", covariates, ignore.case=T)]
        curr_cov = covariates[!covariates == which]
        full = paste0(c("~ 1", "group", curr_cov), collapse=" + ")
        reduced = paste0(c("~ 1", curr_cov), collapse=" + ")
        design = model.matrix(eval(parse(text=full)), data = curr_meta)
        warning(sprintf("Sample_batch and group or sex are confounded for %s %s. Removing Sample_batch as a covariate.", tissue, SEX))
      }else{
        stop(sprintf("Model matrix with design %s is not full rank.", full))
      }
    }else{
      curr_cov = covariates
    }
    curr_voom = voom(curr_counts, design=design, normalize.method="quantile")
    
    limma_model1 = lmFit(curr_voom, design)
    eb_Ftest = eBayes(limma_model1)
    ebayes_list[[SEX]] = eb_Ftest
    res = topTable(eb_Ftest, n=Inf, coef=colnames(design)[grepl('group',colnames(design))], sort.by = "none")
    dt = data.table(feature_ID=rownames(res),
                    assay="ATAC",
                    assay_code='epigen-atac-seq',
                    tissue=tissue,
                    tissue_code=TISSUE_ABBREV_TO_CODE[[tissue]], 
                    removed_samples=ifelse(length(curr_outliers)>0, paste0(curr_outliers, collapse=','), NA_character_),
                    fscore=res$`F`,
                    p_value = res$P.Value,
                    full_model=gsub(' ','',full),
                    reduced_model=gsub(' ','',reduced))
    sex_res[[SEX]] = dt
  }
  
  # save to file 
  if(!is.null(rdata_outfile)){
    if(overwrite | (!overwrite & !file.exists(rdata_outfile))){
      save(ebayes_list, file=rdata_outfile)
      if(verbose) message(sprintf("'ebayes_list' saved in 'rdata_outfile': %s", rdata_outfile))
    }
  }
  
  # merge 
  merged = data.table(merge(sex_res[['male']], sex_res[['female']], 
                            by=c("feature_ID","assay","assay_code","tissue","tissue_code"), 
                            suffixes=c('_male','_female')))
  
  # get a single meta p-value per feature using the male- and female- specific p-values 
  merged[,p_value := sumlog(c(p_value_male, p_value_female))$p, by=seq_len(nrow(merged))]

  # reorder columns
  merged = merged[,.(
    feature_ID,
    assay,
    assay_code,
    tissue,
    tissue_code,
    removed_samples_male, 
    removed_samples_female,
    fscore_male,
    fscore_female,
    p_value_male,
    p_value_female,
    full_model_male,
    full_model_female,
    reduced_model_male,
    reduced_model_female,
    p_value
  )]
  
  merged = as.data.frame(merged)
  
  if(verbose) message("Done.")
  return(merged)
}


#' ATAC-seq timewise differential analysis 
#' 
#' Use limma to perform pairwise contrasts between each group of trained animals
#' and the sex-matched control group for a single tissue. Analysis is performed separately for males and 
#' females. 
#' 
#' @param tissue `r tissue()`
#' @param covariates character vector of covariates that correspond to column names of [MotrpacRatTraining6moData::ATAC_META].
#'   Defaults to covariates that were used for the manuscript. 
#' @param outliers vector of viallabels to exclude during differential analysis. Defaults
#'   to \code{\link[MotrpacRatTraining6moData]{OUTLIERS}$viallabel}
#' @param scratchdir character, local directory in which to download data from 
#'   Google Cloud Storage. Current working directory by default. 
#' @param rdata_outfile NULL or path in which to save eBayes objects in an RData file 
#' @param overwrite boolean, whether to overwrite the file if \code{rdata_outfile} exists
#' @param verbose boolean, whether to print messages
#' @param n_features integer, number of peaks to perform differential analysis on. Useful for troubleshooting. 
#'
#' @return a data frame with one row per feature:
#' \describe{
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{sex}}{`r sex()`}
#'   \item{\code{comparison_group}}{`r comparison_group()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{covariates}}{character, comma-separated list of adjustment variables}
#'   \item{\code{removed_samples}}{character, comma-separated list of outliers (vial labels) removed from differential analysis}
#'   \item{\code{logFC}}{`r logFC()`}
#'   \item{\code{logFC_se}}{`r logFC_se()`}
#'   \item{\code{tscore}}{double, t statistic}
#'   \item{\code{p_value}}{`r p_value_da()`}
#' }
#' 
#' @export
#' @importFrom limma voom lmFit eBayes topTable is.fullrank makeContrasts contrasts.fit
#' @importFrom data.table data.table as.data.table rbindlist
#' @import MotrpacRatTraining6moData
#'
#' @examples
#' \dontrun{
#' # Perform timewise differential analysis for chromatin accessibility peaks measured in hippocampus with default parameters.
#' hippoc_atac_da = atac_timewise_da("HIPPOC")
#' }
atac_timewise_da = function(tissue, 
                            covariates = c("Sample_batch", "peak_enrich.frac_reads_in_peaks.macs2.frip"), 
                            outliers = na.omit(MotrpacRatTraining6moData::OUTLIERS$viallabel[MotrpacRatTraining6moData::OUTLIERS$assay == "ATAC"]),
                            scratchdir = ".",
                            rdata_outfile = NULL,
                            overwrite = FALSE,
                            verbose = TRUE,
                            n_features = Inf){
  
  check_da_args(.tissue, rdata_outfile, overwrite)
  if(verbose) message("Loading data...")
  data = atac_prep_data(tissue, 
                        covariates = covariates,
                        filter_counts = TRUE,
                        return_normalized_data = FALSE, 
                        scratchdir = scratchdir, 
                        outliers = outliers)
  
  meta_df = as.data.frame(data$metadata)
  if("sample_batch" %in% tolower(covariates)){
    which = covariates[grepl("sample_batch", covariates, ignore.case=T)]
    meta_df[,which] = factor(meta_df[,which])
  }
  filt_counts = data$counts
  filt_counts = filt_counts[1:n_features,]
  
  # split by sex 
  sex_res = list()
  ebayes_list = list()
  for(SEX in c('male','female')){
    
    if(verbose) message(sprintf("Performing differential analysis for %s %ss...", tissue, SEX))
    
    curr_meta = meta_df[meta_df$sex==SEX,]
    curr_counts = filt_counts[,curr_meta$viallabel]
    curr_outliers = filter_outliers(TISSUE=tissue, SEX=SEX, outliers=data$outliers)
    
    full = paste0(c("~ 0", "group", covariates), collapse=" + ")
    reduced = paste0(c("~ 0", covariates), collapse=" + ")
    
    # normalize and get voom weights 
    design = model.matrix(eval(parse(text=full)), data = curr_meta)
    # check if full rank
    if(!is.fullrank(design)){
      if("sample_batch" %in% tolower(covariates)){
        which = covariates[grepl("sample_batch", covariates, ignore.case=T)]
        curr_cov = covariates[!covariates == which]
        full = paste0(c("~ 0", "group", curr_cov), collapse=" + ")
        reduced = paste0(c("~ 0", curr_cov), collapse=" + ")
        design = model.matrix(eval(parse(text=full)), data = curr_meta)
        warning(sprintf("Sample_batch and group or sex are confounded for %s %s", tissue, SEX))
      }else{
        stop(sprintf("Model matrix with design %s is not full rank.", full))
      }
    }else{
      curr_cov = covariates
    }
    curr_voom = voom(curr_counts, design=design, normalize.method="quantile")
    fit = lmFit(curr_voom,design)
    
    cont.matrix=makeContrasts(
      '1W'='group1w - groupcontrol',
      '2W'='group2w - groupcontrol',
      '4W'='group4w - groupcontrol',
      '8W'='group8w - groupcontrol',
      levels=design)
    
    fit2=contrasts.fit(fit,cont.matrix)
    e=eBayes(fit2)
    ebayes_list[[SEX]] = e
    
    da_list = list()
    for(tp in c('1W','2W','4W','8W')){
      res = topTable(e, number=nrow(e), coef=tp, confint=TRUE)
      dt = data.table(
        feature_ID = rownames(res),
        sex = SEX, 
        comparison_group = tolower(tp),
        assay = "ATAC",
        assay_code = "epigen-atac-seq",
        tissue = tissue, 
        tissue_code = TISSUE_ABBREV_TO_CODE[[tissue]],
        covariates=paste0(curr_cov, collapse=','),
        removed_samples = ifelse(length(curr_outliers)>0, paste0(curr_outliers, collapse=','), NA_character_), 
        logFC = res$logFC,
        logFC_se = (res$CI.R - res$CI.L)/3.92,
        tscore = res$t,
        p_value = res$P.Value
      )
      da_list[[tp]] = dt
    }
    sex_res[[SEX]] = rbindlist(da_list)
    
  }
  # save to file 
  if(!is.null(rdata_outfile)){
    if(overwrite | (!overwrite & !file.exists(rdata_outfile))){
      save(ebayes_list, file=rdata_outfile)
      if(verbose) message(sprintf("'ebayes_list' saved in 'rdata_outfile': %s", rdata_outfile))
    }
  }
  
  res = rbindlist(sex_res)
  
  # reorder columns
  res = res[,.(
    feature_ID,
    assay,
    assay_code,
    tissue,
    tissue_code,
    covariates,
    removed_samples,
    logFC,
    logFC_se,
    tscore,
    p_value
  )]
  
  if(verbose) message("Done.")
  return(as.data.frame(res))
}
