#' Wrapper for [DESeq2::DESeq()]
#' 
#' TODO
#' 
#' @param counts data frame of raw filtered RNA-seq counts. Row names are gene IDs and 
#'   column names are sample IDs. Column names must correspond to values in \code{meta$viallabel}
#' @param meta data frame of metadata with columns \code{\param{outcome_of_interest}, \param{covar}, 'viallabel'} at a minimum. 
#'   \param{counts} are subset to \code{meta$viallabel}
#' @param covar adjustment variables to include in the DESeq model
#' @param outcome_of_interest outcome of interest to include in the model; 
#'   levels provided in \param{contrasts}
#' @param constrasts list of vectors, where each vector is in the form 
#'   \code(c(outcome_of_interest, numerator, denominator)), e.g. \code{c('sex_group','female.1w','female.control')}
#' @param shrink bool, whether to apply \code{lfcShrink()}
#' @param verbose bool, whether to print the design string 
#' 
#' @return TODO
#' 
#' @export 
#' @import data.table
#' @import DESeq2
#' @import ashr
#' 
#' @examples 
#' TODO
#' 
run_deseq = function(counts, meta, covar, outcome_of_interest, contrasts, shrink=TRUE, verbose=FALSE){
  
  meta = as.data.table(meta)
  meta[,(outcome_of_interest) := as.factor(get(outcome_of_interest))]
  counts = counts[,as.character(meta[,viallabel])]
  
  # coerce to counts (RSEM does something weird)
  counts = as.data.frame(apply(counts, c(1,2), as.integer)) 
  
  # remove missing values; center and scale covariates
  new = fix_covariates(covar, meta, center_scale = TRUE)
  covar = new$covariates
  meta = data.table(new$meta)
  
  # make contrast
  contrast = paste0('~', paste0(c(outcome_of_interest, covar), collapse=' + '))
  if(verbose) message(contrast)
  
  # run DESeq
  dds = DESeqDataSetFromMatrix(countData = counts,
                               colData = meta,
                               design = eval(parse(text=contrast)))
  dds = DESeq(dds, quiet = T)
  
  # get results for each contrast 
  res_list = list()
  for (c in contrasts){
    if(!shrink){
      res = results(dds, contrast = c)
    }else{
      res = lfcShrink(dds, contrast = c, type = 'ashr', quiet = T, 
                      control=list(numiter.em=1000), optmethod = 'mixSQP') # "apeglm" doesn't work with contrasts in this form 
    }
    res_dt = data.table(gene_id = rownames(counts), 
                        log2FoldChange = res$log2FoldChange,
                        lfcSE = res$lfcSE,
                        pvalue = res$pvalue,
                        numerator = c[2],
                        denominator = c[3],
                        covariates = paste0(covar, collapse=','))
    if(!shrink){
      res_dt[,stat := res$stat]
    }
    res_list[[paste0(c, collapse=' ')]] = res_dt
  }
  all_res = rbindlist(res_list)
  return(list(res=as.data.frame(all_res),
              dds=dds))
}


#' RNA-seq timewise differential analysis 
#' 
#' Use DESeq2 to perform pairwise contrasts between each group of trained animals
#' and the sex-matched control group for a single tissue. Analysis is performed separately for males and 
#' females. 
#'
#' @param tissue character, tissue abbreviation, one of [MotrpacRatTraining6moData::TISSUE_ABBREV]
#' @param covariates character vector of covariates that correspond to column names of [MotrpacRatTraining6moData::TRNSCRPT_META].
#'   Defaults to covariates that were used for the manuscript. 
#' @param outliers vector of viallabels to exclude during differential analysis. Defaults
#'   to \code{MotrpacRatTraining6moData::OUTLIERS$viallabel[MotrpacRatTraining6moData::OUTLIERS$assay == "TRNSCRPT"]}
#' @param add_shrunk_logfc boolean, whether to calculate shrunk log fold-changes in addition to standard log fold-changes
#' @param rdata_outfile NULL or path in which to save DESeq2 objects in an RData file 
#' @param overwrite boolean, whether to overwrite the file if \code{rdata_outfile} exists
#' @param verbose boolean, whether to print the DESeq2 design string
#'
#' @return a data frame with one row per gene per contrast (usually 8 rows per gene):
#' \describe{
#'   \item{\code{feature}}{unique gene identifier in the format \code{[ASSAY_ABBREV];[TISSUE_ABBREV];[feature_ID]}}
#'   \item{\code{feature_ID}}{Ensembl gene ID}
#'   \item{\code{sex}}{one of "male" or "female"}
#'   \item{\code{comparison_group}}{time point compared to the sex-matched 
#'       sedentary control animals, one of "1w", "2w", "4w", "8w"}
#'   \item{\code{assay}}{assay abbreviation, one of [MotrpacRatTraining6moData::ASSAY_ABBREV]}
#'   \item{\code{assay_code}}{MoTrPAC assay or "ome" code. "transcript-rna-seq" for RNA-seq datasets.}
#'   \item{\code{tissue}}{tissue abbreviation, one of [MotrpacRatTraining6moData::TISSUE_ABBREV]}
#'   \item{\code{tissue_code}}{MoTrPAC tissue release code. See [MotrpacBicQC::bic_animal_tissue_code] for details.}
#'   \item{\code{covariates}}{comma-separated list of adjustment variables}
#'   \item{\code{removed_samples}}{comma-separated list of outliers (vial labels) removed from differential analysis}
#'   \item{\code{logFC}}{log fold-change of the training group specified by \code{sex} and \code{comparison_group} (e.g., 1-week females) 
#'     relative to the sex-matched sedentary controls}
#'   \item{\code{logFC_se}}{standard error of \code{logFC}}
#'   \item{\code{shrunk_logFC}}{log fold-change shrunk with \code{type = 'ashr'} and \code{optmethod = 'mixSQP'}, only if \code{add_shrunk_logfc = TRUE}}
#'   \item{\code{shrunk_logFC_se}}{standard error of \code{shrunk_logFC}, only if \code{add_shrunk_logfc = TRUE}}
#'   \item{\code{zscore}}{Wald statistic}
#'   \item{\code{p_value}}{nominal p-value corresponding to the contrast between the training group 
#'     (e.g., 1-week females) and the sex-matched sedentary controls}
#'   \item{\code{comparison_average_intensity}}{average normalized RNA-seq counts for samples in the training group (e.g., 1-week females)}
#'   \item{\code{comparison_average_intensity_se}}{standard error of \code{comparison_average_intensity}}
#'   \item{\code{reference_average_intensity}}{average normalized RNA-seq counts for sex-matched sedentary control samples}
#'   \item{\code{reference_average_intensity_se}}{standard error of \code{reference_average_intensity}}
#' }
#' 
#' @export 
#' @import data.table
#' @import MotrpacRatTraining6moData
#'
#' @examples
#' # Perform differential analysis for expressed genes in brown adipose tissue with default parameters, 
#' # i.e., outliers and covariates used for the manuscript; calculate both standard and shrunk log fold-changes
#' dea = transcript_timewise_dea("BAT")
#' 
#' # Same as above but don't calculate shrunk log fold-changes
#' dea = transcript_timewise_dea("BAT", add_shrunk_logfc = FALSE)
#' 
#' # Same as the first example but save the [DESeq2::DESeq2()] DESeqResults objects in an RData file 
#' dea = transcript_timewise_dea("BAT", rdata_outfile = "~/test/BAT_TRNSCRIPT_DA.RData", overwrite = TRUE)
#' 
transcript_timewise_dea = function(tissue, 
                                   covariates = c('pct_globin', 'RIN', 'pct_umi_dup', 'median_5_3_bias'), 
                                   outliers = MotrpacRatTraining6moData::OUTLIERS$viallabel[MotrpacRatTraining6moData::OUTLIERS$assay == "TRNSCRPT"],
                                   add_shrunk_logfc = TRUE, 
                                   rdata_outfile = NULL,
                                   overwrite = FALSE,
                                   verbose = FALSE){
  .tissue = tissue # data.table workaround
  
  check_dea_args(.tissue, rdata_outfile, overwrite)
  
  message("Loading data...")
  data = transcript_prep_data(tissue, covariates = covariates, outliers = outliers, center_scale = TRUE, adjust_covariates = TRUE)
  meta = as.data.table(data$metadata)
  counts = data$filt_counts
  outliers = data$outliers
  covariates = data$covariates
  
  sex_res = list()
  for(SEX in unique(meta[,sex])){
    
    message(sprintf("Performing differential expression analysis for %s %ss...", .tissue, SEX))
    
    # subset counts and meta
    curr_samples = meta[sex == SEX, viallabel]
    curr_meta = meta[sex == SEX]
    curr_counts = counts[,curr_samples]
    
    contrasts = list()
    i = 1
    for (tp in c('1w','2w','4w','8w')){
      contrasts[[i]] = c('group', tp, 'control')
      i = i+1
    }

    # standard results 
    deseq_res = run_deseq(curr_counts, # filtered counts
                          curr_meta, # metadata
                          covariates, # covariates
                          outcome_of_interest = 'group', # outcome of interest
                          contrasts = contrasts, # list of contrasts in format c(outcome_of_interest, numerator_level, denominator_level)
                          shrink = FALSE,
                          verbose = verbose)
    
    # shrunk results 
    if(add_shrunk_logfc){
      deseq_res_shrunk = run_deseq(curr_counts, # filtered counts
                                   curr_meta, # metadata
                                   covariates, # covariates
                                   outcome_of_interest = 'group', # outcome of interest
                                   contrasts = contrasts, # list of contrasts in format c(outcome_of_interest, numerator_level, denominator_level)
                                   shrink = TRUE,
                                   verbose = verbose)
    }

    # save DESeq2 RData
    if(!is.null(rdata_outfile)){
      if(add_shrunk_logfc){
        save(deseq_res, deseq_res_shrunk, file=rdata_outfile)
        message(sprintf("'deseq_res', 'deseq_res_shrunk' saved in 'rdata_outfile': %s", rdata_outfile))
      }else{
        save(deseq_res, file=rdata_outfile)
        message(sprintf("'deseq_res' saved in 'rdata_outfile': %s", rdata_outfile))
      }
    }
    
    # collect res
    if(add_shrunk_logfc){
      res_shrunk = data.table(deseq_res_shrunk$res)
      res_nonshrunk = data.table(deseq_res$res)
      setnames(res_shrunk, c("log2FoldChange", "lfcSE"), c("shrunk_logFC","shrunk_logFC_se"))
      setnames(res_nonshrunk, c("log2FoldChange", "lfcSE", "stat"), c("logFC","logFC_se", "zscore"))
      res_shrunk = res_shrunk[,.(gene_id, shrunk_logFC, shrunk_logFC_se, numerator, denominator)]
      res = merge(res_nonshrunk, res_shrunk, by=c("gene_id","numerator","denominator"))
    }else{
      res = data.table(deseq_res$res)
      setnames(res, c("log2FoldChange", "lfcSE", "stat"), c("logFC","logFC_se", "zscore"))
    }

    setnames(res, c("numerator","pvalue","gene_id"), c("comparison_group","p_value","feature_ID"))
    res[,denominator := NULL]
    
    # add some columns
    res[,sex := SEX]
    res[,removed_samples := paste0(curr_outliers, collapse=',')]
    # res[,covariates := paste0(covariates, collapse=',')] added within run_deseq()
    
    # add average intensities 
    norm_counts = as.data.frame(counts(deseq_res$dds, normalized=T))
    ref_sub = norm_counts[,as.character(curr_meta[group == 'control', viallabel])]
    ref_means = rowMeans(ref_sub, na.rm=T)
    ref_se = apply(ref_sub, 1, function(x) sd(x)/sqrt(sum(!is.na(x))) )
    mlist = list()
    i = 1
    for(tp in unique(res[,comparison_group])){
      # get average values
      counts_sub = norm_counts[,as.character(curr_meta[group == tp, viallabel])]
      counts_means = data.table(sex=SEX,
                                comparison_group=tp,
                                comparison_average_intensity=rowMeans(counts_sub, na.rm=T),
                                comparison_average_intensity_se=apply(counts_sub, 1, 
                                                                      function(x) sd(x)/sqrt(sum(!is.na(x))) ),
                                reference_average_intensity=ref_means,
                                reference_average_intensity_se=ref_se,
                                feature_ID=rownames(counts_sub))
      mlist[[i]] = counts_means
      i = i+1
    }
    cmeans = rbindlist(mlist)
    dt = merge(res, cmeans, by=c('feature_ID', 'sex', 'comparison_group'))
    
    sex_res[[SEX]] = dt
  }
  
  dt = rbindlist(sex_res)
  
  # add columns
  dt[,tissue := .tissue]
  dt[,assay_code := 'transcript-rna-seq']
  dt[,assay := 'TRNSCRPT']
  dt[,feature := sprintf("TRNSCRPT;%s;%s", .tissue, feature_ID)]

  if(add_shrunk_logfc){
    dt = dt[,.(
      feature,
      feature_ID,
      sex,
      comparison_group,
      assay,      
      assay_code,
      tissue,
      tissue_code,
      covariates,
      removed_samples,
      logFC,
      logFC_se,
      shrunk_logFC,
      shrunk_logFC_se,
      zscore,
      p_value,
      comparison_average_intensity,
      comparison_average_intensity_se,
      reference_average_intensity,
      reference_average_intensity_se
    )]
  }else{
    dt = dt[,.(
      feature,
      feature_ID,
      sex,
      comparison_group,
      assay,      
      assay_code,
      tissue,
      tissue_code,
      covariates,
      removed_samples,
      logFC,
      logFC_se,
      zscore,
      p_value,
      comparison_average_intensity,
      comparison_average_intensity_se,
      reference_average_intensity,
      reference_average_intensity_se
    )]
  }
  message("Done.")
  return(as.data.frame(dt))
}


#' TODO
#' Title
#'
#' @param tissue character, tissue abbreviation, one of [MotrpacRatTraining6moData::TISSUE_ABBREV]
#' @param covariates character vector of covariates that correspond to column names of [MotrpacRatTraining6moData::TRNSCRPT_META].
#'   Defaults to covariates that were used for the manuscript. 
#' @param outliers vector of viallabels to exclude during differential analysis. Defaults
#'   to \code{MotrpacRatTraining6moData::OUTLIERS$viallabel[MotrpacRatTraining6moData::OUTLIERS$assay == "TRNSCRPT"]}
#' @param rdata_outfile NULL or path in which to save DESeq2 objects in an RData file 
#' @param overwrite boolean, whether to overwrite the file if \code{rdata_outfile} exists
#' @param verbose boolean, whether to print the DESeq2 design string
#'
#' @return TODO
#' 
#' @export
#' @import metap
#' @import DESeq2
#' @import data.table
#' @import MotrpacRatTraining6moData
#'
#' @examples
#' TODO
#' 
transcript_training_dea = function(tissue, 
                                   covariates = c('pct_globin', 'RIN', 'pct_umi_dup', 'median_5_3_bias'), 
                                   outliers = MotrpacRatTraining6moData::OUTLIERS$viallabel[MotrpacRatTraining6moData::OUTLIERS$assay == "TRNSCRPT"],
                                   rdata_outfile = NULL,
                                   overwrite = FALSE,
                                   verbose = FALSE){
  
  .tissue = tissue # data.table workaround
  
  check_dea_args(.tissue, rdata_outfile, overwrite)
  
  message("Loading data...")
  data = transcript_prep_data(tissue, covariates = covariates, outliers = outliers, center_scale = TRUE, adjust_covariates = TRUE)
  meta = as.data.table(data$metadata)
  counts = data$filt_counts
  outliers = data$outliers
  covariates = data$covariates
  
  sex_res = list()
  for(SEX in unique(meta[,sex])){
    
    # subset counts and meta
    curr_samples = meta[sex == SEX, viallabel]
    curr_meta = meta[sex == SEX]
    curr_counts = counts[,curr_samples]
    
    full = paste0('~', paste0(c(curr_cov, 'group'), collapse=' + '))
    reduced = paste0('~', paste0(curr_cov, collapse=' + '))
    
    dds = DESeqDataSetFromMatrix(countData = curr_counts,
                                 colData = curr_meta,
                                 design = eval(parse(text=full)))
    dds = estimateSizeFactors(dds)
    dds = estimateDispersions(dds)
    dds = nbinomLRT(dds, reduced=eval(parse(text=reduced)), maxit=500)
    
    res = results(dds)
    res_dt = data.table(feature_ID = rownames(res), 
                        lrt = res$stat,
                        p_value = res$pvalue,
                        tissue = tissue_code,
                        removed_samples = paste0(curr_outliers, collapse=','),
                        full_model=gsub(' ','',full),
                        reduced_model=gsub(' ','',reduced))
    # add some columns
    res_dt[,feature := sprintf("TRNSCRPT;%s;%s", .tissue, feature_ID)]
    res_dt[,assay := "TRNSCRPT"]
    res_dt[,assay_code := "transcript-rna-seq"]
    res_dt[,tissue := .tissue]
    res_dt[,tissue_code := TISSUE_ABBREV_TO_CODE[[.tissue]]]
    
    res_dt = res_dt[,.(
      feature,
      feature_ID,
      assay,      
      assay_code,
      tissue,
      tissue_code,
      removed_samples,
      lrt,
      p_value,
      full_model,
      reduced_model
    )]
    
    sex_res[[SEX]] = res_dt
  }
  
  if(length(sex_res) > 1){
    male = sex_res[['male']]
    female = sex_res[['female']]
    merged = merge(male, female, by=c('feature_ID','assay','tissue','removed_samples'),
                   suffixes=c("_male","_female"), all=T)
    missing = merged[is.na(p_value_male) | is.na(p_value_female)]
    complete = merged[!is.na(p_value_male) & !is.na(p_value_female)]
    complete[, p_value := sumlog(c(p_value_male, p_value_female))$p, by=seq(1, nrow(complete))]
    missing[, p_value := ifelse(is.na(p_value_male), p_value_female, p_value_male)]
    res_dt = rbindlist(list(complete, missing))
  }else{
    res_dt = sex_res[[1]]
  }
  
  return(as.data.frame(res_dt))
}
