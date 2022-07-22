#' Format covariates for differential analysis
#' 
#' If fewer than 5% of values are missing from a continuous covariate, replace 
#' missing values with the mean in \code{meta}. Otherwise, remove the covariate from \code{covar}. 
#' If values are missing from a factor covariate, remove the covariate from \code{covar}. 
#' If a covariate has constant values, remove it from \code{covar}. 
#' Note that covariates are only removed from \code{covar}, not \code{meta}.
#'
#' @param covar string or character vector of covariate names that correspond to column names of \code{meta}
#' @param meta sample by variable data frame of metadata
#' @param center_scale boolean, whether to center and scale continuous variables
#'
#' @return named list of two items: 
#' \describe{
#'   \item{\code{meta}}{data frame input \code{meta} with covariates imputed and/or centered and scaled as necessary}
#'   \item{\code{covariates}}{character vector input \code{covar} after removing covariates as necessary}
#' }
#' @export
#' @import data.table
#' 
#' @seealso [transcript_prep_data()]
#'
#' @examples
#' meta = data.frame(V1 = c(rnorm(19), NA),
#'                   V2 = c(rnorm(12), rep(NA, 8)))
#' covar = c("V1","V2")
#' result = fix_covariates(covar, meta)
#' result = fix_covariates(covar, meta, center_scale = TRUE)
#' 
fix_covariates = function(covar, meta, center_scale = FALSE){
  meta = data.table(meta)
  # make sure there are no missing values in covariates
  for(cov in covar){
    if(any(is.na(meta[,get(cov)]))){
      num_missing = nrow(meta[is.na(get(cov))])
      if(is.character(meta[,get(cov)]) | is.factor(meta[,get(cov)])){
        # don't know how to handle missing values for factor. remove.
        warning(sprintf("Categorical variable of interest %s has %s missing values. Removing.\n", cov, num_missing))
        covar = covar[covar != cov]
      }else{
        # if continuous 
        # remove if missing for >5% of samples
        if(num_missing/nrow(meta) > 0.05){
          warning(sprintf("Numeric variable of interest %s has %s missing values. Removing.\n", cov, num_missing))
          covar = covar[covar != cov]
        }else{
          warning(sprintf("Numeric variable of interest %s has %s missing values. Replacing missing values with mean.\n", cov, num_missing))
          meta[is.na(get(cov)), (cov) := mean(meta[,get(cov)], na.rm=T)]
        }
      }
    }
  }
  
  # center and scale continuous variables 
  for (cov in covar){
    # remove if constant
    if(length(unique(meta[,get(cov)])) == 1){
      message(sprintf("Covariate %s is constant. Removing.", cov))
      covar = covar[covar != cov]
    }else{
      # center and scale
      if(is.numeric(meta[,get(cov)])){
        meta[,(cov) := scale(meta[,get(cov)], center = center_scale, scale = center_scale)]
      }
    }
  }
  
  meta = as.data.frame(meta)
  return(list(meta=meta, covariates=covar))
}


#' Preprocess RNA-seq data 
#' 
#' Collect filtered raw counts, normalized sample-level data, phenotypic data, RNA-seq metadata, 
#' covariates, and outliers associated with a given tissue. 
#'
#' @param tissue character, tissue abbreviation, one of [TISSUE_ABBREV]
#' @param sex NULL to return data from both sexes or one of "male" or "female" to return data from the specified sex
#' @param covariates character vector of covariates that correspond to column names of [TRNSCRPT_META].
#'   Defaults to covariates that were used for the manuscript. 
#' @param outliers vector of viallabels to exclude from the returned data. Defaults
#'   to \code{OUTLIERS$viallabel[OUTLIERS$assay == "TRNSCRPT"]}
#' @param adjust_covariates boolean, whether to adjust covariates using [fix_covariates()]. 
#'   Only applies if \code{covariates} is not NULL. 
#' @param center_scale boolean, whether to center and scale continuous covariates within [fix_covariates()]. 
#'   Only applies if \code{adjust_covariates} is also TRUE. 
#'
#' @return named list of five items: 
#' \describe{
#'   \item{\code{metadata}}{data frame of combined [PHENO] and [TRNSCRPT_META], filtered to samples in \code{tissue}.
#'     If \code{adjust_covariates = TRUE}, missing values in \code{covariates} are imputed. 
#'     If also \code{center_scale = TRUE}, continuous variables named by \code{covariates} are centered and scaled.}
#'   \item{\code{covariates}}{character vector of covariates to adjust for during differential analysis. For all tissues except VENACV,
#'     this vector is a (sub)set of the input list of covariates. Covariates are removed from this vector if there are too
#'     many missing values or if all values are constant. See [fix_covariates()] for more details.
#'     If \code{tissue = "VENACV"}, the Ensembl ID for Ucp1 is also added as a covariate.}
#'   \item{\code{counts}}{data frame of raw counts with Ensembl IDs (which are also TRNSCRPT \code{feature_ID}s) 
#'     as row names and vial labels as column names. See [TRNSCRPT_RAW_COUNTS] for details.}
#'   \item{\code{norm_data}}{data frame of TMM-normalized data with Ensembl IDs (which are also TRNSCRPT 
#'     \code{feature_ID}s) as row names and vial labels as column names. See [TRNSCRPT_NORM_DATA] for details.}
#'   \item{\code{outliers}}{subset of \code{outliers} in input removed from the data}
#' }
#' 
#' @seealso [OUTLIERS], [fix_covariates()], [PHENO], [TRNSCRPT_META], [TRNSCRPT_RAW_COUNTS], [TRNSCRPT_NORM_DATA]
#' 
#' @export
#' @import data.table
#' @import MotrpacRatTraining6moData
#'
#' @examples
#' # Process gastrocnemius RNA-seq data with default parameters, i.e., return data from both 
#' # sexes, remove established outliers, impute missing values in default covariates 
#' gastroc_data1 = transcript_prep_data("SKM-GN")
#' 
#' # Same as above but do not remove outliers if they exist 
#' gastroc_data2 = transcript_prep_data("SKM-GN", outliers = NULL)
#' 
#' # Same as above but do not adjust existing variables in the metadata  
#' gastroc_data3 = transcript_prep_data("SKM-GN", covariates = NULL, outliers = NULL)
#' 
#' # Same as above but only return data from male samples
#' gastroc_data4 = transcript_prep_data("SKM-GN", covariates = NULL, outliers = NULL, sex = "male")
#' 
#' # Same as gastroc_data2 but also center and scale default continuous covariates in the returned metadata,
#' # which is also done within [run_deseq()] (called by [transcript_timewise_dea()]) 
#' gastroc_data4 = transcript_prep_data("SKM-GN", outliers = NULL, center_scale = TRUE)
#' 
transcript_prep_data = function(tissue, 
                                sex = NULL, 
                                covariates = c('pct_globin', 'RIN', 'pct_umi_dup', 'median_5_3_bias'), 
                                outliers = na.omit(MotrpacRatTraining6moData::OUTLIERS$viallabel[MotrpacRatTraining6moData::OUTLIERS$assay == "TRNSCRPT"]),
                                adjust_covariates = TRUE,
                                center_scale = FALSE){
  # data.table workaround
  .tissue = tissue
  .sex = sex
  
  if(!is.null(.sex)){
    if(!.sex %in% c('male','female',NULL)){
      stop("'sex' must be one of the following values:\n  'male', 'female', NULL")
    }
  }
  
  accepted_tissues = TISSUE_ABBREV[!TISSUE_ABBREV == "PLASMA"]
  if(!tissue %in% accepted_tissues){
    warning(sprintf("'%s' is not an accepted tissue abbreviation. TRNSCRPT data is available for the following tissues:\n  %s",
                    .tissue,
                    paste0(accepted_tissues, collapse=", ")))
    return()
  }
  
  # load data
  counts = get(sprintf("TRNSCRPT_%s_RAW_COUNTS", gsub("-","",.tissue)))
  tmm = get(sprintf("TRNSCRPT_%s_NORM_DATA", gsub("-","",.tissue)))
  rownames(tmm) = tmm$feature_ID
  tmm$feature_ID = NULL
  
  # filter counts by genes in normalized data 
  counts = counts[rownames(tmm),]
  
  # format metadata
  pheno = as.data.table(PHENO)
  meta = as.data.table(TRNSCRPT_META)
  meta[,viallabel := as.character(vial_label)] 
  meta = merge(pheno, meta, by="viallabel")
  
  # filter by tissue
  meta = meta[tissue == .tissue]
  
  # filter by sex
  if(!is.null(.sex)){
    meta = meta[sex == .sex]
  }

  # remove outliers 
  outliers = outliers[!is.na(outliers)]
  curr_outliers = c()
  if(!is.null(outliers)){
    outliers = as.character(outliers)
    curr_outliers = outliers[outliers %in% meta[,viallabel]]
    meta = meta[!viallabel %in% outliers]
  }
  
  # subset columns
  counts = counts[meta[,viallabel]]
  tmm = tmm[meta[,viallabel]]
  
  if(.tissue == 'VENACV'){
    # add Ucp1 as a covariate
    ucp1 = data.table(viallabel = colnames(counts), ucp1 = unname(unlist(counts['ENSRNOG00000003580',])))
    meta = merge(meta, ucp1, by = 'viallabel')
    covariates = c(covariates, 'ucp1')
  }
  
  # impute missing values and filter covariates
  if(adjust_covariates & !is.null(covariates)){
    new = fix_covariates(covariates, meta, center_scale)
    covariates = new$covariates
    meta = data.table(new$meta)
  }
  
  meta[,sex_group := paste0(sex, ';', group)]
  
  return(list(metadata = meta, 
              covariates = covariates, 
              filt_counts = counts,
              norm_data = tmm, 
              outliers = curr_outliers))
}


#' TODO
#' Prepare ATAC-seq dataset
#' 
#' Retrieve and format ATAC-seq sample-level data and metadata for a given tissue. 
#'
#' @param TISSUE_CODE 
#' @param outliers 
#' @param gsutil_path 
#' @param scratch 
#' @param parallel 
#'
#' @return
#' @export
#' @import data.table
#'
#' @examples
#' TODO
#' 
preprocess_pass1b_atacseq_gcp = function(TISSUE_CODE, # or dataset, if it starts with tissue_code
                                         outliers = NULL, 
                                         gsutil_path='~/google-cloud-sdk/bin/gsutil', 
                                         scratch='/oak/stanford/groups/smontgom/nicolerg/tmp', 
                                         parallel=T){
  
  # ATAC-seq meta 
  wet = dl_read_gcp('gs://motrpac-data-freeze-pass/pass1b-06/v1.0/results/epigenomics/qa-qc/motrpac_pass1b-06_epigen-atac-seq_qa-qc-metrics.csv', sep=',', check_first = T)
  # DMAQC metadata 
  dmaqc_meta = dl_format_pheno(scratch, gsutil_path, parallel=parallel)
  wet[,viallabel := as.character(viallabel)]
  meta = merge(dmaqc_meta, wet, by='viallabel', all.y=T)
  
  colnames(meta) = tolower(colnames(meta))
  
  tissue_number = toupper(gsub("-.*","",TISSUE_CODE))
  tissue = gsub(",.*","",TISSUE_CODE)
  # load raw counts 
  counts = dl_read_gcp(sprintf("gs://motrpac-data-freeze-pass/pass1b-06/v1.0/results/epigenomics/%s/epigen-atac-seq/motrpac_pass1b-06_%s_epigen-atac-seq_counts.txt.gz",TISSUE_CODE,TISSUE_CODE), check_first = T)
  counts = data.frame(counts, check.names=F)
  rownames(counts) = paste0(counts$chrom, ':', counts$start, '-', counts$end)
  counts[,c('chrom','start','end')] = NULL
  
  # remove outliers
  if(!is.null(outliers)){
    counts[,as.character(outliers)] = NULL
  }
  
  # separate data by site
  data_list = list()
  curr_meta = meta[viallabel %in% colnames(counts) & !grepl("^8", viallabel)]
  for(site in unique(curr_meta[,get_site])){
    
    site_meta = curr_meta[get_site==site]
    site_counts = counts[,site_meta[,viallabel]]
    
    # remove non-auto peaks
    site_counts = site_counts[grepl("^chr[0-9]|^chrY|^chrX", rownames(site_counts)),]
    
    # exclude low count peaks in the current dataset
    # at least 10 counts in N samples
    n_samples = 4
    filt_counts = site_counts[rowSums(data.frame(lapply(site_counts, function(x) as.numeric(x >= 10)), check.names=F)) >= n_samples,]
    
    # quantile normalize
    # this takes a couple of minutes given the size of the peak x sample counts matrix
    norm_counts = voom(filt_counts,normalize.method = "quantile")
    sub_norm = round(norm_counts$E,2)
    
    label = sprintf("%s,%s", tissue, tolower(site))
    data_list[[label]] = list(meta = site_meta, 
                              #voom_obj = norm_counts, # this might break some old code, but voom weights need design
                              norm = sub_norm,
                              unfilt_counts = site_counts,
                              filt_counts = filt_counts)
  }
  return(data_list)
}

# TODO
preprocess_pass1b_atacseq_gcp = function(TISSUE_CODE, # or dataset, if it starts with tissue_code
                                         outliers = NULL, 
                                         gsutil_path='~/google-cloud-sdk/bin/gsutil', 
                                         scratch='/oak/stanford/groups/smontgom/nicolerg/tmp', 
                                         parallel=T){
  
  # ATAC-seq meta 
  wet = dl_read_gcp('gs://motrpac-data-freeze-pass/pass1b-06/v1.0/results/epigenomics/qa-qc/motrpac_pass1b-06_epigen-atac-seq_qa-qc-metrics.csv', sep=',', check_first = T)
  # DMAQC metadata 
  dmaqc_meta = dl_format_pheno(scratch, gsutil_path, parallel=parallel)
  wet[,viallabel := as.character(viallabel)]
  meta = merge(dmaqc_meta, wet, by='viallabel', all.y=T)
  
  colnames(meta) = tolower(colnames(meta))
  
  tissue_number = toupper(gsub("-.*","",TISSUE_CODE))
  tissue = gsub(",.*","",TISSUE_CODE)
  # load raw counts 
  counts = dl_read_gcp(sprintf("gs://motrpac-data-freeze-pass/pass1b-06/v1.0/results/epigenomics/%s/epigen-atac-seq/motrpac_pass1b-06_%s_epigen-atac-seq_counts.txt.gz",TISSUE_CODE,TISSUE_CODE), check_first = T)
  counts = data.frame(counts, check.names=F)
  rownames(counts) = paste0(counts$chrom, ':', counts$start, '-', counts$end)
  counts[,c('chrom','start','end')] = NULL
  
  # remove outliers
  if(!is.null(outliers)){
    counts[,as.character(outliers)] = NULL
  }
  
  # separate data by site
  data_list = list()
  curr_meta = meta[viallabel %in% colnames(counts) & !grepl("^8", viallabel)]
  for(site in unique(curr_meta[,get_site])){
    
    site_meta = curr_meta[get_site==site]
    site_counts = counts[,site_meta[,viallabel]]
    
    # remove non-auto peaks
    site_counts = site_counts[grepl("^chr[0-9]|^chrY|^chrX", rownames(site_counts)),]
    
    # exclude low count peaks in the current dataset
    # at least 10 counts in N samples
    n_samples = 4
    filt_counts = site_counts[rowSums(data.frame(lapply(site_counts, function(x) as.numeric(x >= 10)), check.names=F)) >= n_samples,]
    
    # quantile normalize
    # this takes a couple of minutes given the size of the peak x sample counts matrix
    norm_counts = voom(filt_counts,normalize.method = "quantile")
    sub_norm = round(norm_counts$E,2)
    
    label = sprintf("%s,%s", tissue, tolower(site))
    data_list[[label]] = list(meta = site_meta, 
                              #voom_obj = norm_counts, # this might break some old code, but voom weights need design
                              norm = sub_norm,
                              unfilt_counts = site_counts,
                              filt_counts = filt_counts)
  }
  return(data_list)
}

