#' Preprocess RNA-seq data 
#' 
#' Collect filtered raw counts, normalized sample-level data, phenotypic data, RNA-seq metadata, 
#' covariates, and outliers associated with a given tissue. 
#'
#' @param tissue`r tissue()`
#' @param sex `r sex()`
#' @param covariates character vector of covariates that correspond to column names of [MotrpacRatTraining6moData::TRNSCRPT_META].
#'   Defaults to covariates that were used for the manuscript. 
#' @param outliers vector of viallabels to exclude from the returned data. Defaults
#'   to \code{[MotrpacRatTraining6moData::OUTLIERS]$viallabel}
#' @param adjust_covariates boolean, whether to adjust covariates using [fix_covariates()]. 
#'   Only applies if \code{covariates} is not NULL. 
#' @param center_scale boolean, whether to center and scale continuous covariates within [fix_covariates()]. 
#'   Only applies if \code{adjust_covariates} is also TRUE. 
#'
#' @return named list of five items: 
#' \describe{
#'   \item{\code{metadata}}{data frame of combined [MotrpacRatTraining6moData::PHENO] and [MotrpacRatTraining6moData::TRNSCRPT_META], filtered to samples in \code{tissue}.
#'     If \code{adjust_covariates = TRUE}, missing values in \code{covariates} are imputed. 
#'     If also \code{center_scale = TRUE}, continuous variables named by \code{covariates} are centered and scaled.}
#'   \item{\code{covariates}}{character vector of covariates to adjust for during differential analysis. For all tissues except VENACV,
#'     this vector is a (sub)set of the input list of covariates. Covariates are removed from this vector if there are too
#'     many missing values or if all values are constant. See [fix_covariates()] for more details.
#'     If \code{tissue = "VENACV"}, the Ensembl ID for Ucp1 is also added as a covariate.}
#'   \item{\code{counts}}{data frame of raw counts with Ensembl IDs (which are also TRNSCRPT \code{feature_ID}s) 
#'     as row names and vial labels as column names. See [MotrpacRatTraining6moData::TRNSCRPT_RAW_COUNTS] for details.}
#'   \item{\code{norm_data}}{data frame of TMM-normalized data with Ensembl IDs (which are also TRNSCRPT 
#'     \code{feature_ID}s) as row names and vial labels as column names. See [MotrpacRatTraining6moData::TRNSCRPT_NORM_DATA] for details.}
#'   \item{\code{outliers}}{subset of \code{outliers} in input removed from the data}
#' }
#' 
#' @seealso [MotrpacRatTraining6moData::OUTLIERS], 
#'   [fix_covariates()], 
#'   [MotrpacRatTraining6moData::PHENO], 
#'   [MotrpacRatTraining6moData::TRNSCRPT_META], 
#'   [MotrpacRatTraining6moData::TRNSCRPT_RAW_COUNTS], 
#'   [MotrpacRatTraining6moData::TRNSCRPT_NORM_DATA]
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
                                outliers = na.omit(MotrpacRatTraining6moData::OUTLIERS$viallabel),
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
  tmm[,c("feature","feature_ID","tissue","assay")] = NULL
  
  # filter counts by genes in normalized data 
  rownames(counts) = counts$feature_ID
  counts[,c("feature","feature_ID","tissue","assay")] = NULL
  counts = counts[rownames(counts),]
  
  # format metadata
  pheno = as.data.table(PHENO)
  meta = as.data.table(TRNSCRPT_META)
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


#' Prepare ATAC-seq dataset
#' 
#' Retrieve and format ATAC-seq sample-level data and metadata for a given tissue. 
#'
#' @param tissue `r tissue()`
#' @param scratchdir character, local directory in which to download data from 
#'   Google Cloud Storage. Current working directory by default. 
#' @param sex `r sex()`
#' @param covariates character vector of covariates that correspond to column names of [MotrpacRatTraining6moData::ATAC_META].
#'   Defaults to covariates that were used for the manuscript. 
#' @param outliers vector of viallabels to exclude from the returned data. Defaults
#'   to \code{[MotrpacRatTraining6moData::OUTLIERS]$viallabel}
#' @param nrows integer, number of rows to return. Defaults to Inf. Useful to return a subset of a large data frame for tests. 
#' @param filter_counts bool, whether to return filtered raw counts
#' @param return_normalized_data bool, whether to also return normalized data 
#'
#' @return named list of five items: 
#' \describe{
#'   \item{\code{metadata}}{data frame of combined [MotrpacRatTraining6moData::PHENO] and 
#'     [MotrpacRatTraining6moData::ATAC_META], filtered to samples in \code{tissue}.}
#'   \item{\code{covariates}}{character vector of covariates to adjust for during differential analysis, same as input}
#'   \item{\code{raw_counts}}{data frame of raw counts with feature IDs
#'     as row names and vial labels as column names. See [MotrpacRatTraining6moData::ATAC_RAW_COUNTS] for details.}
#'   \item{\code{norm_data}}{data frame of quantile-normalized data with feature IDs as row names and vial labels as column names. 
#'     See [MotrpacRatTraining6moData::ATAC_NORM_DATA] for details.}
#'   \item{\code{outliers}}{subset of \code{outliers} in input removed from the data}
#' }
#' @export
#' @import data.table
#' @import MotrpacRatTraining6moData
#'
#' @examples
#' # Process gastrocnemius ATAC-seq data with default parameters, i.e., return data from both 
#' # sexes, remove established outliers, download data to current working directory
#' gastroc_data1 = atac_prep_data("SKM-GN")
#' 
#' # Same as above but do not remove outliers if they exist 
#' gastroc_data2 = atac_prep_data("SKM-GN", outliers = NULL)
#' 
#' # Same as above but only return data from male samples
#' gastroc_data3 = atac_prep_data("SKM-GN", outliers = NULL, sex = "male")
#' 
atac_prep_data = function(tissue, 
                          sex = NULL, 
                          covariates = c("Sample_batch", "peak_enrich.frac_reads_in_peaks.macs2.frip"), 
                          filter_counts = FALSE,
                          return_normalized_data = FALSE, 
                          scratchdir = ".", 
                          outliers = na.omit(MotrpacRatTraining6moData::OUTLIERS$viallabel[MotrpacRatTraining6moData::OUTLIERS$assay == "ATAC"]),
                          nrows = Inf){
  
  # data.table workaround
  .tissue = tissue
  .sex = sex
  
  if(!is.null(.sex)){
    if(!.sex %in% c('male','female',NULL)){
      stop("'sex' must be one of the following values:\n  'male', 'female', NULL")
    }
  }
  
  # ATAC-seq meta 
  wet = data.table(MotrpacRatTraining6moData::ATAC_META)
  
  # DMAQC metadata 
  dmaqc_meta = data.table(MotrpacRatTraining6moData::PHENO)
  meta = merge(dmaqc_meta, wet, by='viallabel', all.y=T)
  
  # load raw counts 
  counts = load_sample_data(tissue = tissue, 
                            assay = "ATAC", 
                            normalized = FALSE, 
                            training_regulated_only = FALSE, 
                            scratchdir = scratchdir,
                            exclude_outliers = FALSE,
                            nrows = nrows)
  rownames(counts) = counts$feature_ID
  counts[,c("feature","feature_ID","tissue","assay")] = NULL
  
  # filter by tissue
  meta = meta[tissue == .tissue]
  
  # filter by sex
  if(!is.null(.sex)){
    meta = meta[sex == .sex]
  }
  
  # remove other tissues and ref stds
  meta = meta[viallabel %in% colnames(counts) & !grepl("^8", viallabel)]
  
  # remove outliers
  curr_outliers = c()
  if(!is.null(outliers)){
    curr_outliers = curr_outliers[curr_outliers %in% meta[,viallabel]]
    meta = meta[!viallabel %in% as.character(curr_outliers)]
  }
  counts = counts[,meta[,viallabel]]
  
  meta = as.data.frame(meta)
  rownames(meta) = meta$viallabel
  
  # normalized data
  norm = NULL
  if(return_normalized_data){
    # load normalized data 
    norm = load_sample_data(tissue = tissue, 
                            assay = "ATAC", 
                            normalized = TRUE, 
                            training_regulated_only = FALSE, 
                            scratchdir = scratchdir,
                            exclude_outliers = FALSE,
                            nrows = nrows)
    rownames(norm) = norm$feature_ID
    norm[,c("feature","feature_ID","tissue","assay")] = NULL
    norm = norm[,meta[,viallabel]]
  }
  
  # filter counts
  if(filter_counts){
    if(is.null(norm)){
      n_samples = 4
      min_count = 10
      # remove non-auto peaks
      counts = counts[grepl("^chr[0-9]|^chrY|^chrX", rownames(counts)),]
      # exclude low count peaks in the current dataset
      # at least min_count counts in n_samples samples
      counts = counts[rowSums(data.frame(lapply(counts, function(x) as.numeric(x >= min_count)), check.names=F)) >= n_samples,]
    }else{
      counts = counts[rownames(norm),]
    }
  }
  
  return(list(metadata = meta, 
              covariates = covariates, 
              counts = counts,
              norm_data = norm, 
              outliers = curr_outliers))
}


#' Load sample-level data
#' 
#' Load raw counts or normalized sample-level data. Optionally filter by training-regulated features. 
#' For epigenetic data (ATAC and METHYL), if \code{training_regulated_only = FALSE}, sample-level data
#' is downloaded from Google Cloud Storage. 
#'
#' @param tissue `r tissue()`
#' @param assay `r assay()`
#' @param normalized bool, whether to return normalized data. If \code{FALSE}, return raw counts. 
#' @param training_regulated_only bool, whether to filter features down to those training-regulated at 5% FDR
#' @param exclude_outliers bool, whether to remove sample outliers specified by [MotrpacRatTraining6moData::OUTLIERS]
#' @param scratchdir character, local directory in which to download data from 
#'   Google Cloud Storage. Current working directory by default. Only applies if \code{assay} is ATAC or METHYL. 
#' @param nrows integer, number of rows to return. Defaults to Inf. Useful to return a subset of a large data frame for tests. 
#'
#' @return a data.frame where features are in rows and numeric columns correspond to sample identifiers (vial labels)
#' @export
#' @import data.table
#' @import MotrpacRatTraining6moData
#'
#' @examples
#' print("TODO")
load_sample_data = function(tissue, 
                            assay, 
                            normalized = TRUE, 
                            training_regulated_only = FALSE, 
                            exclude_outliers = FALSE, 
                            scratchdir = NULL, 
                            nrows = Inf){
  
  # check inputs 
  if(!tissue %in% MotrpacRatTraining6moData::TISSUE_ABBREV){
    stop(sprintf("'tissue' must be one of TISSUE_ABBREV: \n %s", paste(MotrpacRatTraining6moData::TISSUE_ABBREV, collapse=", ")))
  }
  if(!assay %in% MotrpacRatTraining6moData::ASSAY_ABBREV){
    stop(sprintf("'assay' must be one of ASSAY_ABBREV: \n %s", paste(MotrpacRatTraining6moData::ASSAY_ABBREV, collapse=", ")))
  }
  
  if(normalized){
    # normalized data
    if(assay %in% c("METHYL","ATAC")){
      if(training_regulated_only){
        obj_name = sprintf("%s_%s_NORM_DATA_05FDR", assay, gsub("-","",tissue))
        data = get(obj_name)
      }else{
        # download from GCS
        data = get_rdata_from_url(tissue=tissue, assay=assay, suffix="NORM_DATA", scratchdir=scratchdir, nrows=nrows)
      }
    }else{
      obj_name = sprintf("%s_%s_NORM_DATA", assay, gsub("-","",tissue))
      data = get(obj_name)
    }
  }else{
    # raw data 
    if(assay %in% c("METHYL","ATAC","TRNSCRPT")){
      if(assay %in% c("METHYL","ATAC")){
        # download from GCS
        data = get_rdata_from_url(tissue=tissue, assay=assay, suffix="RAW_COUNTS", scratchdir=scratchdir, nrows=nrows)
      }else{
        obj_name = sprintf("%s_%s_RAW_COUNTS", assay, gsub("-","",tissue))
        data = get(obj_name)
      }
    }else{
      stop(sprintf("Non-normalized data not available for %s. Set 'normalized' to TRUE to get normalized sample-level data.", assay))
    }
  }
  
  # apply optional filters
  if(training_regulated_only){
    data = data[!is.na(data$feature),]
  }
  if(nrows < Inf){
    if(nrows > nrow(data)){
      data = data[1:nrows,]
    }
  }
  if(exclude_outliers){
    outliers = data.table(MotrpacRatTraining6moData::OUTLIERS)
    .assay = assay
    .tissue = .tissue
    curr_outliers = na.omit(outliers[assay == .assay & tissue == .tissue, viallabel])
    if(.tissue == "VENACV"){
      curr_outliers = c(curr_outliers, na.omit(outliers[tissue == .tissue, viallabel]))
    }
    if(length(curr_outliers)>0){
      data[,curr_outliers] = NULL
    }
  }
  
  return(data)
}


#' Load raw METHYL data
#'
#' @param tissue `r tissue()`
#' @param scratchdir character, local directory in which to download data from 
#'   Google Cloud Storage. Current working directory by default.
#'
#' @return named list of three items: 
#' \describe{
#'   \item{\code{counts}}{data frame of raw counts with features in rows and samples (unmethylated or methylated) in columns}
#'   \item{\code{samples}}{data frame of sample-level metadata, including \code{group}, \code{lib.size}, and \code{norm.factors}}
#'   \item{\code{genes}}{data frame of feature-level gene information, including \code{Chr} and \code{Locus}}
#' }
#' @export
#'
#' @examples
#' # Load raw METHYL data for gastrocnemius 
#' data = load_methyl_raw_data("SKM-GN", "/tmp")
#' 
load_methyl_raw_data = function(tissue, scratchdir = "."){
  url = sprintf("https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/%s_raw.RData", tissue)
  data = get_rdata_from_url(url = url, scratchdir = scratchdir)
  return(data)
}


#' Load METHYL feature annotation 
#'
#' @param scratchdir character, local directory in which to download data from 
#'   Google Cloud Storage. Current working directory by default.
#'
#' @return data frame. See [MotrpacRatTraining6moData::METHYL_FEATURE_ANNOT] for details. 
#' @export
#'
#' @examples
#' feature_annot = load_methyl_feature_annotation("/tmp")
#' 
#' @source <https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_FEATURE_ANNOT.rda> 
load_methyl_feature_annotation = function(scratchdir = "."){
  fa = get_rdata_from_url(url = "https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_FEATURE_ANNOT.rda",
                   scratchdir = scratchdir)
  return(fa)
}


#' Load RData from GCS
#' 
#' Function to download RData from the web. Intended to download
#' epigenetic data from Google Cloud Storage but can be used for any public URL
#' that points to an RData file.
#'
#' @param tissue `r tissue()`. Only used if \code{url} is NULL. 
#' @param assay `r assay()`. Only used if \code{url} is NULL. 
#' @param suffix character, object suffix. Only used if \code{url} is NULL. 
#' @param scratchdir character, local directory in which to download data from 
#'   the web. Current working directory by default.
#' @param url character, RData URL. Optional if URL cannot be constructed from 
#'   \code{tissue}, \code{assay}, and \code{suffix}. 
#' @param obj_name character, name of object saved within RData file. Only required
#'   if more than one object is saved in the RData file. 
#' @param nrows integer, number of rows to return. Defaults to Inf. 
#'
#' @return object saved in the RData file, optionally specified by \code{obj_name} 
#' 
#' @export
#' 
#' @examples
#' # return gastrocnemius chromatin accessibility differential analysis results 
#' data = get_rdata_from_url(tissue="SKM-GN", assay="ATAC", suffix="DA", scratchdir="/tmp")
#' 
#' # return brown adipose DNA methylation raw counts 
#' data = get_rdata_from_url(tissue="BAT", assay="METHYL", suffix="RAW_COUNTS", scratchdir="/tmp")
#' 
#' # return DNA methylation feature annotation
#' data = get_rdata_from_url(url="https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_FEATURE_ANNOT.rda", scratchdir="/tmp")
#' 
#' # return raw DNA methylation data for brown adipose
#' data = get_rdata_from_url(url="https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/BAT_raw.RData", scratchdir="/tmp")
#' 
get_rdata_from_url = function(tissue=NULL, assay=NULL, suffix=NULL, scratchdir=".", url=NULL, obj_name=NULL, nrows=Inf){
  # set option if default is low
  if(getOption("timeout") < 1e3){
    options(timeout=1e3)
  }

  # if dir doesn't exist, try to make it
  if(!dir.exists(scratchdir)){
    dir.create(scratchdir, recursive = TRUE)
  }
  
  if(!is.null(assay) & !assay %in% c("ATAC","METHYL")){
    warning("Only ATAC and METHYL data are currently available in GCS.")
  }
  
  possible_tissues = c(
    'BAT',
    'HEART',
    'HIPPOC',
    'KIDNEY',
    'LIVER',
    'LUNG',
    'SKM-GN',
    'WAT-SC'
  )
  if(!is.null(tissue) & !tissue %in% possible_tissues){
    warning(sprintf("Epigenetic data are only available for the following tissues:\n %s", paste(possible_tissues, collapse=", ")))
  }
  
  if(is.null(url)){
    if(is.null(tissue) | is.null(assay) | is.null(suffix)){
      stop("If a URL is not provided, 'tissue', 'assay', and 'suffix' must be non-null.")
    }
    gcs_bucket = "https://storage.googleapis.com/motrpac-rat-training-6mo-extdata"
    url = sprintf("%s/%s_%s_%s.rda", gcs_bucket, assay, gsub("-","",tissue), suffix)
  }
  
  local = sprintf("%s/%s", scratchdir, basename(url))
  download.file(url, local)
  
  if(!is.null(obj_name)){
    load(local)
    data = get(obj_name)
  }else{
    data = get(load(local))
  }
  
  if(nrows < Inf & (is.data.frame(data) | is.matrix(data))){
    data = data[1:nrows,]
  }

  return(data)
}
