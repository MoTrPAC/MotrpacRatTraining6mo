#' Preprocess RNA-seq data 
#' 
#' Collect filtered raw counts, normalized sample-level data, phenotypic data, RNA-seq metadata, 
#' covariates, and outliers associated with a given tissue. 
#'
#' @param tissue `r tissue()`
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
#' # Same as gastroc_data2 but also center and scale default continuous covariates 
#' # in the returned metadata, which is also done within [run_deseq()] 
#' # (called by [transcript_timewise_dea()]) 
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
  
  accepted_tissues = MotrpacRatTraining6moData::TISSUE_ABBREV[!MotrpacRatTraining6moData::TISSUE_ABBREV == "PLASMA"]
  if(!tissue %in% accepted_tissues){
    warning(sprintf("'%s' is not an accepted tissue abbreviation. TRNSCRPT data is available for the following tissues:\n  %s",
                    .tissue,
                    paste0(accepted_tissues, collapse=", ")))
    return()
  }
  
  # load data
  counts = load_sample_data(.tissue, "TRNSCRPT", normalized = FALSE)
  #counts = get(sprintf("TRNSCRPT_%s_RAW_COUNTS", gsub("-","",.tissue)))
  tmm = load_sample_data(.tissue, "TRNSCRPT", normalized = TRUE)
  #tmm = get(sprintf("TRNSCRPT_%s_NORM_DATA", gsub("-","",.tissue)))
  rownames(tmm) = tmm$feature_ID
  tmm[,c("feature","feature_ID","tissue","assay")] = NULL
  
  # filter counts by genes in normalized data 
  rownames(counts) = counts$feature_ID
  counts[,c("feature","feature_ID","tissue","assay")] = NULL
  counts = counts[rownames(tmm),]
  
  # format metadata
  pheno = data.table::as.data.table(MotrpacRatTraining6moData::PHENO)
  meta = data.table::as.data.table(MotrpacRatTraining6moData::TRNSCRPT_META)
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
    ucp1 = data.table::data.table(viallabel = colnames(counts), ucp1 = unname(unlist(counts['ENSRNOG00000003580',])))
    meta = merge(meta, ucp1, by = 'viallabel')
    covariates = c(covariates, 'ucp1')
  }
  
  # impute missing values and filter covariates
  if(adjust_covariates & !is.null(covariates)){
    new = fix_covariates(covariates, meta, center_scale)
    covariates = new$covariates
    meta = data.table::data.table(new$meta)
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
          meta[is.na(get(cov)), (cov) := mean(meta[,get(cov)], na.rm=TRUE)]
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
#' @param tissue character, tissue abbreviation, one of "BAT", "HEART", "HIPPOC", "KIDNEY", "LIVER", "LUNG", "SKM-GN", "WAT-SC" 
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
#'
#' @examples
#' \dontrun{
#' # Process gastrocnemius ATAC-seq data with default parameters, i.e., return data from both 
#' # sexes, remove established outliers, download data to current working directory
#' gastroc_data1 = atac_prep_data("SKM-GN")
#' 
#' # Same as above but do not remove outliers if they exist 
#' gastroc_data2 = atac_prep_data("SKM-GN", outliers = NULL)
#' }
#' # Same as above but only return data from male samples
#' gastroc_data3 = atac_prep_data("SKM-GN", outliers = NULL, sex = "male")
#' 
atac_prep_data = function(tissue, 
                          sex = NULL, 
                          covariates = c("Sample_batch", "peak_enrich.frac_reads_in_peaks.macs2.frip"), 
                          filter_counts = FALSE,
                          return_normalized_data = FALSE, 
                          scratchdir = ".", 
                          outliers = data.table::data.table(
                            MotrpacRatTraining6moData::OUTLIERS)[assay=="ATAC",viallabel],
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
  wet = data.table::data.table(MotrpacRatTraining6moData::ATAC_META)
  
  # DMAQC metadata 
  dmaqc_meta = data.table::data.table(MotrpacRatTraining6moData::PHENO)
  meta = merge(dmaqc_meta, wet, by='viallabel', all.y=TRUE)
  
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
    curr_outliers = outliers[outliers %in% meta[,viallabel]]
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
      counts = counts[rowSums(data.frame(lapply(counts, function(x) as.numeric(x >= min_count)), check.names=FALSE)) >= n_samples,]
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
#' @param warnings bool, whether to print warnings to the console. \code{TRUE} by default. 
#'
#' @return a data.frame where features are in rows and numeric columns correspond to sample identifiers (vial labels)
#' 
#' @seealso [combine_normalized_data()]
#' 
#' @export
#'
#' @examples
#' # Load RNA-seq raw counts for liver
#' data = load_sample_data("LIVER", "TRNSCRPT", normalized = FALSE)
#' 
#' # Load normalized metabolomics data for gastrocnemius
#' data = load_sample_data("SKM-GN", "METAB")
#' 
#' # Load normalized protein abundance data for heart
#' data = load_sample_data("HEART", "PROT")
#' 
#' \dontrun{
#' # Load ATAC-seq raw counts for hippocampus, excluding outliers 
#' data = load_sample_data("HIPPOC", 
#'                         "ATAC", 
#'                         exclude_outliers = TRUE, 
#'                         normalized = FALSE, 
#'                         scratchdir = "/tmp")
#' }
load_sample_data = function(tissue, 
                            assay, 
                            normalized = TRUE, 
                            training_regulated_only = FALSE, 
                            exclude_outliers = FALSE, 
                            scratchdir = ".", 
                            nrows = Inf,
                            warnings = TRUE){
  
  # check inputs 
  if(!tissue %in% MotrpacRatTraining6moData::TISSUE_ABBREV){
    stop(sprintf("'tissue' must be one of TISSUE_ABBREV: \n %s", paste(MotrpacRatTraining6moData::TISSUE_ABBREV, collapse=", ")))
  }
  if(!assay %in% MotrpacRatTraining6moData::ASSAY_ABBREV){
    stop(sprintf("'assay' must be one of ASSAY_ABBREV: \n %s", paste(MotrpacRatTraining6moData::ASSAY_ABBREV, collapse=", ")))
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
  # norm data 
  available_data = c(available_data, paste0(rep(c("ATAC", "METHYL"), each=8), "_", rep(possible_tissues,2), "_NORM_DATA"))
  # counts
  available_data = c(available_data, paste0(rep(c("ATAC", "METHYL"), each=8), "_", rep(possible_tissues,2), "_RAW_COUNTS"))
  
  data = NULL
  if(normalized){
    # normalized data
    if(assay %in% c("METHYL","ATAC")){
      if(training_regulated_only){
        obj_name = sprintf("%s_%s_NORM_DATA_05FDR", assay, gsub("-","",tissue))
        if(obj_name %in% available_data){
          message(obj_name)
          data = fetch_object(obj_name)
        }
      }else{
        # download from GCS
        obj_name = sprintf("%s_%s_NORM_DATA", assay, gsub("-","",tissue))
        if(obj_name %in% available_data){
          message(obj_name)
          data = get_rdata_from_url(tissue=tissue, assay=assay, suffix="NORM_DATA", scratchdir=scratchdir, nrows=nrows)
        }
      }
    }else if(assay == "METAB"){
      # get combined sample-level data
      message(sprintf("METAB %s normalized data from METAB_NORM_DATA_FLAT", tissue))
      data = MotrpacRatTraining6moData::METAB_NORM_DATA_FLAT
      # filter to this tissue
      data = data[data$tissue == tissue,]
      # remove all-NA columns 
      data[,names(colSums(is.na(data))[colSums(is.na(data)) == nrow(data)])] = NULL
    }else if(assay == "IMMUNO"){
      # get combined sample-level data
      message(sprintf("IMMUNO %s normalized data from IMMUNO_NORM_DATA_FLAT", tissue))
      data = MotrpacRatTraining6moData::IMMUNO_NORM_DATA_FLAT
      # filter to this tissue
      data = data[data$tissue == tissue,]
      # remove all-NA columns 
      data[,names(colSums(is.na(data))[colSums(is.na(data)) == nrow(data)])] = NULL
    }else{
      obj_name = sprintf("%s_%s_NORM_DATA", assay, gsub("-","",tissue))
      if(obj_name %in% available_data){
        message(obj_name)
        data = fetch_object(obj_name)
      }
    }
  }else{
    # raw data 
    if(assay %in% c("METHYL","ATAC","TRNSCRPT")){
      if(assay %in% c("METHYL","ATAC")){
        # download from GCS
        obj_name = sprintf("%s_%s_RAW_COUNTS", assay, gsub("-","",tissue))
        if(obj_name %in% available_data){
          message(obj_name)
          data = get_rdata_from_url(tissue=tissue, assay=assay, suffix="RAW_COUNTS", scratchdir=scratchdir, nrows=nrows)
        }
      }else{
        obj_name = sprintf("%s_%s_RAW_COUNTS", assay, gsub("-","",tissue))
        if(obj_name %in% available_data){
          message(obj_name)
          data = fetch_object(obj_name)
        }
      }
    }else{
      if(warnings){
        warning(sprintf("Non-normalized data not available for %s. Set 'normalized' to TRUE to get normalized sample-level data.", assay))
      }
      return()
    }
  }
  
  if(is.null(data)){
    if(warnings){
      warning(sprintf("No data returned for tissue %s and assay %s with current arguments.", tissue, assay))
    }
    return()
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
    outliers = data.table::data.table(MotrpacRatTraining6moData::OUTLIERS)
    .assay = assay
    .tissue = tissue
    curr_outliers = stats::na.omit(outliers[assay == .assay & tissue == .tissue, viallabel])
    if(.tissue == "VENACV"){
      curr_outliers = c(curr_outliers, na.omit(outliers[tissue == .tissue, viallabel]))
    }
    curr_outliers = unique(curr_outliers)
    if(length(curr_outliers)>0){
      data[,curr_outliers] = NULL
    }
  }
  
  return(data)
}


#' Load epigenetic differential analysis results 
#' 
#' Load ATAC or METHYL timewise differential analysis results for a tissue from 
#' Google Cloud Storage. See [MotrpacRatTraining6mo::ATAC_DA] and [MotrpacRatTraining6mo::METHYL_DA]
#' for more details. 
#'
#' @param tissue character, tissue abbreviation, one of "BAT", "HEART", "HIPPOC", "KIDNEY", "LIVER", "LUNG", "SKM-GN", "WAT-SC" 
#' @param assay character, one of "METHYL" (RRBS) or "ATAC"
#' @param scratchdir character, local directory in which to download data from 
#'   Google Cloud Storage. Current working directory by default.
#'
#' @return data frame. See [MotrpacRatTraining6mo::ATAC_DA] and [MotrpacRatTraining6mo::METHYL_DA] for details. 
#' @export
#'
#' @examples
#' \dontrun{
#' # Return ATAC-seq differential analysis results for gastrocnemius skeletal muscle
#' data = load_epigen_da("SKM-GN", "METHYL")
#' }
load_epigen_da = function(tissue, assay, scratchdir="."){
  data = get_rdata_from_url(tissue=tissue, assay=assay, suffix="DA", scratchdir=scratchdir)
  return(data)
}


#' Load raw METHYL data
#' 
#' Load METHYL raw data for a tissue from Google Cloud Storage. 
#' See [MotrpacRatTraining6moData::METHYL_RAW_DATA] for more details. 
#'
#' @param tissue character, tissue abbreviation, one of "BAT", "HEART", "HIPPOC", "KIDNEY", "LIVER", "LUNG", "SKM-GN", "WAT-SC" 
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
#' \dontrun{
#' # Load raw METHYL data for gastrocnemius 
#' data = load_methyl_raw_data("SKM-GN", "/tmp")
#' }
load_methyl_raw_data = function(tissue, scratchdir = "."){
  url = sprintf("https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_%s_RAW_DATA.rda", tissue)
  data = get_rdata_from_url(url = url, scratchdir = scratchdir)
  return(data)
}


#' Load METHYL feature annotation 
#' 
#' Load METHYL feature annotation from Google Cloud Storage. 
#' See [MotrpacRatTraining6mo::METHYL_FEATURE_ANNOT] for more details. 
#'
#' @param scratchdir character, local directory in which to download data from 
#'   Google Cloud Storage. Current working directory by default.
#'
#' @return data frame. See [MotrpacRatTraining6moData::METHYL_FEATURE_ANNOT] for details. 
#' @export
#'
#' @examples
#' \dontrun{
#' feature_annot = load_methyl_feature_annotation("/tmp")
#' }
#' 
#' @source <https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_FEATURE_ANNOT.rda> 
load_methyl_feature_annotation = function(scratchdir = "."){
  fa = get_rdata_from_url(url = "https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_FEATURE_ANNOT.rda",
                   scratchdir = scratchdir)
  return(fa)
}


#' Load ATAC feature annotation 
#' 
#' Load ATAC feature annotation from Google Cloud Storage. 
#' See [MotrpacRatTraining6mo::ATAC_FEATURE_ANNOT] for more details. 
#'
#' @param scratchdir character, local directory in which to download data from 
#'   Google Cloud Storage. Current working directory by default.
#'
#' @return data frame. See [MotrpacRatTraining6moData::ATAC_FEATURE_ANNOT] for details. 
#' @export
#'
#' @examples
#' \dontrun{
#' feature_annot = load_atac_feature_annotation("/tmp")
#' }
#' 
#' @source <https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_FEATURE_ANNOT.rda> 
load_atac_feature_annotation = function(scratchdir = "."){
  fa = get_rdata_from_url(url = "https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_FEATURE_ANNOT.rda",
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
#' @importFrom utils download.file 
#'
#' @details 
#' See the readme for this repository for all available files:
#' <https://github.com/MoTrPAC/MotrpacRatTraining6moData/blob/main/README.md>
#' 
#' @examples
#' \dontrun{
#' # return gastrocnemius chromatin accessibility differential analysis results 
#' data = get_rdata_from_url(tissue="SKM-GN", assay="ATAC", suffix="DA", scratchdir="/tmp")
#' 
#' # return brown adipose DNA methylation raw counts 
#' data = get_rdata_from_url(tissue="BAT", assay="METHYL", suffix="RAW_COUNTS", scratchdir="/tmp")
#' 
#' # return DNA methylation feature annotation
#' data = get_rdata_from_url(
#'   url="https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_FEATURE_ANNOT.rda", 
#'   scratchdir="/tmp"
#' )
#' 
#' # return raw DNA methylation data for brown adipose
#' data = get_rdata_from_url(
#'   url="https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/BAT_raw.RData", 
#'   scratchdir="/tmp"
#' )
#' }
get_rdata_from_url = function(tissue=NULL, assay=NULL, suffix=NULL, scratchdir=".", url=NULL, obj_name=NULL, nrows=Inf){
  # set option if default is low
  if(getOption("timeout") < 1e3){
    options(timeout=1e3)
  }

  # if dir doesn't exist, try to make it
  if(!dir.exists(scratchdir)){
    dir.create(scratchdir, recursive = TRUE)
  }
  
  if(!is.null(assay)){
    if(!assay %in% c("ATAC","METHYL")){
      warning("Only ATAC and METHYL data are currently available in GCS.")
    }
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
  if(!is.null(tissue)){
    if(!tissue %in% possible_tissues){
      warning(sprintf("Epigenetic data are only available for the following tissues:\n %s", paste(possible_tissues, collapse=", ")))
    }
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
  
  # delete local copy
  file.remove(local)

  return(data)
}


#' Filter outliers
#' 
#' Filter a list of outliers to those belonging to the specified dataset. 
#' Used to specify sex-specific outliers within differential analysis functions. 
#'
#' @param tissue optional `r tissue()` 
#' @param sex optional `r sex()` 
#' @param outliers vector of vial labels to consider as outliers.
#'   Defaults to vial labels in [MotrpacRatTraining6moData::OUTLIERS].
#'
#' @return character vector, subset of \code{outliers} that correspond to the 
#'   specified tissue and sex
#'   
#' @export
#'
#' @examples
#' curr_outliers = filter_outliers(tissue="HIPPOC")
#' curr_outliers = filter_outliers(tissue="HIPPOC", sex="male")
filter_outliers = function(tissue=NULL, sex=NULL, outliers=MotrpacRatTraining6moData::OUTLIERS$viallabel){
  TISSUE = tissue 
  SEX = sex
  outliers = as.character(outliers)
  if(length(outliers) == 0){
    return(character(0))
  }
  pheno = data.table::as.data.table(MotrpacRatTraining6moData::PHENO)
  if(!is.null(TISSUE)){
    pheno = pheno[tissue == TISSUE]
  }
  if(!is.null(SEX)){
    pheno = pheno[sex == SEX]
  }
  curr_outliers = pheno[viallabel %in% outliers, viallabel]
  return(curr_outliers)
}


#' Combine normalized sample-level data
#' 
#' Combine normalized sample-level data from the specified tissues and omes(s)/assay(s). 
#' If no tissues or omes are specified, all data is returned. 
#' In order to combine data from different tissues and data types, sample-specific vial labels 
#' are converted to animal-specific Participant IDs (PIDs).
#' 
#' @param tissues optional character vector of tissue abbreviations, 
#'   one of [MotrpacRatTraining6moData::TISSUE_ABBREV].
#' @param assays optional character vector of assay abbreviations, 
#'   one of [MotrpacRatTraining6moData::ASSAY_ABBREV]
#' @param include_epigen bool, whether to include the full ATAC or METHYL 
#'   differential analysis results from Google Cloud Storage. 
#'   Only relevant if \code{assays} includes "ATAC" or "METHYL".
#'   \code{FALSE} by default. 
#' @param scratchdir character, local directory in which to download data from the web. 
#'   Current working directory by default. Only relevant if \code{assays} includes "ATAC" or "METHYL".
#' @param training_regulated_only bool, whether to filter features down to those training-regulated at 5% FDR
#' @param exclude_outliers bool, whether to remove sample outliers specified by [MotrpacRatTraining6moData::OUTLIERS]
#' @param nrows integer, number of rows to return from each dataset. Defaults to Inf. 
#'   Useful to return a subset of a large data frame for tests. 
#' 
#' @export
#' 
#' @seealso [load_sample_data()], [viallabel_to_pid()]
#' 
#' @return data frame with features in rows and Participant IDs (PIDs) in columns
#' 
#' @examples
#' # Return all normalized RNA-seq data
#' data = combine_normalized_data(assays = "TRNSCRPT")
#' 
#' # Return all normalized proteomics data. Exclude outliers 
#' data = combine_normalized_data(assays = c("PROT","UBIQ","PHOSPHO","ACETYL"),
#'                                exclude_outliers = TRUE)
#' 
#' # Return normalized ATAC-seq data for training-regulated features 
#' data = combine_normalized_data(assays = "ATAC", training_regulated_only = TRUE)
#' 
#' # Return normalized ATAC-seq data for the first 1000 features in each tissue 
#' data = combine_normalized_data(assays = "ATAC", 
#'                                nrows = 1000, 
#'                                scratchdir = "/tmp", 
#'                                include_epigen = TRUE)
#'
#' # Return all normalized metabolomics data 
#' data = combine_normalized_data(assays = "METAB")
combine_normalized_data = function(tissues = MotrpacRatTraining6moData::TISSUE_ABBREV, 
                                   assays = MotrpacRatTraining6moData::ASSAY_ABBREV, 
                                   include_epigen = FALSE,
                                   scratchdir = ".", 
                                   training_regulated_only = FALSE, 
                                   exclude_outliers = FALSE, 
                                   nrows = Inf){
  
  if( ("ATAC" %in% assays | "METHYL" %in% assays) & !include_epigen & !training_regulated_only){
    warning("'include_epigen' is FALSE. Excluding ATAC and METHYL data.")
    assays = assays[!assays %in% c("ATAC","METHYL")]
  }
  
  # get current list of objects
  current_env = ls()
  
  res = list()
  for(a in assays){
    for(t in tissues){
      data = load_sample_data(t, a, 
                              normalized=TRUE, 
                              training_regulated_only=training_regulated_only, 
                              exclude_outliers=exclude_outliers, 
                              nrows=nrows, 
                              scratchdir=scratchdir,
                              warnings=FALSE)
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
      # add to result
      res[[sprintf("%s_%s",a,t)]] = data
    }
  }
  if(length(res)==0){
    warning(sprintf("No normalized data returned for tissues %s and assays %s.",
                    paste(tissues, collapse=", "),
                    paste(assays, collapse=", ")))
    return()
  }
  
  ## this doesn't work within the function 
  # # get current list of objects
  # new_env = ls()
  # 
  # # remove new ones to save space
  # remove = setdiff(new_env, current_env)
  # remove = remove[grepl(paste(assays,collapse="|"), remove)]
  # rm(list=remove)
  
  res = as.data.frame(data.table::rbindlist(res, fill=TRUE))
  return(res)
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
#' data = combine_da_results()
#' 
#' # Return all global proteomics differential analysis results
#' data = combine_da_results(assays="PROT")
#' 
#' \dontrun{
#' # Return METHYL and ATAC differential analysis results for gastrocnemius 
#' data = combine_da_results(tissues="SKM-GN", 
#'                           assays=c("ATAC","METHYL"),
#'                           include_epigen=TRUE)
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
