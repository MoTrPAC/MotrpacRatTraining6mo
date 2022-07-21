#' Handle missing values in covariates for differential analysis
#' 
#' If fewer than 5% of values are missing from a continuous covariate, replace 
#' missing values with the mean in \code{meta}. Otherwise, remove the covariate from \code{covar}. 
#' If values are missing from a factor covariate, remove the covariate from \code{covar}. 
#' Note that covariates are only removed from \code{covar}, not \code{meta}.
#'
#' @param covar string or character vector of covariate names that correspond to column names of \code{meta}
#' @param meta sample by variable data frame of metadata
#'
#' @return named list of two items: 
#' \describe{
#'   \item{\code{meta}}{data frame input \code{meta} with covariates imputed as necessary}
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
#' 
fix_covariates = function(covar, meta){
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
  meta = as.data.frame(meta)
  return(list(meta=meta, covariates=covar))
}


#' Preprocess RNA-seq data 
#'
#' @param tissue character, tissue abbreviation, one of [TISSUE_ABBREV]
#' @param sex NULL to return data from both sexes or one of "male" or "female" to return data from the specified sex
#' @param covariates character vector of covariates that correspond to column names of [TRNSCRPT_META].
#'   Defaults to covariates that were used for the manuscript. 
#' @param outliers vector of viallabels to exclude during differential analysis. Defaults
#'   to \code{OUTLIERS$viallabel[OUTLIERS$assay == "TRNSCRPT"]}
#'
#' @return named list of five items: 
#' \describe{
#'   \item{\code{metadata}}{TODO}
#'   \item{\code{covariates}}TODO
#'   \item{\code{filt_counts}}{TODO}
#'   \item{\code{norm_data}}{TODO}
#'   \item{\code{outliers}}{TODO}
#' }
#' 
#' @seealso [OUTLIERS]
#' 
#' @export
#' @import data.table
#' @import MotrpacRatTraining6moData
#'
#' @examples
#' gastroc_data = transcript_prep_data("SKM-GN")
#' gastroc_data = transcript_prep_data("SKM-GN", outliers = NULL)
#' gastroc_data = transcript_prep_data("SKM-GN", covariates = NULL, outliers = NULL)
#' gastroc_data = transcript_prep_data("SKM-GN", covariates = NULL, outliers = NULL, sex = "male")
#' 
transcript_prep_data = function(tissue, 
                                sex = NULL, 
                                covariates = c('pct_globin', 'RIN', 'pct_umi_dup', 'median_5_3_bias'), 
                                outliers = OUTLIERS$viallabel[OUTLIERS$assay == "TRNSCRPT"]){
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
  
  # impute missing values
  new = fix_covariates(covariates, meta)
  covariates = new$covariates
  meta = data.table(new$meta)
  
  # center and scale continuous covariates
  for (cov in covariates){
    if(is.numeric(meta[,get(cov)])){
      meta[,(cov) := scale(meta[,get(cov)], center = T, scale = T)]
    }
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
