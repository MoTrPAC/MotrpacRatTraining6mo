#' Filter and normalize raw RNA-seq counts
#' 
#' Load raw counts and return the filtered and normalized data. 
#' Alternatively, the user can provide a numeric data frame of raw RNA-seq counts. 
#'
#' @param tissue `r tissue()`
#' @param min_cpm double, retain genes with more than \code{min_cpm} counts per million in at least \code{min_num_samples} samples
#' @param min_num_samples double, retain genes with more than \code{min_cpm} counts per million in at least \code{min_num_samples} samples
#' @param norm_method character, one of \code{c("TMM","TMMwsp","RLE","upperquartile","none")}. "TMM" by default.
#' @param counts optional user-supplied numeric data frame or matrix where row names are 
#'   gene IDs and column names are sample identifiers
#'
#' @return data frame where row names are feature_ID and column names are viallabel 
#' @export
#' 
#' @importFrom edgeR calcNormFactors cpm DGEList
#' 
#' @seealso [MotrpacRatTraining6moData::TRNSCRPT_NORM_DATA]
#'
#' @examples
#' norm_data = transcript_normalize_counts("LUNG")
#' 
#' # Simulate "user-supplied data"
#' counts = load_sample_data("LUNG", "TRNSCRPT", normalized=FALSE)
#' counts = df_to_numeric(counts)
#' norm_data = transcript_normalize_counts(counts = counts)
#' 
#' @details 
#' Note that while this function is identical to the code used to generate the 
#' normalized RNA-seq data tables ([MotrpacRatTraining6moData::TRNSCRPT_NORM_DATA])
#' and the normalized RNA-seq data available through the MoTrPAC Data Hub,
#' \code{transcript_normalize_counts(tissue)} yields slightly fewer genes than its
#' corresponding [MotrpacRatTraining6moData::TRNSCRPT_NORM_DATA] object. 
#' Investigation of this discrepancy suggests minor functional differences in the 
#' version of [edgeR::cpm()] used ~2.5 years apart. Find more details in 
#' [this GitHub issue](https://github.com/MoTrPAC/MotrpacRatTraining6moData/issues/41).  
#' 
transcript_normalize_counts = function(tissue, min_cpm = 0.5, min_num_samples = 2, norm_method = "TMM", counts = NULL){

  if(!is.null(counts)){
    counts = as.data.frame(counts)
    if(!all(apply(counts, c(1,2), is.numeric))){
      stop("'counts' must be a numeric table with features as row names and samples as column names.")
    }
    if(all(rownames(counts) == as.character(1:nrow(counts)))){
      stop("'counts' should have meaningful row names, i.e., feature IDs.")
    }
  }else{
    obj_name = sprintf("TRNSCRPT_%s_RAW_COUNTS", gsub("-","",tissue))
    counts = .get(obj_name)
    counts = df_to_numeric(counts)
  }
  
  raw_dge = edgeR::DGEList(counts=counts) 
  keep = rowSums(edgeR::cpm(raw_dge) > min_cpm) >= min_num_samples
  filt_dge = raw_dge[keep, , keep.lib.sizes=FALSE]
  
  # filt --> tmm 
  dge = edgeR::calcNormFactors(filt_dge, method=norm_method)
  tmm = edgeR::cpm(dge,log=TRUE)
  
  return(tmm)
}


#' Filter and normalize raw ATAC-seq counts
#' 
#' Download raw counts from Google Cloud Storage and return the filtered and 
#' quantile-normalized data. Alternatively, the user can provide a numeric
#' data frame of raw ATAC-seq counts. 
#'
#' @param tissue `r tissue()`
#' @param scratchdir character, local directory in which to download data from 
#'   Google Cloud Storage. Current working directory by default. 
#' @param n_samples integer, retain features with at least \code{min_count} counts in at least \code{n_samples} samples 
#' @param min_count integer, retain features with at least \code{min_count} counts in at least \code{n_samples} samples 
#' @param counts optional user-supplied numeric data frame or matrix where row names are 
#'   feature IDs and column names are sample identifiers 
#'
#' @return data frame where row names are feature_ID and column names are viallabel 
#' @export
#' @importFrom limma voom
#'
#' @seealso [MotrpacRatTraining6moData::ATAC_NORM_DATA]
#'
#' @examples
#' \dontrun{
#' norm_data = atac_normalize_counts("BAT", scratchdir = "/tmp")
#' }
#' 
#' @details 
#' Non-autosomal peaks are removed, i.e., peak IDs that don't begin with "chrX", "chrY",
#' or "chr"+number. If you are providing your own counts matrix, ensure that peak IDs
#' follow the same naming convention. 
#' 
atac_normalize_counts = function(tissue, scratchdir = ".", n_samples = 4, min_count = 10, counts = NULL){
  
  if(!is.null(counts)){
    counts = as.data.frame(counts)
    if(!all(apply(counts, c(1,2), is.numeric))){
      stop("'counts' must be a numeric table with features as row names and samples as column names.")
    }
    if(all(rownames(counts) == as.character(1:nrow(counts)))){
      stop("'counts' should have meaningful row names, i.e., feature IDs.")
    }
    if(!any(grepl("^chr[0-9]|^chrY|^chrX", rownames(counts)))){
      stop(paste("No autosomal features found in input 'counts'. Did you use the expected chromosome naming convention?",
                 "i.e., chr1, chr2, ..., chr20, chrX, chrY."))
    }
  }else{
    # load raw counts 
    counts = load_sample_data(tissue = tissue, 
                              assay = "ATAC", 
                              normalized = FALSE, 
                              training_regulated_only = FALSE, 
                              scratchdir = scratchdir,
                              exclude_outliers = FALSE)
    counts = df_to_numeric(counts)
  }
  
  # remove non-auto peaks
  counts = counts[grepl("^chr[0-9]|^chrY|^chrX", rownames(counts)),]

  # exclude low count peaks in the current dataset
  # at least min_count counts in n_samples samples
  filt_counts = counts[rowSums(data.frame(lapply(counts, function(x) as.numeric(x >= min_count)), check.names=FALSE)) >= n_samples,]

  # quantile normalize
  # this takes a couple of minutes given the size of the peak x sample counts matrix
  norm_counts = limma::voom(filt_counts,normalize.method = "quantile")
  sub_norm = round(norm_counts$E,2)
  
  return(as.data.frame(sub_norm))
}
