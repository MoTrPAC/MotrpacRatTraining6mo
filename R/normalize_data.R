#' Filter and normalize raw RNA-seq counts
#' 
#' Get raw counts and return the filtered and normalized data. 
#'
#' @param tissue `r tissue()`
#' @param min_cpm double, retain genes with more than \code{min_cpm} counts per million in at least \code{min_num_samples} samples
#' @param min_num_samples double, retain genes with more than \code{min_cpm} counts per million in at least \code{min_num_samples} samples
#' @param norm_method character, one of \code{c("TMM","TMMwsp","RLE","upperquartile","none")}. "TMM" by default.
#'
#' @return data frame where row names are feature_ID and column names are viallabel 
#' @export
#' 
#' @import MotrpacRatTraining6moData
#' @importFrom edgeR calcNormFactors cpm DGEList
#' 
#' @seealso [MotrpacRatTraining6moData::TRNSCRPT_NORM_DATA]
#'
#' @examples
#' norm_data = transcript_normalize_counts("LUNG")
#' 
transcript_normalize_counts = function(tissue, min_cpm = 0.5, min_num_samples = 2, norm_method="TMM"){

  counts = get(sprintf("TRNSCRPT_%s_RAW_COUNTS", gsub("-","",tissue)))
  rownames(counts) = counts$feature_ID
  counts[,c("feature","feature_ID","tissue","assay")] = NULL
  
  raw_dge = edgeR::DGEList(counts=counts) 
  keep = rowSums(cpm(raw_dge) > min_cpm) >= min_num_samples
  filt_dge = raw_dge[keep, , keep.lib.sizes=FALSE]
  
  # filt --> tmm 
  dge = edgeR::calcNormFactors(filt_dge, method=norm_method)
  tmm = edgeR::cpm(dge,log=TRUE)
  
  return(tmm)
}


#' Filter and normalize raw ATAC-seq counts
#' 
#' Download raw counts from Google Cloud Storage and return the filtered and 
#' quantile-normalized data. 
#'
#' @param tissue `r tissue()`
#' @param scratchdir character, local directory in which to download data from 
#'   Google Cloud Storage. Current working directory by default. 
#' @param n_samples integer, retain features with at least \code{min_count} counts in at least \code{n_samples} samples 
#' @param min_count integer, retain features with at least \code{min_count} counts in at least \code{n_samples} samples 
#'
#' @return data frame where row names are feature_ID and column names are viallabel 
#' @export
#' @importFrom limma voom
#'
#' @seealso [MotrpacRatTraining6moData::ATAC_NORM_DATA]
#'
#' @examples
#' norm_data = atac_normalize_counts("BAT", scratchdir = "/tmp")
atac_normalize_counts = function(tissue, scratchdir = ".", n_samples = 4, min_count = 10){
  
  # load raw counts 
  counts = load_sample_data(tissue = tissue, 
                            assay = "ATAC", 
                            normalized = FALSE, 
                            training_regulated_only = FALSE, 
                            scratchdir = scratchdir,
                            exclude_outliers = FALSE)
  rownames(counts) = counts$feature_ID
  counts[,c("feature","feature_ID","tissue","assay")] = NULL
  
  # remove non-auto peaks
  counts = counts[grepl("^chr[0-9]|^chrY|^chrX", rownames(counts)),]

  # exclude low count peaks in the current dataset
  # at least min_count counts in n_samples samples
  filt_counts = counts[rowSums(data.frame(lapply(counts, function(x) as.numeric(x >= min_count)), check.names=F)) >= n_samples,]

  # quantile normalize
  # this takes a couple of minutes given the size of the peak x sample counts matrix
  norm_counts = voom(filt_counts,normalize.method = "quantile")
  sub_norm = round(norm_counts$E,2)
  
  return(as.data.frame(sub_norm))
}
