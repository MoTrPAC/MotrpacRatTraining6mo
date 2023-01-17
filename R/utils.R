#' Map vial labels to PIDs
#' 
#' Map sample identifiers (vial labels) to participant IDs (PIDs).
#' Vial labels are unique per sample; PIDs are unique per animal. 
#' Multiple vial labels correspond to different aliquots of the same tissue sample. 
#'
#' @param viallabels vector, list of vial labels  
#'
#' @return named list where names are vial labels and values are PIDs 
#' @export
#'
#' @examples
#' map = viallabel_to_pid(c("90416015402", "90416015403", "90416015302"))
#' map = viallabel_to_pid(c(90416015402, 90416015403, 90416015302))
#' map 
viallabel_to_pid = function(viallabels){
  pheno = data.table::as.data.table(MotrpacRatTraining6moData::PHENO)
  pheno = unique(pheno[,.(viallabel, pid)])
  pheno = pheno[viallabel %in% as.character(viallabels)]
  vl_to_pid = pheno[,pid]
  names(vl_to_pid) = pheno[,viallabel]
  # if possible, make order same as input
  if(length(viallabels) == length(vl_to_pid)){
    vl_to_pid = vl_to_pid[viallabels]
  }
  return(vl_to_pid)
}


#' List available data
#' 
#' List available data, including lazily-loaded data. Useful to gather
#' split data frames.
#'
#' @param package optional string to specify a package 
#'
#' @return character vector of names of data objects available to load
#' @export
#'
#' @examples
#' list_available_data()
#' list_available_data("MotrpacRatTraining6moData")
#' 
list_available_data = function(package=NULL){
  res = utils::data(package=package)
  obj = res$results[,3]
  # remove objects that can't be called directly
  obj = obj[!grepl("\\(", obj)]
  return(obj)
}


#' Check arguments for DA functions 
#' 
#' Internal function used to check arguments for differential analysis functions
#'
#' @param tissue `r tissue()`
#' @param outfile character, output file 
#' @param overwrite bool, whether to overwrite \code{outfile} if it exists
#' @param outfile_is_rdata bool, whether \code{outfile} is intended to save RData
#'
#' @keywords internal
#' 
check_da_args = function(tissue, outfile, overwrite, outfile_is_rdata = TRUE){
  # check arguments 
  if(length(tissue)>1){
    stop("Please specify a single tissue, e.g., 'BAT' for brown adipose tissue. See 'TISSUE_ABBREV' for options.")
  }
  if(!is.null(outfile)){
    if(!overwrite & file.exists(outfile)){
      stop(sprintf("'%s' already exists and 'overwrite' = FALSE. Specify a new outfile or set 'overwrite' to TRUE to generate new results.",
                   outfile))
    }
    if(outfile_is_rdata){
      if(!any(unlist(lapply(c("\\.rda$", "\\.rdata"), function(pattern){
        grepl(pattern, outfile, ignore.case = TRUE)
      })))){
        stop(sprintf("Outfile '%s' should end with '.rda' or '.RData' to indicate an RData file (not case-sensitive).", outfile))
      }
    }
    # check if path is valid
    dir = dirname(outfile)
    if(dir != "" & !dir.exists(dir)){
      message(sprintf("Creating directory for outfile: %s", dir))
      dir.create(dir, recursive = TRUE) # this will throw an error if it's not a valid path/uncreatable 
    }
  }
  return()
}


#' Get genomic peak annotations
#' 
#' Get and fix peak annotations from [ChIPseeker::annotatePeak()]
#' 
#' @param counts data table or data frame with columns 'chrom','start','end', and 'feature_ID' OR 
#'   for ATAC only, data frame with row names as feature IDs in the format format 'chrom:start-end'
#' @param species character, scientific name of the organism for which to create the 
#'   Ensembl database with [GenomicFeatures::makeTxDbFromEnsembl()]
#' @param release integer, the Ensembl release to query. 
#'   If set to NA, the current release is used.  
#' @param txdb optional [GenomicFeatures::TxDb-class] object, in which case 
#'   the database is not regenerated 
#' 
#' @export
#' 
#' @return data frame with formatted [ChIPseeker::annotatePeak()] output 
#'   and additional columns "custom_annotation" and "relationship_to_gene". 
#'   
#' @details 
#'   "relationship_to_gene" is the shortest distance between the feature and the start or end of the closest gene. 
#'   It is 0 if the feature has any overlap with the gene. 
#'   "custom_annotation" fixes many issues with the \code{ChIPseeker} annotation (v1.22.1). 
#'   
#' @examples 
#' \dontrun{
#' # Load raw ATAC-seq counts for one tissue 
#' raw_counts = load_sample_data("SKM-GN",
#'                               "ATAC",
#'                               normalized = FALSE,
#'                               scratchdir = "/tmp")
#' pa = get_peak_annotations(raw_counts)
#' }
get_peak_annotations = function(counts, species="Rattus norvegicus", release=96, txdb=NULL){

  for (pkg in c('GenomeInfoDb','GenomicRanges','ChIPseeker','IRanges')){
    if (!requireNamespace(pkg, quietly = TRUE)){
      stop(
        sprintf("Package '%s' must be installed to generate a TxDb object.", pkg),
        call. = FALSE
      )
    }
  }
  
  if(!"feature_ID" %in% colnames(counts) & !data.table::is.data.table(counts)){
    genomic_peaks = data.table::data.table(
      feature_ID = rownames(counts),
      chrom = gsub(":.*","",rownames(counts)),
      start = as.numeric(gsub(".*:|-.*","",rownames(counts))),
      end = as.numeric(gsub(".*-","",rownames(counts)))
    )
  }else if(!"feature_ID" %in% colnames(counts) & data.table::is.data.table(counts)){
    counts = counts
    counts[,feature_ID := paste0(chrom,':',start,'-',end)]
    genomic_peaks = counts[,.(chrom, start, end, feature_ID)]
  }else if("feature_ID" %in% colnames(counts) & data.table::is.data.table(counts)){
    counts = counts
    genomic_peaks = counts[,.(chrom, start, end, feature_ID)]
  }else if("feature_ID" %in% colnames(counts) & data.table::is.data.table(counts)){
    counts = data.table::as.data.table(counts)
    genomic_peaks = counts[,.(chrom, start, end, feature_ID)]
  }else{
    stop("Incorrect input format. 'counts' should be a data frame or data table with columns 'chrom', 'start', 'end', and 'feature_ID'.")
  }
  
  if(is.null(txdb)){

    if (!requireNamespace("GenomicFeatures", quietly = TRUE) | !requireNamespace("RMariaDB", quietly = TRUE)){
      stop(
        "Packages 'GenomicFeatures' and 'RMariaDB' must be installed to generate a TxDb object with 'makeTxDbFromEnsembl()'.",
        call. = FALSE
      )
    }

    # make TxDb object 
    txdb = GenomicFeatures::makeTxDbFromEnsembl(organism=species,
                                                release=release)
  }
  
  accepted_chrom = GenomeInfoDb::seqlevels(txdb)
  # remove contigs
  accepted_chrom = accepted_chrom[!grepl("\\.",accepted_chrom)]
  
  # remove contigs from input 
  genomic_peaks = genomic_peaks[!grepl("\\.",chrom)]
  
  genomic_peaks[,chrom := gsub("^chr","",as.character(chrom))]
  if(!all(unique(genomic_peaks[,chrom]) %in% accepted_chrom)){
    stop(sprintf("The following chromosomes are found in the input but not in the txdb object: %s", paste0(unique(!genomic_peaks[,chrom] %in% accepted_chrom), collapse=', ')))
  }
  
  # annotate peaks 
  peak = GenomicRanges::GRanges(seqnames = genomic_peaks[,chrom], 
                                ranges = IRanges::IRanges(as.numeric(genomic_peaks[,start]), as.numeric(genomic_peaks[,end])))
  peakAnno = ChIPseeker::annotatePeak(peak, 
                                      level = "gene",
                                      tssRegion=c(-2000,1000), 
                                      TxDb=txdb,
                                      overlap = "all")
  pa = data.table::as.data.table(peakAnno@anno)
  
  # add back feature_ID
  if(nrow(pa)==nrow(genomic_peaks)){
    pa[,feature_ID := genomic_peaks[,feature_ID]]
  }else{
    cols=c('seqnames','start','end')
    pa[,(cols) := lapply(.SD, as.character), .SDcols=cols]
    cols=c('chrom','start','end')
    genomic_peaks[,(cols) := lapply(.SD, as.character), .SDcols=cols]
    pa = merge(pa, genomic_peaks, by.x=c('seqnames','start','end'), by.y=c('chrom','start','end'), all.y=T)
  }
  
  # add simpler annotation
  pa[,short_annotation := annotation]
  pa[grepl('Exon', short_annotation), short_annotation := 'Exon']
  pa[grepl('Intron', short_annotation), short_annotation := 'Intron']
  
  pa[,c('geneChr','strand') := NULL]
  
  # add a "relationship_to_gene" (from gene start or end, whichever is closer)
  cols=c('start','end','geneStart','geneEnd','geneStrand')
  pa[,(cols) := lapply(.SD, as.numeric), .SDcols=cols]
  pa[,dist_upstream := ifelse(end-geneStart <=0, end-geneStart, NA_real_)]
  pa[,dist_downstream := ifelse(start-geneEnd >= 0, start-geneEnd, NA_real_)]
  # if feature overlaps with gene body, say dist_to_gene is 0
  pa[end >= geneStart & start <= geneEnd, dist_downstream := 0]
  pa[end >= geneStart & start <= geneEnd, dist_upstream := 0]
  pa[, relationship_to_gene := ifelse(is.na(dist_downstream), dist_upstream, dist_downstream)]
  pa[,c('dist_upstream','dist_downstream') := NULL]
  
  # fix "Downstream" and "Distal Intergenic" annotations
  pa[relationship_to_gene == 0 & grepl("Downstream|Intergenic", short_annotation), short_annotation := "Overlaps Gene"]
  # Downstream
  pa[geneStrand == 1 & relationship_to_gene > 0 & relationship_to_gene < 5000, short_annotation := "Downstream (<5kb)"]
  pa[geneStrand == 2 & relationship_to_gene < 0 & relationship_to_gene > -5000, short_annotation := "Downstream (<5kb)"]
  # Upstream (promoter excluded)
  pa[geneStrand == 1 & relationship_to_gene > -5000 & relationship_to_gene < 0 & 
       grepl("Downstream|Intergenic", short_annotation), short_annotation := "Upstream (<5kb)"]
  pa[geneStrand == 2 & relationship_to_gene < 5000 & relationship_to_gene > 0 & 
       grepl("Downstream|Intergenic", short_annotation), short_annotation := "Upstream (<5kb)"]
  pa[abs(relationship_to_gene) >= 5000, short_annotation := "Distal Intergenic"]

  data.table::setnames(pa, 
                       c('short_annotation','annotation','seqnames','geneId'), 
                       c('custom_annotation','chipseeker_annotation','chrom','ensembl_gene'))
  
  return(as.data.frame(pa))
}


#' Extract standard errors
#' 
#' Extract standard errors from \code{limma} results
#'
#' @param limma_res data frame returned by [limma::topTable()]
#' @param effect_col character, name of column containing the effect size
#' @param t_col character, column name containing the t statistic
#'
#' @return numeric vector of standard errors for effect sizes
#' 
#' @keywords internal
#'
limma_res_extract_se = function(limma_res,
                                effect_col="logFC",
                                t_col="t"){
  effects = limma_res[[effect_col]]
  ts = limma_res[[t_col]]
  ses1 = effects/ts
  return(ses1)
}


#' Make data frame numeric only
#' 
#' Set row names and remove all non-numeric columns. This is useful for
#' reformatting data objects in `MotrpacRatTraining6moData`,
#' e.g., [MotrpacRatTraining6moData::PROT_HEART_NORM_DATA]
#'
#' @param df data frame with at least one numeric column
#' @param rownames either NULL, a character specifying a column name in \code{df},
#'   or a vector. Specifies target row names. Default: "feature_ID" 
#'
#' @return numeric data frame
#' @export
#'
#' @examples
#' df = MotrpacRatTraining6moData::PROT_HEART_NORM_DATA
#' head(df)
#' head(df_to_numeric(df))
#' 
#' df = load_sample_data("SKM-GN", "TRNSCRPT")
#' head(df)
#' head(df_to_numeric(df))
#' 
#' df = MotrpacRatTraining6moData::METAB_NORM_DATA_FLAT 
#' head(df)
#' rn = paste(df$assay, df$tissue, df$feature_ID, df$dataset, sep=";")
#' head(df_to_numeric(df, rownames = rn))
df_to_numeric = function(df, rownames="feature_ID"){
  
  df = as.data.frame(df)
  
  # assign rownames
  if(is.null(rownames)){
    rownames(df) = NULL
  }else{
    if(length(rownames) == 1){
      if(!rownames %in% colnames(df)){
        stop(sprintf("'%s' not found in column names of input data frame. 'rownames' must specify a valid column.", rownames))
      }
      if(any(duplicated(df[,rownames]))){
        stop(sprintf("Row names must be unique, but df$%s has duplicate values.", rownames))
      }
      rownames(df) = df[,rownames]
    }else{
      if(any(duplicated(rownames))){
        stop("Row names must be unique, but 'rownames' has duplicate values.")
      }
      if(length(rownames) != nrow(df)){
        stop("Length of 'rownames' does not equal number of rows in 'df'.")
      }
      rownames(df) = as.character(rownames)
    }
  }
  
  # remove non-numeric colnames 
  non_num_cols = colnames(df)[unlist(lapply(df, function(x) !is.numeric(x)))]
  if(length(non_num_cols) == ncol(df)){
    stop("Input data frame has no non-numeric columns.")
  }
  df[,non_num_cols] = NULL

  return(df)  
}


#' Get object from MotrpacRatTraining6moData
#' 
#' For internal use only. Users can use \code{get()} or \code{data()} directly
#' after attaching the \code{MotrpacRatTraining6moData} package. 
#'
#' @param object_name_as_string character, name of data object in \code{MotrpacRatTraining6moData}
#'   R package 
#'
#' @return specified object
#'
#' @keywords internal
.get = function(object_name_as_string){
  data = get(object_name_as_string, envir=as.environment("package:MotrpacRatTraining6moData"))
  return(data)
}
