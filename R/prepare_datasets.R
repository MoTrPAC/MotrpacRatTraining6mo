#' Title
#' 
#' Description
#'
#' @param tissue character, tissue abbreviation, one of [TISSUE_ABBREV], excluding "PLASMA"
#' @param sex character, "male", "female", or "all"
#' @param outliers optional, vector of vial labels to exclude 
#'
#' @return
#' @export
#' 
#' @import data.table
#' @import MotrpacRatTraining6moData
#'
#' @examples
#' TODO
#' 
prep_data_TRNSCRPT = function(tissue, sex, outliers = NULL){
  
  if(!sex %in% c('male','female','all')){
    stop('"sex" must be one of the following values:\n  "male","female","all"\n')
  }
  
  if(!tissue %in% MotrpacRatTraining6moData::TISSUE_ABBREV){
    stop(sprintf('"tissue" must be one of the following values:\n  %s\n',
                 paste0(MotrpacRatTraining6moData::TISSUE_ABBREV, collapse=", ")))
  }
  
  dmaqc_metadata = as.data.table(MotrpacRatTraining6moData::PHENO)
  
  # load counts 
  # will return error if data doesn't exist
  counts = get(sprintf("TRNSCRPT_%s_RAW_COUNTS", gsub("-","",tissue)))
  raw_counts = counts
  
  # get normalized data
  tmm = MotrpacRatTraining6moData::TRNSCRPT_SAMPLE_DATA
  # select this tissue
  tmm = tmm[grepl(sprintf("TRNSCRPT;%s", tissue), tmm$feature),]
  # convert PID to viallabel 
  
  
  # filter genes
  # TODO
  
  # get nor
  
  # filter by genes in normalized data 
  tmm = dl_read_gcp(sprintf('gs://motrpac-data-freeze-pass/pass1b-06/v1.0/analysis/transcriptomics/transcript-rna-seq/normalized-data/motrpac_pass1b-06_%s_transcript-rna-seq_normalized-log-cpm.txt', TISSUE_CODE), 
                    tmpdir = scratch, 
                    GSUTIL_PATH = gsutil_path,
                    check_first = parallel)
  # remove pid and bid rows
  tmm = tmm[3:nrow(tmm)]
  counts = counts[tmm[,viallabel],]
  
  # convert tmm to data.frame
  tmm = as.data.frame(tmm)
  rownames(tmm) = tmm$viallabel
  tmm$viallabel = NULL
  
  # subset metadata 
  meta_data = dmaqc_metadata[as.character(viallabel) %in% colnames(tmm)]
  # remove reference standards from counts 
  counts = counts[,as.character(meta_data[,viallabel])]
  
  # coerce counts to integer (RSEM has fractional values for counts sometimes)
  counts_round = as.data.frame(apply(counts, c(1,2), as.integer)) 
  raw_counts_round = as.data.frame(apply(raw_counts, c(1,2), as.integer)) 
  
  # RNA-seq metadata 
  qa_qc = dl_read_gcp('gs://motrpac-data-freeze-pass/pass1b-06/v1.0/results/transcriptomics/qa-qc/motrpac_pass1b-06_transcript-rna-seq_qa-qc-metrics.csv', 
                      sep=',', 
                      tmpdir = scratch, 
                      GSUTIL_PATH = gsutil_path,
                      check_first = parallel)
  # adjust column names
  colnames(qa_qc) = tolower(gsub(' .*','',colnames(qa_qc)))
  # remove duplicate columns
  qa_qc[,c('pid','bid') := NULL]
  meta_data[,viallabel := as.character(viallabel)]
  qa_qc[,vial_label := as.character(vial_label)]
  meta = merge(meta_data, qa_qc, by.x = 'viallabel', by.y = 'vial_label')
  
  if(SEX != 'all'){
    # subset to sex 
    meta = meta[sex == SEX]
    tmm = tmm[as.character(meta[,viallabel])]
    counts_round = counts_round[as.character(meta[,viallabel])]
  }
  
  # remove outliers 
  if(!is.null(outliers)){
    meta = meta[,viallabel := as.character(viallabel)]
    meta = meta[!viallabel %in% outliers]
    tmm = tmm[meta[,viallabel]]
    counts_round = counts_round[meta[,viallabel]]
  }
  
  return(list(meta = meta,
              raw_counts = raw_counts_round, 
              counts = counts_round,
              norm = tmm))
}


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
