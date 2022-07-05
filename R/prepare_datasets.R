#' Retrieve phenotypic data 
#' 
#' Download and format phenotypic data from the endurance exercise training study in 6-month-old rats.
#' Codes are replaced with human-readable strings. 
#' Note that this currently does not handle variables that are recorded multiple times. 
#' 
#' @param scratch directory in which to download files  
#' @param gsutil_path path to \code{gsutil} binary. If you aren't sure where this is, run \code{which gsutil} in the Terminal and copy the result. 
#' @param parallel whether or not to re-download the file from GCP if it already exists locally 
#' 
#' @return data.frame with formatted phenotypic data; row names are sample identifiers (vial labels)
#' 
#' @export
#' @import data.table
#' 
#' @examples 
#' dl_format_pheno("~/Documents/motrpac_data", "/usr/local/bin/google-cloud-sdk/bin/gsutil")
#' dl_format_pheno("~/Documents/motrpac_data", "/usr/local/bin/google-cloud-sdk/bin/gsutil", parallel = FALSE)
#' 
dl_format_pheno = function(scratch, gsutil_path, parallel=TRUE){
  dmaqc_metadata = 'gs://motrpac-data-freeze-pass/pass1b-06/v1.0/results/phenotype/pass1b_6m_viallabel_data.txt'
  dmaqc_dict = 'gs://motrpac-data-freeze-pass/pass1b-06/v1.0/results/phenotype/merged_dictionary.txt'
  
  # download and format phenotypic data 
  dmaqc_metadata = dl_read_gcp(dmaqc_metadata, tmpdir = scratch, GSUTIL_PATH = gsutil_path, check_first = parallel)
  cols = dl_read_gcp(dmaqc_dict, tmpdir = scratch, GSUTIL_PATH = gsutil_path, check_first = parallel)
  old_cols = colnames(dmaqc_metadata)
  new_cols = tolower(cols[match(old_cols, BICUniqueID), FullName]) 
  colnames(dmaqc_metadata) = new_cols # this isn't perfect, but we don't care about the columns it doesn't work for for now 
  
  # make some variables human-readable
  # create new variables "protocol", "agegroup", "intervention", "sacrificetime", "sex" with readable strings 
  for (var in c('key.protocol','key.agegroup','key.intervention','key.sacrificetime','registration.sex')){
    d = cols[Field.Name == gsub('.*\\.','',var)]
    keys=unname(unlist(strsplit(d[,Categorical.Values],'\\|')))
    values=tolower(unname(unlist(strsplit(d[,Categorical.Definitions],'\\|'))))
    names(values) = keys
    # match keys to values; create new column 
    new_var = gsub(".*\\.","",var)
    dmaqc_metadata[,(new_var) := unname(values)[match(get(var), names(values))]]
  }
  dmaqc_metadata[,time_to_freeze := calculated.variables.frozetime_after_train - calculated.variables.deathtime_after_train]
  
  # clean up "sacrificetime"
  dmaqc_metadata[,sacrificetime := sapply(sacrificetime, function(x) gsub(' week.*','w',x))]
  
  # clean up 'intervention'
  dmaqc_metadata[grepl('training',intervention), intervention := 'training']
  
  # make "group" - "1w", "2w", "4w", "8w", "control"
  dmaqc_metadata[,group := sacrificetime]
  dmaqc_metadata[intervention == 'control', group := 'control']
  
  # make tech ID a string
  dmaqc_metadata[,specimen.processing.techid := paste0('tech',specimen.processing.techid)]
  
  # make viallabel char
  dmaqc_metadata[,viallabel := as.character(viallabel)]
  
  # convert to data.frame
  dmaqc_metadata_df = as.data.frame(dmaqc_metadata)
  rownames(dmaqc_metadata_df) = dmaqc_metadata_df$viallabel
  
  return(dmaqc_metadata)
}


#' Title
#'
#' @param TISSUE_CODE 
#' @param SEX 
#' @param gsutil_path 
#' @param scratch 
#' @param outliers 
#' @param parallel 
#'
#' @return
#' @export
#'
#' @examples
#' TODO
#' 
preprocess_pass1b_rnaseq_gcp = function(TISSUE_CODE, SEX, gsutil_path='~/google-cloud-sdk/bin/gsutil', scratch='/oak/stanford/groups/smontgom/nicolerg/tmp', outliers=NULL, parallel=F){
  
  if(!SEX %in% c('male','female','all')){
    stop('"SEX" must take on one of the following values:\n  "male","female","all"\n')
  }
  
  # fix some things
  if(TISSUE_CODE == 't54-hypothalmus'){
    TISSUE_CODE = 't54-hypothalamus'
  }
  
  dmaqc_metadata = dl_format_pheno(scratch, gsutil_path, parallel=parallel)
  
  # load counts 
  counts = dl_read_gcp(sprintf('gs://motrpac-data-freeze-pass/pass1b-06/v1.0/results/transcriptomics/%s/transcript-rna-seq/motrpac_pass1b-06_%s_transcript-rna-seq_rsem-genes-count.txt', TISSUE_CODE, TISSUE_CODE), 
                       tmpdir = scratch, 
                       GSUTIL_PATH = gsutil_path,
                       check_first = parallel)
  if(is.null(counts)){return()}
  counts = as.data.frame(counts)
  rownames(counts) = counts$gene_id
  counts$gene_id = NULL
  raw_counts = counts
  
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
