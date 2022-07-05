#' Title
#'
#' @param tissue_code 
#' @param meta 
#' @param counts 
#' @param tmm 
#' @param covariates 
#' @param outliers 
#' @param gsutil_path 
#' @param parallel 
#'
#' @return TODO
#' 
#' @export
#' @import data.table
#'
#' @examples
#' TODO
#' 
transcript_prep_data = function(tissue_code, 
                                meta, 
                                counts, 
                                tmm,
                                covariates=c('pct_globin', 'rin', 'pct_umi_dup', 'median_5_3_bias'), 
                                outliers=NULL, 
                                gsutil_path='~/google-cloud-sdk/bin/gsutil',
                                parallel=F){
  
  if(tissue_code %in% c("t31-plasma", "t57-tibia")){return(NULL)}
  
  # fix some inconsistencies 
  if(tissue_code == 't54-hypothalmus'){
    tissue_code = 't54-hypothalamus'
  }
  
  # load outliers
  if(is.null(outliers)){
    rna_outliers = dl_read_gcp('gs://mawg-data/pass1b-06/transcript-rna-seq/dea/pass1b-06_transcript-rna-seq_removed-outliers_20201028.txt', sep='\t', check_first=parallel)
    outliers = as.character(rna_outliers[,viallabel])
  }
  
  cat(tissue_code, sep = '\n')
  
  # remove outliers 
  curr_outliers = outliers[outliers %in% as.character(meta[,viallabel])]
  if(length(curr_outliers)>0){
    meta = meta[,viallabel := as.character(viallabel)]
    meta = meta[!viallabel %in% outliers]
    counts = counts[meta[,viallabel]]
    tmm = tmm[meta[,viallabel]]
  }
  if(tissue_code == 't65-aorta'){
    # add Ucp1 as a covariate
    ucp1 = data.table(viallabel = colnames(counts), ucp1 = unname(unlist(counts['ENSRNOG00000003580',])))
    meta = merge(meta, ucp1, by = 'viallabel')
    covariates = c(covariates, 'ucp1')
  }
  
  # impute missing values
  new = fix_missing(covariates, meta)
  covariates = new$covariates
  meta = new$meta
  
  # center and scale continuous covariates
  for (cov in covariates){
    if(is.numeric(meta[,get(cov)])){
      meta[,(cov) := scale(meta[,get(cov)], center = T, scale = T)]
    }
  }
  
  meta[,sex_group := paste0(sex, ';', group)]
  
  return(list(fixed_meta = meta, 
              fixed_covariates = covariates, 
              fixed_counts = counts,
              fixed_tmm = tmm, 
              curr_outliers = curr_outliers))
}


#' Title
#'
#' @param tissue_code 
#' @param meta 
#' @param counts 
#' @param covariates 
#' @param curr_outliers 
#' @param date 
#' @param save_rdata 
#' @param write 
#'
#' @return TODO
#' 
#' @export 
#' @import data.table
#'
#' @examples
#' TODO
#' 
transcript_timewise_dea_each_sex = function(tissue_code, meta, counts, covariates, curr_outliers, date, save_rdata=T, write=T){
  
  # fix some inconsistencies 
  if(tissue_code == 't54-hypothalmus'){
    tissue_code = 't54-hypothalamus'
  }
  
  outfile = sprintf('dea/pass1b-06_%s_transcript-rna-seq_timewise-dea_%s.txt',
                    tissue_code,date)
  if(file.exists(outfile)){
    dt = fread(outfile, sep='\t', header=T)
    return(dt)
  }
  
  sex_res = list()
  for(SEX in unique(meta[,sex])){
    
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
    
    # shrunk results 
    # function in pi1_cook_fx.R
    deseq_res_shrunk = run_deseq(curr_counts, # filtered counts
                                 curr_meta, # metadata
                                 covariates, # covariates
                                 'group', # outcome of interest
                                 contrasts, # list of contrasts in format c(outcome_of_interest, numerator_level, denominator_level)
                                 shrink = T)
    
    # non-shrunk results 
    deseq_res = run_deseq(curr_counts, # filtered counts
                          curr_meta, # metadata
                          covariates, # covariates
                          'group', # outcome of interest
                          contrasts, # list of contrasts in format c(outcome_of_interest, numerator_level, denominator_level)
                          shrink = F)
    
    if(save_rdata){
      save(deseq_res, deseq_res_shrunk, file=sprintf('rdata/%s_%s_timewise-dea_%s.RData', tissue_code, SEX, date))
    }
    
    # collect res
    res_shrunk = data.table(deseq_res_shrunk$res)
    res_nonshrunk = data.table(deseq_res$res)
    setnames(res_shrunk, c("log2FoldChange", "lfcSE"), c("shrunk_logFC","shrunk_logFC_se"))
    setnames(res_nonshrunk, c("log2FoldChange", "lfcSE", "stat"), c("logFC","logFC_se", "zscore"))
    res_shrunk = res_shrunk[,.(gene_id, shrunk_logFC, shrunk_logFC_se, numerator, denominator)]
    res = merge(res_nonshrunk, res_shrunk, by=c("gene_id","numerator","denominator"))
    res[,sex := SEX]
    
    setnames(res, c("numerator","pvalue","gene_id"), c("comparison_group","p_value","feature_ID"))
    res[,denominator := NULL]
    
    # add some columns
    res[,tissue := tissue_code]
    res[,assay := 'transcript-rna-seq']
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
  
  dt = dt[,.(
    feature_ID,
    sex,
    comparison_group,
    assay,
    tissue,
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
  
  # if aorta, remove 1w, 2w F
  if(tissue_code=="t65-aorta"){
    dt = dt[!(sex=='female' & comparison_group %in% c('1w','2w'))]
  }
  
  if(write){
    write.table(dt, file=outfile, sep='\t', col.names=T, row.names=F, quote=F)
  }
  return(dt)
  
}


#' Title
#'
#' @param tissue_code 
#' @param meta 
#' @param counts 
#' @param covariates 
#' @param curr_outliers 
#' @param date 
#' @param write 
#'
#' @return TODO
#' 
#' @export
#' @import metap
#' @import DESeq2
#' @import data.table
#'
#' @examples
#' TODO
#' 
transcript_training_dea_each_sex = function(tissue_code, meta, counts, covariates, curr_outliers, date, write=T){
  
  # fix some inconsistencies 
  if(tissue_code == 't54-hypothalmus'){
    tissue_code = 't54-hypothalamus'
  }
  
  if(write){
    outfile = sprintf('dea/pass1b-06_%s_transcript-rna-seq_training-dea_%s.txt',
                      tissue_code,date)
    
    if(file.exists(outfile)){
      dt = fread(outfile, sep='\t', header=T)
      return(dt)
    }
  }
  
  # add vena cava outliers
  if(tissue_code=="t65-aorta"){
    curr_outliers = unique(c(curr_outliers, as.character(meta[sex=='female' & group %in% c('1w','2w'), viallabel])))
    meta = meta[,viallabel := as.character(viallabel)]
    meta = meta[!viallabel %in% curr_outliers]
    counts = counts[meta[,viallabel]]
  }
  
  sex_res = list()
  for(SEX in unique(meta[,sex])){
    
    # subset counts and meta
    curr_samples = meta[sex == SEX, viallabel]
    curr_meta = meta[sex == SEX]
    curr_counts = counts[,curr_samples]
    
    # center and scale continuous variables
    curr_cov = covariates
    for (cov in curr_cov){
      # remove if constant
      if(length(unique(curr_meta[,get(cov)])) == 1){
        message(sprintf("Covariate %s is constant for %s. Removing.", cov, SEX))
        curr_cov = curr_cov[curr_cov != cov]
      }else{
        # center and scale
        if(is.numeric(curr_meta[,get(cov)])){
          curr_meta[,(cov) := scale(curr_meta[,get(cov)], center = T, scale = T)]
        }
      }
    }
    
    full = paste0('~', paste0(c(curr_cov, 'group'), collapse=' + '))
    reduced = paste0('~', paste0(curr_cov, collapse=' + '))
    # # custom contrast for vena cava
    # if(tissue_code=="t65-aorta"){
    #   meta[,group := factor(group, levels=c('control','1w','2w','4w','8w'))]
    #   coldata = data.frame(model.matrix(eval(parse(text=contrast)), data=meta))
    #   coldata$`group1w.sexmale` = NULL
    #   coldata$`group2w.sexmale` = NULL
    #   coldata$X.Intercept. = NULL
    #   contrast = paste0('~', paste0(colnames(coldata), collapse=' + '))
    #   meta = coldata
    #   cols = colnames(meta)[grepl('group|sex', colnames(meta))]
    #   meta[cols] = sapply(meta[cols],as.factor)
    #   reduced = paste0('~', paste0(colnames(meta)[!grepl('group', colnames(meta))], collapse=' + '))
    # }
    
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
                        assay = 'transcript-rna-seq',
                        removed_samples = paste0(curr_outliers, collapse=','),
                        full_model=gsub(' ','',full),
                        reduced_model=gsub(' ','',reduced))
    res_dt = res_dt[,.(
      feature_ID,
      assay,
      tissue,
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
  
  if(write){
    write.table(res_dt, file=outfile, sep='\t', col.names=T, row.names=F, quote=F)
  }
  return(res_dt)
}
