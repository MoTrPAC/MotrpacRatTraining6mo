# TODO
atac_training_dea_each_sex = function(filt_counts, meta, covariates, curr_outliers=c(), label=NULL, return_resid=F){
  
  meta_df = as.data.frame(meta, check.names=F)
  if("sample_batch" %in% tolower(covariates)){
    which = covariates[grepl("sample_batch", covariates, ignore.case=T)]
    meta_df[,which] = factor(meta_df[,which])
  }
  meta_df$group = factor(meta_df$group, levels=c('control','1w','2w','4w','8w')) # IMPORTANT
  
  # remove outliers 
  if(length(curr_outliers) > 0){
    meta_df = meta_df[!as.character(meta_df$viallabel) %in% as.character(curr_outliers),]
  }
  filt_counts = filt_counts[,as.character(meta_df$viallabel)]
  
  # split by sex 
  sex_res = list()
  resids = list()
  for(SEX in c('male','female')){
    
    curr_meta = meta_df[meta_df$sex==SEX,]
    curr_counts = filt_counts[,curr_meta$viallabel]
    
    full = paste0(c("~ 1", "group", covariates), collapse=" + ")
    reduced = paste0(c("~ 1", covariates), collapse=" + ")
    
    # normalize and get voom weights 
    design = model.matrix(eval(parse(text=full)), data = curr_meta)
    # check if full rank
    if(!is.fullrank(design)){
      if("sample_batch" %in% tolower(covariates)){
        which = covariates[grepl("sample_batch", covariates, ignore.case=T)]
        curr_cov = covariates[!covariates == which]
        full = paste0(c("~ 1", "group", curr_cov), collapse=" + ")
        reduced = paste0(c("~ 1", curr_cov), collapse=" + ")
        design = model.matrix(eval(parse(text=full)), data = curr_meta)
        warning(sprintf("Sample_batch and group or sex are confounded for %s %s", label, SEX))
      }else{
        stop(sprintf("Model matrix with design %s is not full rank.", full))
      }
    }else{
      curr_cov = covariates
    }
    curr_voom = voom(curr_counts, design=design, normalize.method="quantile")
    
    limma_model1 = lmFit(curr_voom, design)
    eb_Ftest = eBayes(limma_model1)
    res = topTable(eb_Ftest, n=Inf, coef=colnames(design)[grepl('group',colnames(design))], sort.by = "none")
    dt = data.table(assay='epigen-atac-seq',
                    dataset=label, 
                    tissue=gsub(",.*","",label),
                    feature_ID=rownames(res),
                    fscore=res$`F`,
                    p_value = res$P.Value,
                    adj_p_value = p.adjust(res$P.Value, method='BH'),
                    full_model=gsub(' ','',full),
                    reduced_model=gsub(' ','',reduced),
                    removed_samples=paste0(curr_outliers, collapse=','),
                    covariates=paste0(curr_cov, collapse=','))
    sex_res[[SEX]] = dt
    
    # residuals
    residual_mat = residuals(limma_model1, curr_voom)
    resids[[SEX]] = residual_mat
  }
  
  # merge 
  merged = data.table(merge(sex_res[['male']], sex_res[['female']], 
                            by=c("assay","tissue","feature_ID","removed_samples"), 
                            suffixes=c('_male','_female')))
  
  # get a single meta p-value per feature using the male- and female- specific p-values 
  merged[,p_value := sumlog(c(p_value_male, p_value_female))$p, by=seq_len(nrow(merged))]
  merged[,adj_p_value := p.adjust(p_value, method='BH')]
  
  if(!return_resid){
    return(merged)
  }
  
  residual_mat = data.frame(cbind(resids[['male']], resids[['female']]))
  
  return(list(res=merged, 
              residuals=residual_mat))
}


# TODO
atac_timewise_dea_each_sex = function(filt_counts, meta, covariates, curr_outliers=c(), label=NULL){
  
  meta_df = as.data.frame(meta, check.names=F)
  if("sample_batch" %in% tolower(covariates)){
    which = covariates[grepl("sample_batch", covariates, ignore.case=T)]
    meta_df[,which] = factor(meta_df[,which])
  }
  # remove outliers 
  if(length(curr_outliers) > 0){
    meta_df = meta_df[!as.character(meta_df$viallabel) %in% as.character(curr_outliers),]
  }
  filt_counts = filt_counts[,as.character(meta_df$viallabel)]
  
  # split by sex 
  sex_res = list()
  for(SEX in c('male','female')){
    
    curr_meta = meta_df[meta_df$sex==SEX,]
    curr_counts = filt_counts[,curr_meta$viallabel]
    
    full = paste0(c("~ 0", "group", covariates), collapse=" + ")
    reduced = paste0(c("~ 0", covariates), collapse=" + ")
    
    # normalize and get voom weights 
    design = model.matrix(eval(parse(text=full)), data = curr_meta)
    # check if full rank
    if(!is.fullrank(design)){
      if("sample_batch" %in% tolower(covariates)){
        which = covariates[grepl("sample_batch", covariates, ignore.case=T)]
        curr_cov = covariates[!covariates == which]
        full = paste0(c("~ 0", "group", curr_cov), collapse=" + ")
        reduced = paste0(c("~ 0", curr_cov), collapse=" + ")
        design = model.matrix(eval(parse(text=full)), data = curr_meta)
        warning(sprintf("Sample_batch and group or sex are confounded for %s %s", label, SEX))
      }else{
        stop(sprintf("Model matrix with design %s is not full rank.", full))
      }
    }else{
      curr_cov = covariates
    }
    curr_voom = voom(curr_counts, design=design, normalize.method="quantile")
    fit = lmFit(curr_voom,design)
    
    cont.matrix=makeContrasts(
      '1W'='group1w - groupcontrol',
      '2W'='group2w - groupcontrol',
      '4W'='group4w - groupcontrol',
      '8W'='group8w - groupcontrol',
      levels=design)
    
    fit2=contrasts.fit(fit,cont.matrix)
    e=eBayes(fit2)
    
    dea_list = list()
    for(tp in c('1W','2W','4W','8W')){
      res = topTable(e, number=nrow(e), coef=tp, confint=TRUE)
      dt = data.table(assay='epigen-atac-seq',
                      dataset=label, 
                      tissue=gsub(",.*","",label),
                      feature_ID=rownames(res),
                      sex=SEX,
                      comparison_group=tolower(tp),
                      logFC=res$logFC,
                      logFC_se=(res$CI.R - res$CI.L)/3.92,
                      tscore=res$t,
                      p_value = res$P.Value,
                      removed_samples=paste0(curr_outliers, collapse=','),
                      covariates=paste0(curr_cov, collapse=','))
      
      dea_list[[tp]] = dt
    }
    sex_res[[SEX]] = rbindlist(dea_list)
  }
  res = rbindlist(sex_res)
  res[,adj_p_value := p.adjust(p_value, method='BY')]
  return(res)
}
