#' Function to extract standard errors from limma results
#'
#' @param limma_res Result table as produced by limma
#' @param e_fit eBayes fit object as produced by limma
#' @param effect_col The column containing the effect size
#' @param t_col The column containing the t statistic
#'
#' @return A vector of standard errors for effect sizes
#' @export
#'
limma_res_extract_se<-function(limma_res,
                               e_fit,
                               effect_col="logFC",
                               t_col="t"){
  # First approach
  effects = limma_res[[effect_col]]
  ts = limma_res[[t_col]]
  ses1 = effects/ts
  if(is.null(colname)){
    return(ses1)
  }
}




## Proteomics Timewise DEA ------------------------------------------------------------------

#' Differential abundance analysis of proteomics, phosphoproteomics, acetylome, ubiquitylome
#'
#' @param assay_abbrev Abbreviation for proteomics assay to be analyzed as defined by ASSAY_ABBREV. 
#' One of the following: PROT, PHOSPHO, ACETYL, UBIQ
#' 
#' @param tissue_abbrev Abbreviation for proteomics tissue to be analyzed as defined by TISSUE_ABBREV. 
#' One of the following: PROT, PHOSPHO, ACETYL, UBIQ
#'
#' @return A data frame containing differential enrichment results
#' 
#' @import MotrpacRatTraining6moData
#' @import dplyr
#' @import limma
#' @import metap
#' @import data.table
#' @export
#'
#' @examples
#' proteomics_timewise_dea("PROT","HEART")
proteomics_timewise_dea  = function(assay_abbrev, tissue_abbrev){
  
  assay = assay_abbrev
  tissue = tissue_abbrev
  tpDA_split_sex = c() # keep the timewise results
  
  
  #Extract current dataset and metadata
  x = 
    get(sprintf("%s_%s_NORM_DATA", assay,gsub("-","",tissue))) %>%
    column_to_rownames(var = "feature_ID") %>%
    as.matrix()
  
  #Specify covariates
  covs <-  
    PHENO %>%
    filter(viallabel %in% colnames(x)) %>%
    transmute(
      sex = if_else(sex == "male","M","F"),
      group,
      tr = factor(case_when(
        group == "control" ~ "8_1",
        group == "1w" ~ "1_0",
        group == "2w" ~ "2_0",
        group == "4w" ~ "4_0",
        group == "8w" ~ "8_0",
      ),
      levels = c("8_1","1_0","2_0","4_0","8_0")
      )
    )
  
  covs <- covs[colnames(x),]
  
  #Variables that save sex-specific results
  sex_res = list() #F-test result
  sex_ttest_res = list() #T-test result
  
  for(SEX in c('M','F')){
    curr_meta = covs %>% filter(sex == SEX) 
    curr_counts = x[,rownames(curr_meta)]
    
    #Extract treatment covariate
    tr <- curr_meta$tr
    
    design = model.matrix(~0+tr)
    
    #Set contrasts
    cont.matrix = makeContrasts(
      tr1_0 - tr8_1, tr2_0 - tr8_1, tr4_0 - tr8_1, tr8_0 - tr8_1, 
      levels = design
    )
    
    colnames(cont.matrix) = c("1w","2w","4w","8w")
    
    # Fit the new model
    limma_model2 = lmFit(curr_counts,design)
    lmfit.cont <- contrasts.fit(limma_model2, cont.matrix)
    lmfit.cont.ebayes <- eBayes(lmfit.cont)
    
    #Extract results for each timepoint
    for(curr_tp in colnames(lmfit.cont.ebayes$t)){
      
      #Extract results
      limma_res = topTable(lmfit.cont.ebayes,
                           coef = curr_tp,number = Inf,sort.by = "none")
      
      curr_res = data.frame(
        feature_ID = rownames(limma_res),
        tissue=tissue,
        assay=assay,
        sex = if_else(SEX == "M","male","female"),
        logFC_se = limma_res_extract_se(limma_res,lmfit.cont.ebayes),
        logFC = limma_res$logFC,
        tscore = limma_res$t,
        covariates = NA,
        comparison_group = curr_tp,
        p_value = limma_res$P.Value
      )
      # Add group average intensities
      case_samps = covs %>% 
        filter(sex == SEX &
                 group == curr_tp) %>%
        rownames()
      curr_res$comparison_average_intensity = apply(x[,case_samps],1,mean,na.rm=TRUE)
      
      control_samps = covs %>% 
        filter(sex == SEX &
                 group == "control") %>%
        rownames()
      curr_res$reference_average_intensity = apply(x[,control_samps],1,mean,na.rm=TRUE)
      
      # Add NA counts
      curr_res$numNAs = rowSums(is.na(x[,c(case_samps,control_samps)]))
      
      # Add the results
      tpDA_split_sex = rbind(tpDA_split_sex,curr_res)
    }
  }
  
  return(tpDA_split_sex)
}










## Proteomics Training DEA ------------------------------------------------------------------

#' Title
#'
#' @param assay_abbrev Abbreviation for proteomics assay to be analyzed as defined by ASSAY_ABBREV. 
#' One of the following: PROT, PHOSPHO, ACETYL, UBIQ
#' 
#' @param tissue_abbrev Abbreviation for proteomics tissue to be analyzed as defined by TISSUE_ABBREV. 
#' One of the following: PROT, PHOSPHO, ACETYL, UBIQ
#'
#' @return A data frame containing differential enrichment results
#' 
#' @import MotrpacRatTraining6moData
#' @import dplyr
#' @import limma
#' @import metap
#' @import data.table
#' @export
#'
#' @examples
#' proteomics_training_dea("PROT","HEART")
proteomics_training_dea  <- function(assay_abbrev, tissue_abbrev){
  
  assay = assay_abbrev
  tissue = tissue_abbrev
  ftest_res_split_sex = c() # keeps all ftest results
  
  
  #Extract current dataset and metadata
  x = 
    get(sprintf("%s_%s_NORM_DATA", assay,gsub("-","",tissue))) %>%
    column_to_rownames(var = "feature_ID") %>%
    as.matrix()
  
  #Specify covariates
  covs <-  
    PHENO %>%
    filter(viallabel %in% colnames(x)) %>%
    transmute(
      sex = if_else(sex == "male","M","F"),
      group,
      tr = factor(case_when(
        group == "control" ~ "8_1",
        group == "1w" ~ "1_0",
        group == "2w" ~ "2_0",
        group == "4w" ~ "4_0",
        group == "8w" ~ "8_0",
      ),
      levels = c("8_1","1_0","2_0","4_0","8_0")
      )
    )
  
  covs <- covs[colnames(x),]
  
  #Variables that save sex-specific results
  sex_res = list() #F-test result
  sex_ttest_res = list() #T-test result
  
  for(SEX in c('M','F')){
    curr_meta = covs %>% filter(sex == SEX) 
    curr_counts = x[,rownames(curr_meta)]
    
    ###################################################################
    # F-test analysis - training-dea table
    
    #Extract treatment covariate
    tr <- curr_meta$tr
    #Generate the experimental model
    design <- model.matrix(~ 1+tr)
    fit <- lmFit(curr_counts, design)
    fit.eb <- eBayes(fit)
    
    
    #Extract results
    res = topTable(fit.eb, coef = 2:5, n=nrow(fit.eb))
    dt = data.table(tissue=tissue,
                    assay=assay,
                    feature_ID=rownames(res),
                    fscore=res$`F`,
                    p_value = res$P.Value,
                    full_model = "~1+group",
                    reduced_model = "~1")
    sex_res[[SEX]] = dt
    
  }
  
  #Merge F-test results
  merged = data.frame(merge(sex_res[['M']], sex_res[['F']], 
                            by=c("tissue","assay","feature_ID","full_model","reduced_model"), 
                            suffixes=c('_male','_female'))) %>%
    mutate(p_value_male = replace_na(p_value_male,1),
           p_value_female = replace_na(p_value_female,1)) %>%
    mutate(p_value = map2_dbl(p_value_male,p_value_female,function(x,y){sumlog(c(x,y))$p}))
  
  ftest_res_split_sex <- rbind(ftest_res_split_sex,merged)
  return(ftest_res_split_sex)
}
