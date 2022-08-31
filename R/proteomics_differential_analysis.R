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
limma_res_extract_se = function(limma_res,
                               e_fit,
                               effect_col="logFC",
                               t_col="t"){
  effects = limma_res[[effect_col]]
  ts = limma_res[[t_col]]
  ses1 = effects/ts
  return(ses1)
}


#' Proteomics timewise differential analysis
#' 
#' Timewise differential analysis for the proteome, phosphoproteome, acetylome, and ubiquitylome. 
#' Use limma to perform pairwise contrasts between each group of trained animals
#' and the sex-matched control group for a single tissue and proteomics assay. 
#' Analysis is performed separately for males and females. 
#'
#' @param assay character, abbreviation for proteomics assay to be analyzed as defined by [MotrpacRatTraining6moData::ASSAY_ABBREV]. 
#'   One of the following: PROT, PHOSPHO, ACETYL, UBIQ
#' @param tissue `r tissue()`
#' @param exclude_outliers bool, whether to remove sample outliers specified in 
#'   [MotrpacRatTraining6moData::OUTLIERS]. \code{TRUE} by default. 
#'
#' @return a data frame with one row per proteomics feature per contrast (usually 8 rows per gene):
#' \describe{
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{sex}}{`r sex()`}
#'   \item{\code{comparison_group}}{`r comparison_group()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{removed_samples}}{character, comma-separated list of outliers (vial labels) removed from differential analysis}
#'   \item{\code{logFC}}{`r logFC()`}
#'   \item{\code{logFC_se}}{`r logFC_se()`}
#'   \item{\code{tscore}}{double, t statistic}
#'   \item{\code{p_value}}{`r p_value_da()`}
#'   \item{\code{comparison_average_intensity}}{`r comparison_average_intensity()`}
#'   \item{\code{reference_average_intensity}}{`r reference_average_intensity()`}
#'   \item{\code{numNAs}}{`r numNAs()`}
#' }
#' 
#' @importFrom dplyr select filter transmute if_else case_when
#' @importFrom tibble column_to_rownames
#' @importFrom limma lmFit eBayes topTable makeContrasts contrasts.fit
#' @importFrom stats model.matrix
#' 
#' @export
#'
#' @examples
#' # Run timewise differential analysis for heart proteins
#' proteomics_timewise_da("PROT","HEART")
proteomics_timewise_da  = function(assay, tissue, exclude_outliers=TRUE){
  
  tpDA_split_sex = c() # keep the timewise results
  outliers = data.table::data.table(MotrpacRatTraining6moData::OUTLIERS)
  # data.table workaround 
  .tissue = tissue
  .assay = assay
  
  #Extract current dataset and metadata
  x = 
    load_sample_data(tissue, assay, exclude_outliers=exclude_outliers) %>%
    tibble::column_to_rownames(var = "feature_ID") %>%
    dplyr::select(-c("feature","tissue","assay")) %>%
    as.matrix()
  
  #Specify covariates
  covs <-  
    MotrpacRatTraining6moData::PHENO %>%
    dplyr::filter(viallabel %in% colnames(x)) %>%
    dplyr::transmute(
      sex = dplyr::if_else(sex == "male","M","F"),
      group,
      tr = factor(dplyr::case_when(
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
    curr_meta = covs %>% dplyr::filter(sex == SEX) 
    curr_counts = x[,rownames(curr_meta)]
    
    #Extract treatment covariate
    tr <- curr_meta$tr
    
    design = stats::model.matrix(~0+tr)
    
    #Set contrasts
    cont.matrix = limma::makeContrasts(
      tr1_0 - tr8_1, tr2_0 - tr8_1, tr4_0 - tr8_1, tr8_0 - tr8_1, 
      levels = design
    )
    
    colnames(cont.matrix) = c("1w","2w","4w","8w")
    
    # Fit the new model
    limma_model2 = limma::lmFit(curr_counts,design)
    lmfit.cont <- limma::contrasts.fit(limma_model2, cont.matrix)
    lmfit.cont.ebayes <- limma::eBayes(lmfit.cont)
    
    #Extract results for each timepoint
    for(curr_tp in colnames(lmfit.cont.ebayes$t)){
      
      #Extract results
      limma_res = limma::topTable(lmfit.cont.ebayes,
                                  coef = curr_tp,
                                  number = Inf,
                                  sort.by = "none")
      
      curr_res = data.frame(
        feature_ID = rownames(limma_res),
        tissue=tissue,
        assay=assay,
        sex = dplyr::if_else(SEX == "M","male","female"),
        logFC_se = limma_res_extract_se(limma_res,lmfit.cont.ebayes),
        logFC = limma_res$logFC,
        tscore = limma_res$t,
        comparison_group = curr_tp,
        p_value = limma_res$P.Value
      )
      # Add group average intensities
      case_samps = covs %>% 
        dplyr::filter(sex == SEX &
                 group == curr_tp) %>%
        rownames()
      curr_res$comparison_average_intensity = apply(x[,case_samps],1,mean,na.rm=TRUE)
      
      control_samps = covs %>% 
        dplyr::filter(sex == SEX &
                 group == "control") %>%
        rownames()
      curr_res$reference_average_intensity = apply(x[,control_samps],1,mean,na.rm=TRUE)
      
      # Add NA counts
      curr_res$numNAs = rowSums(is.na(x[,c(case_samps,control_samps)]))
      
      # Add outliers
      curr_outliers = NA_character_
      if(exclude_outliers){
        curr_outliers = outliers[tissue==.tissue & 
                                   assay==.assay & 
                                   grepl(sprintf("control|%s", curr_tp), group) & 
                                   grepl(sprintf("^%s", ifelse(SEX=="M", "male", "female")), group)]
        if(nrow(curr_outliers)>0){
          curr_outliers = paste(curr_outliers[,viallabel], collapse=",")
        }else{
          curr_outliers = NA_character_
        }
      }
      curr_res$removed_samples = curr_outliers
      
      # Add the results
      tpDA_split_sex = rbind(tpDA_split_sex,curr_res)
    }
  }
  
  # add some columns and reorder
  tpDA_split_sex$assay_code = MotrpacRatTraining6moData::ASSAY_ABBREV_TO_CODE[[assay]]
  tpDA_split_sex$tissue_code = MotrpacRatTraining6moData::TISSUE_ABBREV_TO_CODE[[tissue]]
  tpDA_split_sex = tpDA_split_sex[,c('feature_ID',
                                     'sex',
                                     'comparison_group',
                                     'assay',
                                     'assay_code',
                                     'tissue',
                                     'tissue_code',
                                     'removed_samples',
                                     'logFC',
                                     'logFC_se',
                                     'tscore',
                                     'p_value',
                                     'comparison_average_intensity',
                                     'reference_average_intensity',
                                     'numNAs')]
  
  return(tpDA_split_sex)
}


#' Proteomics training differential analysis 
#' 
#' Training differential analysis for the proteome, phosphoproteome, acetylome, and ubiquitylome.
#' Use limma to perform an F-test to test the effect of training
#' across time points. Analysis is performed separately for males and females. 
#'
#' @param assay character, abbreviation for proteomics assay to be analyzed as defined by [MotrpacRatTraining6moData::ASSAY_ABBREV]. 
#'   One of the following: PROT, PHOSPHO, ACETYL, UBIQ
#' @param tissue `r tissue()`
#' @param exclude_outliers bool, whether to remove sample outliers specified in
#'   [MotrpacRatTraining6moData::OUTLIERS]. \code{TRUE} by default. 
#'
#' @return a data frame with one row per proteomics feature:
#' \describe{
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{removed_samples}}{character, comma-separated list of outliers (vial labels) removed from differential analysis}
#'   \item{\code{fscore_male}}{double, F statistic for males}
#'   \item{\code{fscore_female}}{double, F statistic for females}
#'   \item{\code{p_value_male}}{double, nominal F-test p-value for males}
#'   \item{\code{p_value_female}}{double, nominal F-test p-value for females}
#'   \item{\code{full_model}}{character, full model used in F-test}
#'   \item{\code{reduced_model}}{character, reduced model used in F-test}
#'   \item{\code{p_value}}{double, combined male and female nominal p-value using the sum of logs}
#' }
#' 
#' @importFrom dplyr select filter transmute case_when if_else mutate
#' @importFrom tibble column_to_rownames
#' @importFrom purrr map2_dbl
#' @importFrom tidyr replace_na
#' @importFrom limma lmFit eBayes topTable
#' @importFrom metap sumlog
#' @importFrom stats model.matrix
#' 
#' @export
#'
#' @examples
#' # Run training differential analysis for heart proteins
#' proteomics_training_da("PROT","HEART")
proteomics_training_da = function(assay, tissue, exclude_outliers=TRUE){
  
  ftest_res_split_sex = c() # keeps all ftest results
  outliers = data.table::data.table(MotrpacRatTraining6moData::OUTLIERS)
  # data.table workaround 
  .tissue = tissue
  .assay = assay
  
  #Extract current dataset and metadata
  x = 
    load_sample_data(tissue, assay, exclude_outliers=exclude_outliers) %>%
    tibble::column_to_rownames(var = "feature_ID") %>%
    dplyr::select(-c("feature","tissue","assay")) %>%
    as.matrix()
  
  #Specify covariates
  covs <-  
    MotrpacRatTraining6moData::PHENO %>%
    dplyr::filter(viallabel %in% colnames(x)) %>%
    dplyr::transmute(
      sex = dplyr::if_else(sex == "male","M","F"),
      group,
      tr = factor(dplyr::case_when(
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
    curr_meta = covs %>% dplyr::filter(sex == SEX) 
    curr_counts = x[,rownames(curr_meta)]
    
    ###################################################################
    # F-test analysis - training-da table
    
    #Extract treatment covariate
    tr <- curr_meta$tr
    #Generate the experimental model
    design <- stats::model.matrix(~ 1+tr)
    fit <- limma::lmFit(curr_counts, design)
    fit.eb <- limma::eBayes(fit)
    
    #Extract results
    res = limma::topTable(fit.eb, coef = 2:5, n=nrow(fit.eb))
    dt = data.table::data.table(tissue=tissue,
                    assay=assay,
                    feature_ID=rownames(res),
                    fscore=res$`F`,
                    p_value = res$P.Value,
                    full_model = "~1+group",
                    reduced_model = "~1")

    # Add outliers
    curr_outliers = NA_character_
    if(exclude_outliers){
      curr_outliers = outliers[tissue==.tissue & 
                                 assay==.assay & 
                                 grepl(sprintf("^%s", ifelse(SEX=="M", "male", "female")), group)]
      if(nrow(curr_outliers)>0){
        curr_outliers = paste(curr_outliers[,viallabel], collapse=",")
      }else{
        curr_outliers = NA_character_
      }
    }
    dt[,removed_samples := curr_outliers]
    
    sex_res[[SEX]] = dt
    
  }
  
  # Merge F-test results
  merged = data.frame(merge(sex_res[['M']], sex_res[['F']], 
                            by=c("tissue","assay","feature_ID","full_model","reduced_model"), 
                            suffixes=c('_male','_female'))) %>%
    dplyr::mutate(p_value_male = tidyr::replace_na(p_value_male,1),
           p_value_female = tidyr::replace_na(p_value_female,1)) %>%
    dplyr::mutate(p_value = purrr::map2_dbl(p_value_male,p_value_female,function(x,y){metap::sumlog(c(x,y))$p}))
  
  ftest_res_split_sex <- rbind(ftest_res_split_sex,merged)
  
  # Add columns and change order
  ftest_res_split_sex$tissue_code = MotrpacRatTraining6moData::TISSUE_ABBREV_TO_CODE[[tissue]]
  ftest_res_split_sex$assay_code = MotrpacRatTraining6moData::ASSAY_ABBREV_TO_CODE[[assay]]
  ftest_res_split_sex = ftest_res_split_sex[,c(
    "feature_ID",
    "assay",
    "assay_code",
    "tissue",
    "tissue_code",
    "removed_samples_male",
    "removed_samples_female",
    "fscore_male",
    "fscore_female",
    "p_value_male",
    "p_value_female",
    "full_model",
    "reduced_model",
    "p_value"
  )]
  
  return(ftest_res_split_sex)
}
