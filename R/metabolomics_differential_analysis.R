#' Function to extract standard errors from limma results
#'
#' @param limma_res Result table as produced by limma
#' @param e_fit eBayes fit object as produced by limma
#' @param effect_col The column containing the effect size
#' @param t_col The column containing the t statistic
#'
#' @return A vector of standard errors for effect sizes
#' 
#' @export
#'
limma_res_extract_se<-function(limma_res,
                               e_fit,
                               effect_col="logFC",
                               t_col="t"){
  effects = limma_res[[effect_col]]
  ts = limma_res[[t_col]]
  ses1 = effects/ts
  return(ses1)
}

#' Metabolomics Training Differential Analysis
#' 
#' Use limma to test the effect of training across timepoints. The analysis
#' is performed separately for Male and Female rats and combined into a single
#' p-values using sum of logs.
#'
#' @param tissue_abbrev Abbreviation for proteomics tissue to be analyzed as defined by TISSUE_ABBREV.
#'
#' @return a data frame with one row per metabolite:
#' \describe{
#'   \item{\code{feature_ID}}{Metabolite name}
#'   \item{\code{dataset}}{The metabolomics assaay in which the metabolite is detected.}
#'   \item{\code{groups_tested_female}}{The timepoints used to perform the F-test in females. Some tissues or assays are missing timepoints}
#'   \item{\code{groups_tested_male}}{The timepoints used to perform the F-test in males. Some tissues or assays are missing timepoints}
#'   \item{\code{fscore_male}}{F-statistic for males}
#'   \item{\code{fscore_female}}{F-statistic for females}
#'   \item{\code{p_value_male}}{nominal p-value for males}
#'   \item{\code{p_value_female}}{nominal  p-value for females}
#'   \item{\code{full_model}}{full model used in limma}
#'   \item{\code{reduced_model}}{Representation of the reduced model, although not used by limma}
#'   \item{\code{p_value}}{combined male and female nominal p-value using the sum of logs}
#' }
#' 
#' @export
#' @import MotrpacRatTraining6moData
#' @import dplyr
#' @import limma
#' @import metap
#' @import data.table
#' @import magrittr
#' 
#' @examples
#' metab_training_dea("HEART")
metab_training_dea <- function(tissue_abbrev){
  
  tissue = tissue_abbrev
  
  #Save F-test results
  f_results = list()
  
  for(metab_assay in names(MotrpacRatTraining6moData::METAB_NORM_DATA_NESTED)){
      
    #Extract current dataset and metadata
    x = MotrpacRatTraining6moData::METAB_NORM_DATA_NESTED[[metab_assay]][[tissue]]
    if(is.null(x)){
      next
    }
    
    x <- as.matrix(x)
    
    
    #Specify covariates
    covs <-  
      MotrpacRatTraining6moData::PHENO %>%
      dplyr::filter(viallabel %in% colnames(x)) %>%
      dplyr::transmute(
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
    
    #Variables that save sex-specific results
    sex_res = list() #F-test result
    
    for(SEX in c('M','F')){
      curr_meta = covs %>% filter(sex == SEX) 
      if(nrow(curr_meta) > 0){
        curr_counts = x[,rownames(curr_meta)]
        
        #Extract treatment covariate
        tr <- factor(curr_meta$tr,
                     levels= c("8_1",
                               intersect(c("1_0","2_0","4_0","8_0"),
                                         unique(curr_meta$tr))))
        
        #Generate the experimental model
        design <- model.matrix(~ 1+tr)
        fit <- limma::lmFit(curr_counts, design)
        fit.eb <- limma::eBayes(fit)
        
        #Extract results
        id_coefs <- colnames(design)[2:ncol(design)]
        res = limma::topTable(fit.eb, coef = id_coefs, n=nrow(fit.eb))
        dt = data.table(feature_ID=rownames(res),
                        fscore=res$`F`,
                        p_value = res$P.Value,
                        full_model = "~1+group",
                        reduced_model = "~1",
                        groups_tested = paste(id_coefs,collapse=";"))
        sex_res[[SEX]] = dt
      } else {
        #Handle cases in which there are no samples for one sex
        curr_meta = covs %>% filter(sex != SEX) 
        curr_counts = x[,rownames(curr_meta)]
        dt = data.table(feature_ID=rownames(curr_counts),
                        fscore=NA,
                        p_value = NA,
                        full_model = "~1+group",
                        reduced_model = "~1",
                        groups_tested = NA)
        sex_res[[SEX]] = dt
      }
      
    }
    
    #Merge F-test results
    
    merged <- data.frame(merge(sex_res[['M']], sex_res[['F']], 
                               by=c("feature_ID","full_model","reduced_model"), 
                               suffixes=c('_male','_female'))) 
    
    #Special handling for ovaries and testes which only have samples for one sex
    if(sum(is.na(merged$p_value_female)) == nrow(merged)){
      merged <- merged %>%
        mutate(p_value_male = replace_na(p_value_male,1),
               p_value = p_value_male)
    } else if (sum(is.na(merged$p_value_male)) == nrow(merged)){
      merged <- merged %>% 
        mutate(p_value_female = replace_na(p_value_female,1),
               p_value = p_value_female)
    } else {
      merged <- merged %>%
        mutate(p_value_male = replace_na(p_value_male,1),
               p_value_female = replace_na(p_value_female,1)) %>%
        mutate(p_value = map2_dbl(p_value_male,p_value_female,function(x,y){sumlog(c(x,y))$p}))
    }
    
    merged$dataset = metab_assay
    
    f_results = c(f_results,list(merged))
  }
  return(dplyr::bind_rows(f_results))
}



#' Metabolomics Training Differential Analysis
#' 
#' Use limma to test the effect of training for each exercised time point vs the control. The analysis
#' is performed separately for Male and Female rats 
#'
#' @param tissue_abbrev Abbreviation for proteomics tissue to be analyzed as defined by TISSUE_ABBREV.
#'
#' @return a data frame with one row per metabolite, time, and sex combination:
#' \describe{
#'   \item{\code{feature_ID}}{Metabolite name}
#'   \item{\code{sex}}{one of "male" or "female"}
#'   \item{\code{cv}}{Coefficient of variation for case sample}
#'   \item{\code{control_cv}}{Coefficient of variation for control samples}
#'   \item{\code{logFC}}{log fold-change of the training group specified by \code{sex} and \code{comparison_group} (e.g., 1-week females) 
#'     relative to the sex-matched sedentary controls}
#'   \item{\code{logFC_se}}{standard error of \code{logFC}}
#'   \item{\code{tscore}}{t-statistic of the training group specified by \code{sex} and \code{comparison_group} (e.g., 1-week females) 
#'     relative to the sex-matched sedentary controls}
#'   \item{\code{covariates}}{Additional covariates used in the analysis. NA for none.}
#'   \item{\code{comparison_group}}{time point compared to the sex-matched 
#'       sedentary control animals, one of "1w", "2w", "4w", "8w"}
#'   \item{\code{p_value}}{nominal p-value corresponding to the contrast between the training group 
#'     (e.g., 1-week females) and the sex-matched sedentary controls}
#'   \item{\code{comparison_average_intensity}}{average normalized RNA-seq counts for samples in the training group (e.g., 1-week females)}
#'   \item{\code{reference_average_intensity}}{average normalized RNA-seq counts for sex-matched sedentary control samples}
#'   \item{\code{dataset}}{The metabolomics assaay in which the metabolite is detected.}
#'}
#' 
#' @export
#' @import MotrpacRatTraining6moData
#' @import dplyr
#' @import limma
#' @import metap
#' @import data.table
#' @import magrittr
#' 
#' @examples
#' # Perform differential analysis for metabolites in heart tissue.
#' metab_training_dea("HEART")
metab_timewise_dea <- function(tissue_abbrev){
  
  tissue = tissue_abbrev
  
  #Save F-test results
  t_results = list()
  
  for(metab_assay in names(MotrpacRatTraining6moData::METAB_NORM_DATA_NESTED)){
    
    #Extract current dataset and metadata
    x = MotrpacRatTraining6moData::METAB_NORM_DATA_NESTED[[metab_assay]][[tissue]]
    if(is.null(x)){
      next
    }
    
    x <- as.matrix(x)
    
    
    #Specify covariates
    covs <-  
      MotrpacRatTraining6moData::PHENO %>%
      dplyr::filter(viallabel %in% colnames(x)) %>%
      dplyr::transmute(
        sex = if_else(sex == "male","M","F"),
        group,
        tr = factor(case_when(
          group == "control" ~ "8_1",
          group == "1w" ~ "1_0",
          group == "2w" ~ "2_0",
          group == "4w" ~ "4_0",
          group == "8w" ~ "8_0",
        ),
        levels = c("8_1","1_0","2_0","4_0","8_0")),
        is_control = dplyr::if_else(tr == "8_1", 1, 0),
        viallabel
      )
    
    
    sex_ttest_res = list() #T-test result
    
    #Calculate CV
    x_cv <- apply(x,1,sd)/apply(x,1,mean)
    
    #Calculate CV for control
    control_viallabels <- covs %>% filter(is_control == 1) %>% .$viallabel
    x_control = x[,control_viallabels]
    if(length(control_viallabels)>3){
      control_cvs = apply(x_control,1,sd)/apply(x_control,1,mean)
    }else{
      control_cvs = rep(NA,length(x_cvs))
    }
    
    for(SEX in c('M','F')){
      curr_meta = covs %>% filter(sex == SEX) 
      
      if(nrow(curr_meta) > 0){
        
        curr_counts = x[,rownames(curr_meta)]
        
        
        #Extract treatment covariate
        tr <- factor(curr_meta$tr,
                     levels= c("8_1",
                               intersect(c("1_0","2_0","4_0","8_0"),
                                         unique(curr_meta$tr))))
        
        design = model.matrix(~0+tr)
        
        #Set contrasts
        cont.matrix = limma::makeContrasts(
          tr1_0 - tr8_1, tr2_0 - tr8_1, tr4_0 - tr8_1, tr8_0 - tr8_1, 
          levels = c("tr8_1","tr1_0","tr2_0","tr4_0","tr8_0")
        )
        
        #Remove contrasts for which samples are not available in a group
        cont.matrix <- cont.matrix[rownames(cont.matrix) %in% 
                                     paste0("tr",levels(tr)),]
        
        cont.matrix <- cont.matrix[,apply(cont.matrix,2,sum)== 0]
        
        #Rename column of contrast matrix to something readable
        colnames(cont.matrix) = case_when(
          colnames(cont.matrix) == "tr1_0 - tr8_1" ~ "1w",
          colnames(cont.matrix) == "tr2_0 - tr8_1" ~ "2w",
          colnames(cont.matrix) == "tr4_0 - tr8_1" ~ "4w",
          colnames(cont.matrix) == "tr8_0 - tr8_1" ~ "8w",
        )
        
        # Fit the new model
        limma_model2 = limma::lmFit(curr_counts,design)
        lmfit.cont <- limma::contrasts.fit(limma_model2, cont.matrix)
        lmfit.cont.ebayes <- limma::eBayes(lmfit.cont)
        
        #Extract results for each timepoint
        for(curr_tp in colnames(lmfit.cont.ebayes$t)){
          
          #Extract results
          limma_res = limma::topTable(lmfit.cont.ebayes,
                               coef = curr_tp,number = Inf,sort.by = "none")
          
          curr_res = data.frame(
            feature_ID = rownames(limma_res),
            sex = if_else(SEX == "M","male","female"),
            cv = x_cv,
            control_cv = control_cvs, 
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
          
          curr_res$comparison_average_intensity = apply(data.frame(x[,case_samps]),1,mean,na.rm=TRUE)
          
          
          control_samps = covs %>% 
            filter(sex == SEX &
                     is_control == 1) %>%
            rownames()
          curr_res$reference_average_intensity = apply(x[,control_samps],1,mean,na.rm=TRUE)
          
          
          # Add the results
          sex_ttest_res = rbind(sex_ttest_res,curr_res)
        }
      }
    }
    sex_ttest_res$dataset = metab_assay
    t_results = c(t_results,list(sex_ttest_res))
  }
  return(dplyr::bind_rows(t_results))
}
