#' Metabolomics Training Differential Analysis
#' 
#' Use limma to test the effect of training across timepoints. The analysis
#' is performed separately for Male and Female rats and combined into a single
#' p-values using sum of logs.
#'
#' @param tissue `r tissue()`
#'
#' @return a data frame with one row per metabolite:
#' \describe{
#'   \item{\code{feature_ID}}{Metabolite name}
#'   \item{\code{platform}}{The metabolomics assaay in which the metabolite is detected.}
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
#' @importFrom dplyr filter transmute bind_rows if_else case_when mutate
#' @importFrom tidyr replace_na
#' @importFrom limma lmFit eBayes topTable 
#' @importFrom purrr map2_dbl
#' @import metap sumlog
#' 
#' @examples
#' metab_training_da("HEART")
metab_training_da <- function(tissue){
  
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
    
    #Variables that save sex-specific results
    sex_res = list() #F-test result
    
    for(SEX in c('M','F')){
      curr_meta = covs %>% dplyr::filter(sex == SEX) 
      if(nrow(curr_meta) > 0){
        curr_counts = x[,rownames(curr_meta)]
        
        #Extract treatment covariate
        tr <- factor(curr_meta$tr,
                     levels= c("8_1",
                               intersect(c("1_0","2_0","4_0","8_0"),
                                         unique(curr_meta$tr))))
        
        #Generate the experimental model
        design <- stats::model.matrix(~ 1+tr)
        fit <- limma::lmFit(curr_counts, design)
        fit.eb <- limma::eBayes(fit)
        
        #Extract results
        id_coefs <- colnames(design)[2:ncol(design)]
        res = limma::topTable(fit.eb, coef = id_coefs, n=nrow(fit.eb))
        dt = data.table::data.table(feature_ID=rownames(res),
                        fscore=res$`F`,
                        p_value = res$P.Value,
                        full_model = "~1+group",
                        reduced_model = "~1",
                        groups_tested = paste(id_coefs,collapse=";"))
        sex_res[[SEX]] = dt
      } else {
        #Handle cases in which there are no samples for one sex
        curr_meta = covs %>% dplyr::filter(sex != SEX) 
        curr_counts = x[,rownames(curr_meta)]
        dt = data.table::data.table(feature_ID=rownames(curr_counts),
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
        dplyr::mutate(p_value_male = tidyr::replace_na(p_value_male,1),
               p_value = p_value_male)
    } else if (sum(is.na(merged$p_value_male)) == nrow(merged)){
      merged <- merged %>% 
        dplyr::mutate(p_value_female = tidyr::replace_na(p_value_female,1),
               p_value = p_value_female)
    } else {
      merged <- merged %>%
        dplyr::mutate(p_value_male = tidyr::replace_na(p_value_male,1),
               p_value_female = tidyr::replace_na(p_value_female,1)) %>%
        dplyr::mutate(p_value = purrr::map2_dbl(p_value_male,p_value_female,function(x,y){metap::sumlog(c(x,y))$p}))
    }
    
    merged$platform = metab_assay
    
    f_results = c(f_results,list(merged))
  }
  return(dplyr::bind_rows(f_results))
}



#' Metabolomics Training Differential Analysis
#' 
#' Use limma to test the effect of training for each exercised time point vs the control. The analysis
#' is performed separately for Male and Female rats 
#'
#'  @param tissue `r tissue()`
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
#'   \item{\code{platform}}{The metabolomics platform in which the metabolite is detected.}
#'}
#' 
#' @export
#' @import dplyr filter transmute bind_rows if_else case_when mutate
#' @importFrom limma lmFit eBayes topTable makeContrasts contrasts.fit
#' 
#' @examples
#' # Perform differential analysis for metabolites in heart tissue.
#' metab_training_da("HEART")
metab_timewise_da <- function(tissue){
  
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
        sex = dplyr::if_else(sex == "male","M","F"),
        group,
        tr = factor(dplyr::case_when(
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
    control_viallabels <- covs %>% dplyr::filter(is_control == 1) %>% .$viallabel
    x_control = x[,control_viallabels]
    if(length(control_viallabels)>3){
      control_cvs = apply(x_control,1,sd)/apply(x_control,1,mean)
    }else{
      control_cvs = rep(NA,length(x_cvs))
    }
    
    for(SEX in c('M','F')){
      curr_meta = covs %>% dplyr::filter(sex == SEX) 
      
      if(nrow(curr_meta) > 0){
        
        curr_counts = x[,rownames(curr_meta)]
        
        
        #Extract treatment covariate
        tr <- factor(curr_meta$tr,
                     levels= c("8_1",
                               intersect(c("1_0","2_0","4_0","8_0"),
                                         unique(curr_meta$tr))))
        
        design = stats::model.matrix(~0+tr)
        
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
        colnames(cont.matrix) = dplyr::case_when(
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
            sex = dplyr::if_else(SEX == "M","male","female"),
            cv = x_cv,
            control_cv = control_cvs, 
            logFC_se = limma_res_extract_se(limma_res),
            logFC = limma_res$logFC,
            tscore = limma_res$t,
            covariates = NA,
            comparison_group = curr_tp,
            p_value = limma_res$P.Value
            
          )
          # Add group average intensities
          case_samps = covs %>% 
            dplyr::filter(sex == SEX &
                     group == curr_tp) %>%
            rownames()
          
          curr_res$comparison_average_intensity = apply(data.frame(x[,case_samps]),1,mean,na.rm=TRUE)
          
          
          control_samps = covs %>% 
            dplyr::filter(sex == SEX &
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
