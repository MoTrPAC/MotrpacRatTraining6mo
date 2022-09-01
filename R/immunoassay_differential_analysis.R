#' Immunoassay timewise differential analysis 
#' 
#' For each panel and tissue, perform pairwise contrasts between each group of 
#' trained animals and the sex-matched control group. 
#' Analysis is performed separately for males and females.
#'
#' @return a data frame with one row per immunoassay feature per contrast (usually 8 rows per feature):
#' \describe{
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{sex}}{`r sex()`}
#'   \item{\code{comparison_group}}{`r comparison_group()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{dataset}}{character, LUMINEX panel}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{covariates}}{character, comma-separated list of adjustment variables}
#'   \item{\code{logFC}}{`r logFC()`}
#'   \item{\code{logFC_se}}{`r logFC_se()`}
#'   \item{\code{p_value}}{`r p_value_da()`}
#'   \item{\code{comparison_average_intensity}}{`r comparison_average_intensity()`}
#'   \item{\code{reference_average_intensity}}{`r reference_average_intensity()`}
#' }
#' @export
#' 
#' @importFrom multcomp glht adjusted
#' @importFrom dplyr as_tibble mutate filter tibble arrange relocate rename
#' @importFrom tibble column_to_rownames
#'
#' @examples
#' res = immuno_timewise_da()
immuno_timewise_da = function(){
  
  data = fetch_object("IMMUNO_NORM_DATA_NESTED")
  meta = fetch_object("IMMUNO_META")
  meta = meta[,c("sex", "group", "tissue", "log2_CHEX4", "panel_name", "viallabel")]
  
  alldat = lapply(names(data), function(mydataset){

    rv2 = lapply(names(data[[mydataset]]), function(thisTissue){

      thisDat = data[[mydataset]][[thisTissue]]

      log2_mfi = thisDat %>% tibble::column_to_rownames("viallabel")
      mymetadata = meta[meta$tissue == thisTissue & meta$panel_name == mydataset,]
      
      # make sure viallabels match
      rownames(mymetadata) = mymetadata$viallabel
      mymetadata = mymetadata[rownames(log2_mfi),]
      stopifnot(all(rownames(mymetadata) == rownames(log2_mfi)))
      stopifnot(nrow(mymetadata) == nrow(log2_mfi))

      tmpans = lapply(1:ncol(log2_mfi), function(i){
        ans = dplyr::as_tibble(mymetadata[,c("group","sex","tissue","log2_CHEX4")])
        ans$mfi = log2_mfi[,i]
        ans$cytokine = colnames(log2_mfi)[i]
        ans$assay = mydataset
        ans
      })
      rv = do.call(rbind, tmpans)
      return(rv)
    })

    rv2 = do.call(rbind, rv2)
    return(rv2)
  })
  
  testdat = dplyr::as_tibble(do.call(rbind, alldat)) %>% 
    dplyr::mutate(cytokine_by_assay=paste(assay,cytokine,sep=":"))
  
  testdat$group = factor(testdat$group)
  testdat$group = stats::relevel(testdat$group,"control")
  testdat$sex = factor(testdat$sex)
  testdat$sex = stats::relevel(testdat$sex, "female")
  
  testdat2 = testdat %>% dplyr::mutate(test_grouping=paste(assay,cytokine,tissue,sex,sep=":"))

  testresList = lapply(unique(testdat2$test_grouping), function(thisTestGrouping){
      
      subdat = testdat2 %>% dplyr::filter(test_grouping == thisTestGrouping)
      
      # calculate averages before removing groups with only 1 replicate 
      myavg = tapply(subdat$mfi, subdat$group, mean)
      
      # check to see if there are at least two non-na samples per group
      groups_to_discard = c()
      group_table = table(subdat$group)
      if(any(group_table < 2)){
        groups_to_discard = names(group_table)[group_table < 2]
        if("control" %in% groups_to_discard){
          # if the timepoint is control, discard the whole analyte
          message(sprintf("Group 'control' in %s has less than 2 values. Removing analyte.", 
                          thisTestGrouping))
          return()
        }else{
          # if not, discard that timepoint
          message(sprintf("Group %s in %s has less than 2 values. Removing remaining samples in this group.",
                          paste0(groups_to_discard, collapse=', '), thisTestGrouping))
          subdat = subdat[!subdat$group %in% groups_to_discard,]
        }
      }
      
      # ANOVA 
      testaov = stats::aov(mfi ~ 0+ group + log2_CHEX4 , data = subdat)
      
      # Simultaneous Tests for General Linear Hypotheses - no multiple comparisons
      # only write contrasts for remaining groups after removing single-sample groups 
      groups_to_test = sort(unique(as.character(subdat$group[subdat$group!='control'])))
      contrast_strings = unname(unlist(sapply(groups_to_test, function(g){
        return(sprintf("group%s - groupcontrol = 0", g))
      })))
      glhtres = multcomp::glht(testaov, linfct = contrast_strings)
      
      # # Multiple Comparisons of Means - Dunnett Contrasts:
      # dunnettres = glht(testaov, linfct = mcp(group="Dunnett"))
      
      glht_unadj = summary(glhtres, test=multcomp::adjusted("none"))
      
      myansglht = dplyr::tibble(
        test_grouping = thisTestGrouping,
        terms = names(glht_unadj$test$tstat),
        coef = glht_unadj$test$coefficients,
        coef_se = glht_unadj$test$sigma,
        pval = as.vector(glht_unadj$test$pvalues),
        # padj = as.vector(summary(dunnettres)$test$pvalues), # Dunnett adjusted p-values 
        reference_average_intensity = myavg["control"],
        comparison_average_intensity = myavg[groups_to_test]
      )
      # add "NA" for discarded groups
      if(length(groups_to_discard)>0){
        filler_terms = unname(unlist(sapply(groups_to_discard, function(g){
          return(sprintf("group%s - groupcontrol", g))
        })))
        filler = dplyr::tibble(
          test_grouping = thisTestGrouping,
          terms = filler_terms,
          coef = NA_real_,
          coef_se = NA_real_,
          pval = NA_real_,
          reference_average_intensity = myavg["control"],
          comparison_average_intensity = myavg[groups_to_discard]
        )
        myansglht = rbind(myansglht, filler)
        # sort by group
        myansglht = dplyr::arrange(myansglht, by=terms)
      }
      
      return(list(glhtres=myansglht))
    })
  
  glhtres = do.call(rbind, lapply(testresList, "[[","glhtres"))
  
  reportTab = glhtres %>% dplyr::mutate(
    panel = sapply(strsplit(test_grouping,":"),"[",1),
    feature_ID = sapply(strsplit(test_grouping,":"),"[",2),
    assay = "immunoassay",
    tissue = sapply(strsplit(test_grouping,":"),"[",3),
    sex = sapply(strsplit(test_grouping,":"),"[",4),
    covariates = "log2_CHEX4") %>% 
    dplyr::rename(logFC=coef, logFC_se=coef_se, p_value=pval, comparison_group=terms) %>% 
    dplyr::select(-test_grouping) %>% 
    dplyr::relocate(feature_ID,panel,assay,tissue)
  
  reportTab = reportTab %>% dplyr::mutate(comparison_group=sub("^group(.+?) .*","\\1",comparison_group)) %>% as.data.frame()
  
  # reorder columns for consistency
  reportTab$tissue_code = MotrpacRatTraining6moData::TISSUE_ABBREV_TO_CODE[reportTab$tissue]
  reportTab$assay_code = reportTab$assay
  reportTab$assay = MotrpacRatTraining6moData::ASSAY_CODE_TO_ABBREV[reportTab$assay_code]
  reportTab$dataset = reportTab$panel
  reportTab = reportTab[,c(
    "feature_ID",
    "sex",
    "comparison_group",
    "assay",
    "assay_code",
    "dataset", 
    "tissue",
    "tissue_code",
    "covariates",
    'logFC',
    'logFC_se',
    'p_value',
    'comparison_average_intensity',
    'reference_average_intensity'
  )]
  
  return(reportTab)
}


#' Immunoassay training differential analysis 
#' 
#' For each panel and tissue, perform a likelihood ratio test to test the effect of training
#' across time points. Analysis is performed separately for males and females. 
#' 
#' @return a data frame with one row per feature:
#' \describe{
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{dataset}}{character, LUMINEX panel}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{removed_samples_male}}{character, comma-separated list of male samples (vial labels) removed from differential analysis}
#'   \item{\code{removed_samples_female}}{character, comma-separated list of female samples (vial labels) removed from differential analysis}
#'   \item{\code{lrt_male}}{double, likelihood ratio test statistic for males}
#'   \item{\code{lrt_female}}{double, likelihood ratio test statistic for females}
#'   \item{\code{p_value_male}}{double, nominal LRT p-value for males}
#'   \item{\code{p_value_female}}{double, nominal LRT p-value for females}
#'   \item{\code{full_model}}{character, full model used in LRT for males and females}
#'   \item{\code{reduced_model}}{character, reduced model used in LRT for males and females}
#'   \item{\code{p_value}}{double, combined male and female nominal p-value using the sum of logs}
#' }
#' @export
#' 
#' @importFrom lmtest lrtest
#'
#' @examples
#' res = immuno_training_da()
immuno_training_da = function(){
  
  imputed_data = fetch_object("IMMUNO_NORM_DATA_NESTED")
  metadata = data.table::data.table(fetch_object("IMMUNO_META"))[,c("sex", "group", "tissue", "log2_CHEX4", "panel_name", "viallabel"), with=F]
  
  training_res = list()
  i=1
  for(PANEL in names(imputed_data)){
    for(TISSUE in names(imputed_data[[PANEL]])){
      
      data = imputed_data[[PANEL]][[TISSUE]]
      all_meta = metadata[tissue == TISSUE & panel_name == PANEL]
      
      # melt data
      mdata = as.data.frame(data.table::melt(as.data.table(data), 
                                             id.vars=c('viallabel'),
                                             variable.name='feature_ID', 
                                             value.name='value'))
      # merge with meta
      mdata = merge(mdata, all_meta, by="viallabel")
      
      # for each analyte, perform sex-split tests
      for(analyte in unique(mdata$feature_ID)){
        training_sex = list()
        for(SEX in unique(mdata$sex)){
          curr_data = mdata[mdata$feature_ID == analyte & mdata$sex == SEX,]
          # check to see if there are at least two non-na samples per group
          curr_data = curr_data[!is.na(curr_data$value),]
          group_table = table(curr_data$group)
          if(any(group_table < 2)){
            groups_to_discard = names(group_table)[group_table < 2]
            if("control" %in% groups_to_discard){
              # if the timepoint is control, discard the whole analyte
              message(sprintf("Group 'control' in %s %s %s %s has less than 2 values. Removing analyte.", 
                              PANEL, TISSUE, SEX, analyte))
              next
            }else{
              # if not, discard that timepoint
              message(sprintf("Group %s in %s %s %s %s has less than 2 values. Removing remaining samples in this group.",
                              paste0(groups_to_discard, collapse=', '), PANEL, TISSUE, SEX, analyte))
              curr_data = curr_data[!curr_data$group %in% groups_to_discard,]
            }
          }
          
          all_samples = all_meta[sex == SEX, viallabel]
          removed_samples = all_samples[!all_samples %in% as.character(curr_data$viallabel)]
          
          # training-dea
          lm_full = stats::lm(value ~ group + log2_CHEX4, data=curr_data)
          lm_reduced = stats::lm(value ~ log2_CHEX4, data=curr_data)
          lrt = lmtest::lrtest(lm_full, lm_reduced)
          
          # save results to table 
          training_sex[[SEX]] = data.frame(assay="IMMUNO",
                                           assay_code='immunoassay',
                                           dataset=PANEL,
                                           feature_ID=analyte,
                                           tissue=TISSUE,
                                           lrt=lrt$LogLik[2],
                                           p_value=lrt$`Pr(>Chisq)`[2],
                                           removed_samples=ifelse(length(removed_samples)==0, 
                                                                  NA_character_, 
                                                                  paste0(removed_samples, collapse=',')),
                                           full_model="~group+log2_CHEX4",
                                           reduced_model="~log2_CHEX4")
        }
        # meta-analysis of p-values
        if(length(training_sex) > 1){
          # merge male and female results for this feature
          male = training_sex[['male']]
          female = training_sex[['female']]
          merged = merge(male, female, by=c('feature_ID','dataset','assay','assay_code', 'tissue','reduced_model','full_model'),
                         suffixes=c("_male","_female"), all=T)
          # calculate meta-analysis p-value (one per feature)
          pvals = c(merged$p_value_male, merged$p_value_female)
          pvals = pvals[!is.na(pvals)]
          merged$p_value = sumlog(pvals)$p
          # save results to list 
          training_res[[i]] = merged
          i=i+1
        }else{
          sex_label = names(training_sex)[1]
          df = training_sex[[1]]
          colnames(df)[colnames(df) == 'lrt'] = sprintf('lrt_%s', sex_label)
          colnames(df)[colnames(df) == 'removed_samples'] = sprintf('removed_samples_%s', sex_label)
          training_res[[i]] = df
          i=i+1
        }
      }
    }
  }
  training = as.data.frame(data.table::rbindlist(training_res, fill=TRUE))
  training$removed_samples_male = as.character(training$removed_samples_male)
  training$removed_samples_female = as.character(training$removed_samples_female)
  training$tissue_code = MotrpacRatTraining6moData::TISSUE_ABBREV_TO_CODE[training$tissue]

  # reorder columns for consistency 
  training = training[,c(
    "feature_ID",
    "assay",
    "assay_code",
    "dataset",
    "tissue",
    "tissue_code",
    "removed_samples_male",
    "removed_samples_female",
    "lrt_male",
    "lrt_female",
    "p_value_male",
    "p_value_female",
    "full_model",
    "reduced_model",
    "p_value"
  )]
  
  return(training)
}
