## usethis namespace: start
#' @importFrom data.table as.data.table data.table rbindlist setnames := copy
#' @import MotrpacRatTraining6moData
#' @importFrom magrittr %>%
#' @importFrom ggplot2 guides scale_size theme annotate guide_legend aes element_rect ggplot
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats na.omit prcomp sd median poly
#' @importFrom utils download.file data globalVariables
## usethis namespace: end
NULL

## Internal assignment to avoid data.table note about `.`
`.` = list

## Initialize variables used within tidyr or data.table scope to avoid note
# See: https://github.com/Rdatatable/data.table/issues/850
utils::globalVariables(
  c("PC",
    "denominator",
    "fscore_female",
    "fscore_male",
    "full_model",
    "full_model_female",
    "full_model_male",
    "gene_id",
    "group",
    "group_tp1_0",
    "group_tp2_0",
    "group_tp4_0",
    "group_tp8_0",
    "group_tp8_1",
    "lrt",
    "numerator",
    "p_value",
    "p_value_female",
    "p_value_male",
    "pid",
    "reduced_model",
    "reduced_model_female",
    "reduced_model_male",
    "removed_samples_female",
    "removed_samples_male",
    "sex_group",
    "shrunk_logFC",
    "shrunk_logFC_se",
    "stat",
    "tr1_0",
    "tr2_0",
    "tr4_0",
    "tr8_0",
    "tr8_1",
    "size", # immuno_timewise_da
    "cytokine", # immuno_timewise_da
    "test_grouping", # immuno_timewise_da
    "terms", # immuno_timewise_da
    "coef", # immuno_timewise_da
    "coef_se", # immuno_timewise_da
    "pval", # immuno_timewise_da
    "panel" # immuno_timewise_da
  ))
