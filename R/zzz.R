## usethis namespace: start
#' @import MotrpacRatTraining6moData
#' @importFrom data.table as.data.table data.table rbindlist setnames := copy melt setkey dcast
#' @importFrom magrittr %>%
#' @importFrom ggplot2 guides scale_size theme annotate guide_legend aes element_rect ggplot geom_bar theme_classic facet_wrap scale_fill_manual coord_flip labs scale_x_discrete element_blank geom_line
#' @importFrom stats na.omit prcomp sd median poly phyper cor p.adjust
#' @importFrom utils download.file data globalVariables combn
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
    "size",
    "is_control",
    "x_cvs",
    "size", # immuno_timewise_da
    "cytokine", # immuno_timewise_da
    "test_grouping", # immuno_timewise_da
    "terms", # immuno_timewise_da
    "coef", # immuno_timewise_da
    "coef_se", # immuno_timewise_da
    "pval", # immuno_timewise_da
    "panel", # immuno_timewise_da
    "panel_name", #immuno_training_da
    ########################################### enrichment_network_vis()
    "gene_symbol", #replace_ensembl_with_symbol
    "ensembl_gene", #replace_ensembl_with_symbol
    "adj_p_value",
    "term_name",
    "term_id",
    "intersection",
    "intersection_original",
    "n_dataset",
    "n_datasets",
    "dataset",
    "datasets",
    "ome",
    "tissue",
    "tissues",
    "symbols",
    "n_genes",
    "genes",
    "computed_p_value",
    "term_size",
    "query_size",
    "enriched",
    "intersection_formatted",
    "V1",
    "V2",
    "similarity_score",
    "edges",
    "points", 
    "genes_Var1",
    "genes_Var2",
    "pathway_class",
    "pathway_subclass",
    "Var",
    "title_br", 
    "sumlog_p",
    "enriched_br",
    "intersection_br",
    "id",
    "group_number",
    "value",
    "color.background",
    "intersection_label",
    "from",
    "to",
    "color.highlight.border",
    "sumlog_p_log10",
    ########################################### cluster_pathway_enrichment()
    "cluster",
    "effective_domain_size",
    "gost_adj_p_value",
    "i",
    "intersection_size",
    "kegg_id",
    ########################################### get_peak_annotations()
    ".SD", 
    "annotation",
    "chrom",
    "dist_downstream",
    "dist_upstream",
    "end",
    "geneEnd",
    "geneStart",
    "geneStrand",
    "relationship_to_gene",
    "short_annotation",
    "start",
    "intersection_ensembl" # pathway_hypergeom_test()
  ))
