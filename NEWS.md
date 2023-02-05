# MotrpacRatTraining6mo 1.6.1 (2023-02-05)  

* Accommodate RefMet IDs as feature IDs in `plot_feature_normalized_data()` and `plot_feature_logfc()`.  
* Pass the user-supplied version of the feature identifier when `plot_feature_logfc()` is called within `plot_feature_normalized_data()`.  

# MotrpacRatTraining6mo 1.6.0 (2023-01-26)

* Add functions for GSEA and PTM-SEA: `ssGSEA2_wrapper()`, `prepare_gsea_input()`, 
`prepare_ptmsea_input()`, `find_flanks()`, `load_uniprot_human_fasta()`  
* Add `extract_top_trajectories()` function. 

# MotrpacRatTraining6mo 1.5.2 (2023-01-20)

* For `enrichment_network_vis()` with argument `return_html=TRUE`, don't require GNU sed.  

# MotrpacRatTraining6mo 1.5.1 (2023-01-19)

* In `enrichment_network_vis()`, include pathway names in determination of node labels.  
* In `enrichment_network_vis()`, use `MotrpacRatTraining6moData::FEATURE_TO_GENE_FILT` as
`feature_to_gene` default instead of the unfiltered version.  
* Add more details to `enrichment_network_vis()` docs.  

# MotrpacRatTraining6mo 1.5.0 (2023-01-18)

* Add `load_feature_annotation()`.  
* Add `counts` argument to `transcript_normalize_counts()` and `atac_normalize_counts()` 
to allow user-supplied data.  
* Add `tissues`, `assays`, and `cluster` arguments to `enrichment_network_vis()` to 
provide an alternate way for users to specify results from `MotrpacRatTraining6moData::GRAPH_PW_ENRICH`.    
* Add tests for `enrichment_network_vis()`.  
* Use `signif()` instead of `round()` to display values in plot titles. 
* Add note about reproducibility issues to documentation for `transcript_normalize_counts()`.  
* Move `MotrpacRatTraining6moData` from `Imports` to `Depends`. This means `MotrpacRatTraining6moData` is also attached
(not just loaded) when `MotrpacRatTraining6mo` is attached.  
* Replace `fetch_object(name_as_string)` with `.get`, which internally uses `get(name_as_string, envir=as.environment("package:MotrpacRatTraining6moData"))`.  
* Change URLs for `load_methyl_feature_annotation()` and `load_atac_feature_annotation()`.  
* Fix bug in `plot_feature_logfc()` that prevented epigenetic features from being plotted.   
* Speed up `plot_feature_logfc()` for differential epigenetic features.  
* In `get_rdata_from_url()`, calculate log fold-change standard errors for METHYL 
differential analysis results: `logFC_se = logFC/zscore`  

# MotrpacRatTraining6mo 1.4.3 (2023-01-06)

* In `plot_feature_normalized_data()` and `plot_feature_logfc()`, handle input features that don't exist in the data.  

# MotrpacRatTraining6mo 1.4.2 (2023-01-06)

* Dodge points in `plot_feature_normalized_data()` and `plot_feature_logfc()` for readability.  
* Make `plot_feature_normalized_data()` and `plot_feature_logfc()` more accommodating of features with multiple measurements.  
* Add `add_adj_p` parameter to `plot_feature_normalized_data()`.  
* Add `return_data` parameter to `plot_feature_logfc()`. 
* Return NULL instead of error if data doesn't exist for plotting `plot_feature_normalized_data()` in loops. 

# MotrpacRatTraining6mo 1.4.1 (2023-01-05)

* Fix `custom_cluster_pathway_enrichment()` to accept a list of lists as input. 
* Use filtered feature-to-gene map to perform pathway enrichment. 
* In `cluster_pathway_enrichment()`, remove pathway enrichments driven by a single gene
for consistency with the landscape paper.   

# MotrpacRatTraining6mo 1.4.0 (2023-01-05)

* Add `gene_pathway_enrichment()` (wrapper for `gprofiler2::gost()`). 
* Add example to "Get Started" vignette. 

# MotrpacRatTraining6mo 1.3.0 (2023-01-04)

* Add `df_to_numeric()` to easily format data frames. 
* Add `plot_feature_normalized_data()` to plot normalized sample-level data for a single feature. 
* Add `plot_feature_logfc()` to plot differential analysis results for a single feature. 
* Add more examples and screenshots to the "Get Started" vignette.  

# MotrpacRatTraining6mo 1.2.0 (2022-12-27)

* Add functions to perform meta-analysis and meta-regression of redundant metabolites. 
* Add `load_training_da()` to load training summary statistics from GSC.  
* Retain order of input vial labels in result of `viallabel_to_pid()`.  
 
# MotrpacRatTraining6mo 1.1.1 (2022-11-09)

* Use `MotrpacRatTraining6moData::TRAINING_REGULATED_NORM_DATA` 
and add a `training_regulated_only` argument to speed up `plot_feature_trajectories()` 
* Fix bugs in `get_tree_plot_for_tissue()` and `get_trajectory_sizes_from_edge_sets()`
* Add `scale` argument to `call_pca_outliers()` (default `TRUE`)  
* Set default values for arguments `plot` and `verbose` in `call_pca_outliers()` (default `TRUE`)  

# MotrpacRatTraining6mo 1.1.0 (2022-10-10)

* Add `atac_call_outliers()` 
* Bug fix in `atac_prep_data()` 

# MotrpacRatTraining6mo 1.0.2 (2022-10-08)

* Use a hashmap to speed up calculation of similarity scores in `enrichment_network_viz()` (thanks, Jay Yu)
* Improve `visInteraction` settings in `enrichment_network_viz()`, i.e., show edges on zoom and slow down zoom
* Don't run memory-intensive examples and vignette code blocks to prevent GitHub Actions OOM error

# MotrpacRatTraining6mo 1.0.1 (2022-09-21)

* Bug fixes in `transcript_timewise_da()` and `transcript_prep_data()`   
* Keep examples from printing too much output   
* Rename `tutorial.Rmd` to `MotrpacRatTraining6mo`   
* Add `@keywords internal` to functions that are not exported   

# MotrpacRatTraining6mo 1.0.0 (2022-09-20)

First release
