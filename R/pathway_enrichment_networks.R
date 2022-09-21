#' Pathway enrichment network
#' 
#' Plot an interactive network of pathway enrichments. 
#' 
#' @param pw_enrich_res data frame of enrichment results. 
#'   May be a subset of [MotrpacRatTraining6moData::GRAPH_PW_ENRICH] 
#'   or a data frame returned from [cluster_pathway_enrichment()].
#'   Columns must include "adj_p_value", "ome", "tissue", "intersection", "computed_p_value", 
#'   "term_name", "term_id". 
#' @param feature_to_gene data frame, map between \code{intersection_id_type} and gene symbols. 
#'   Columns must include "gene_symbol" and "ensembl_gene" if \code{intersection_id_type == "ensembl_gene"}.
#'   [MotrpacRatTraining6moData::FEATURE_TO_GENE] by default. 
#' @param intersection_id_type character, type of gene identifier used to define 
#'   \code{pw_enrich_res$intersection}, either "ensembl_gene" or "gene_symbol"
#' @param similarity_cutoff numeric, edges are drawn between pairs of pathways if they have 
#'   a similarity metric above this value. 0.375 by default. 
#' @param adj_pval_cutoff numeric, pathway enrichments are only considered for 
#'   the network if the corresponding adjusted p-value (\code{pw_enrich_res$adj_p_value}) 
#'   is less than this value. 0.1 (10% FDR) by default. 
#' @param title character, plot title. \code{NULL} by default. 
#' @param label_nodes boolean, whether to label groups at all. \code{TRUE} by default.
#'   If \code{label_nodes = FALSE}, \code{add_group_label_nodes} is ignored. 
#' @param add_group_label_nodes boolean, whether to label groups with nodes. 
#'   If \code{FALSE}, use a standard legend instead. 
#' @param parent_pathways named list, map of KEGG and REAC pathway term ID to parent pathway.
#'   Used to create labels for clusters of pathway enrichments. 
#'   List names must correspond to values in \code{pw_enrich_res$term_id}. 
#'   [MotrpacRatTraining6moData::PATHWAY_PARENTS] by default. 
#'   If NULL, then only the pathway names are used to create cluster labels, which
#'   is quite meaningless is all pathway names are unique. In that case,
#'   it is recommended to set \code{label_nodes} to \code{FALSE}.
#' @param multitissue_pathways_only boolean. If \code{TRUE}, only include pathways 
#'   in the network if they are significantly enriched in more than one tissue
#'   in \code{pw_enrich_res}.
#' @param include_metab_singletons boolean. If \code{TRUE}, include pathways enriched only 
#'   by metabolites (ome METAB) as singleton nodes
#' @param similarity_scores_file character or \code{NULL}, path in which to save 
#'   pairwise pathway similarity scores as an RDS file. If this file exists, the pairwise 
#'   similarity scores are loaded from the file instead of being recalculated. 
#'   \code{NULL} by default. 
#' @param out_html character, output file for HTML. "/dev/null" by default. 
#' @param overwrite_html boolean, whether to overwrite \code{out_html} if it already exists.
#'   \code{TRUE} by default. If \code{out_html} exists and \code{overwrite_html=FALSE},
#'   then the path \code{out_html} is returned. 
#' @param return_html boolean, whether to return the path to \code{out_html}. If
#'   \code{FALSE}, return the \code{visNetwork} object instead. \code{FALSE} by default. 
#' @param return_graph_for_cytoscape boolean. If \code{TRUE}, return an \code{igraph}
#'   graph object that can be exported to Cytoscape.
#' @param verbose boolean, whether to print verbose output
#' 
#' @return 
#' If \code{return_graph_for_cytoscape=FALSE} and \code{return_html=FALSE}, return a 
#' \code{visNetwork} graph. If \code{return_graph_for_cytoscape=TRUE}, return an 
#' igraph object. If \code{return_html=TRUE}, return the path to the HTML file in
#' which the interactive \code{visNetwork} was saved. 
#' 
#' @export
#' 
#' @importFrom igraph graph_from_data_frame cluster_louvain
#' @importFrom visNetwork visNetwork visGroups visInteraction visIgraphLayout visLegend visSave
#' @importFrom RColorBrewer brewer.pal
#' 
#' @seealso Ancillary functions include [calc_similarity_metric()], [replace_ensembl_with_symbol()],
#'   [format_gene_symbols()], [edge_intersection()], [collapse_p()], [add_line_breaks()], 
#'   [format_gene_lists()], and [cleanup()]. 
#'   The input \code{pw_enrich_res} can also be generated with [cluster_pathway_enrichment()]. 
#' 
#' @examples 
#' # Example 1: Plot an interactive network of pathway enrichments from the HEART:8w_F1_M1 node,
#' # i.e., features that are up-regulated in both males and females after 8 weeks of training
#' enrich_res = MotrpacRatTraining6moData::GRAPH_PW_ENRICH
#' enrich_res = enrich_res[enrich_res$cluster == "HEART:8w_F1_M1",]
#' enrichment_network_vis(enrich_res, add_group_label_nodes = TRUE) 
#' 
#' \dontrun{
#' # Example 2: Export the above network to Cytoscape.
#' # Cytoscape must be running locally for this to work. 
#' library(RCy3)
#' g = enrichment_network_vis(enrich_res, return_graph_for_cytoscape = T) 
#' RCy3::cytoscapePing()
#' RCy3::createNetworkFromIgraph(g, new.title='HEART:8w_F1_M1')
#' }
#' 
#' # Example 3: Plot an interactive network of pathway enrichments corresponding to 
#' # features that are up-regulated in both sexes at 8 weeks in any of the 3 muscle tissues. 
#' # Only include pathways that are enriched in at least 2 muscle tissues.
#' # Plot the color legend on the left instead of using colored labels pointing to nodes 
#' enrich_res = MotrpacRatTraining6moData::GRAPH_PW_ENRICH
#' enrich_res = enrich_res[enrich_res$tissue %in% c("SKM-GN","SKM-VL","HEART") &
#'                          grepl(":8w_F1_M1",enrich_res$cluster),]
#' enrichment_network_vis(enrich_res,
#'                        similarity_cutoff = 0.3,
#'                        title = "Muscle 8w_F1_M1",
#'                        add_group_label_nodes = FALSE,
#'                        multitissue_pathways_only = TRUE) 
#' 
enrichment_network_vis = function(pw_enrich_res, 
                                  feature_to_gene = MotrpacRatTraining6moData::FEATURE_TO_GENE,
                                  intersection_id_type = 'ensembl_gene',
                                  similarity_cutoff = 0.375,
                                  adj_pval_cutoff = 0.1,
                                  title = NULL,
                                  label_nodes = TRUE, 
                                  add_group_label_nodes = FALSE,
                                  parent_pathways = MotrpacRatTraining6moData::PATHWAY_PARENTS,
                                  multitissue_pathways_only = FALSE,
                                  include_metab_singletons = TRUE,
                                  similarity_scores_file = NULL,
                                  out_html = "/dev/null",
                                  overwrite_html = TRUE,
                                  return_html = FALSE,
                                  return_graph_for_cytoscape = FALSE,
                                  verbose = TRUE){
  
  ptm = proc.time()
  set.seed(123)
  
  if(!overwrite_html & file.exists(out_html) & out_html != "/dev/null" & return_html){
    message(sprintf("HTML file %s already exists. Returning path", out_html))
    return(out_html)
  }
  
  if(return_graph_for_cytoscape & return_html){
    stop("Only one of 'return_html' and 'return_graph_for_cytoscape' can be set to TRUE.")
  }
  
  if(verbose) message("Formatting inputs...")
  
  # check format of enrich_res 
  sub_enrich = data.table::data.table(pw_enrich_res)
  req_cols = c("adj_p_value", "ome", "tissue", "intersection", "computed_p_value", "term_name", "term_id")
  if(!all(req_cols %in% colnames(sub_enrich))){
    stop(sprintf("The input is missing at least one of the required columns:\n $s", paste(req_cols, collapse=", ")))
  }
  
  # check format of feature_to_gene
  feature_to_gene = data.table::data.table(feature_to_gene)
  if(!"gene_symbol" %in% colnames(feature_to_gene)){
    stop("The feature-to-gene map is missing the required column 'gene_symbol'.")
  }
  if(intersection_id_type == 'ensembl_gene' & !"ensembl_gene" %in% colnames(feature_to_gene)){
    stop("The feature-to-gene map is missing the column 'ensembl_gene'. This column is required when `intersection_id_type == 'ensembl_gene'`.")
  }
  
  curr_tissues = unique(sub_enrich[,tissue])
  if(length(curr_tissues) < 2 & multitissue_pathways_only){
    message("'multitissue_pathways_only' is set to TRUE but only 1 tissue is included in the input. Setting 'multitissue_pathways_only' to FALSE.")
    multitissue_pathways_only = F
  }
  
  # only consider significant results 
  clust1sig = sub_enrich[adj_p_value < adj_pval_cutoff]
  
  # skip if there are no pathways
  if(length(unique(clust1sig[,term_id])) == 0){
    message("No significant enrichments.")
    return()
  } 
  
  # skip if there is only one significant pathway
  if(length(unique(clust1sig[,term_id])) < 2){
    message("Fewer than 2 significant enrichments. Not suitable for a network view.")
    return()
  }
  
  if(verbose) message("Subsetting feature-to-gene map...")
  # filter map to make string replacement faster - BIG time saver 
  all_genes = unique(unname(unlist(sapply(clust1sig[,intersection], function(x) unname(unlist(strsplit(x, ',')))))))
  feature_to_gene = feature_to_gene[get(intersection_id_type) %in% all_genes]
  
  if(verbose) message("Calculating similarity metric between pairs of enriched pathways...")
  ## calculate similarity metric between all pairs of results 
  # merge features from different omes and tissues 
  clust1sig[,n_dataset := 1]
  clust1sig[,dataset := sprintf('%s:%s',ome,tissue)]
  
  if(intersection_id_type=="ensembl_gene"){
    # convert ensembl to gene symbol now 
    map = unique(feature_to_gene[,.(ensembl_gene, gene_symbol)])
    map = map[!is.na(ensembl_gene)]
    data.table::setkey(map, ensembl_gene)
    clust1sig[,symbols := replace_ensembl_with_symbol(intersection, map), by = 1:nrow(clust1sig)]
  }else if(intersection_id_type=="gene_symbol"){
    clust1sig[,symbols := format_gene_symbols(intersection), by = 1:nrow(clust1sig)] 
  }else{
    stop("'intersection_id_type' must be one of c('ensembl_gene','gene_symbol').")
  }
  
  clust1sig[,n_genes := as.numeric(gsub(":.*","",symbols))]
  clust1sig[,genes := gsub(".*:","",symbols)]
  clust1sig[ome == "METAB", genes := "unknown (METAB)"]
  
  # keep omes separate
  clust1sig[,enriched := sprintf('<b>%s:</b> %s',dataset,genes)]
  
  clust1sig_collapsed = clust1sig[,list(intersection_formatted = paste0(genes, collapse=', '),
                                        genes = paste0(genes, collapse=', '),
                                        intersection_original = paste0(intersection, collapse=','), 
                                        sumlog_p = collapse_p(computed_p_value),
                                        #term_size = sum(term_size),
                                        #query_size = sum(query_size),
                                        omes = paste0(unique(ome), collapse=', '),
                                        tissues = paste0(unique(tissue), collapse=', '),
                                        datasets = paste0(dataset, collapse='; '),
                                        n_datasets = sum(n_dataset),
                                        enriched = paste0(enriched, collapse="<br>")), 
                                  by=.(term_name, term_id)]
  
  stopifnot(!any(duplicated(clust1sig_collapsed[,term_id])))
  
  if(!include_metab_singletons){
    # remove pathways (nodes) if they have no intersection (i.e. METAB-only enrichments)
    clust1sig_collapsed = clust1sig_collapsed[intersection_formatted != "NA"]
    clust1sig_collapsed = clust1sig_collapsed[!is.na(intersection_formatted)]
    clust1sig_collapsed = clust1sig_collapsed[intersection_formatted != ""]
    clust1sig_collapsed = clust1sig_collapsed[genes != "NA"]
    clust1sig_collapsed = clust1sig_collapsed[!is.na(genes)]
    clust1sig_collapsed = clust1sig_collapsed[genes != ""]
    clust1sig_collapsed = clust1sig_collapsed[genes != "unknown (METAB)"]
    
    # skip if there is only one significant pathway
    if(length(unique(clust1sig_collapsed[,term_id])) < 2){
      message("Fewer than 2 significant enrichments after removing metabolomics enrichments.")
      return()
    }
  }
  
  # if multitissue_pathways_only == T, remove enrichments driven by a single tissue 
  if(multitissue_pathways_only){
    clust1sig_collapsed = clust1sig_collapsed[grepl(",", tissues)]
    # skip if there is only one significant pathway
    if(length(unique(clust1sig_collapsed[,term_id])) < 2){
      message("Fewer than 2 significant enrichments after removing single-tissue enrichments.")
      return()
    }
  }
  
  if(length(unique(clust1sig_collapsed[,term_id])) < 2){
    message("Fewer than 2 significant enrichments after filtering.")
    return()
  }
  
  # calculate similarity metric 
  # make pairs 
  allpw = unique(clust1sig_collapsed[,term_id])
  pairs = data.table::data.table(t(utils::combn(allpw, 2, simplify = T)))
  pairs = pairs[V1 != V2]
  pairs[,similarity_score := NA_real_]
  
  generate_scores = T
  if(!is.null(similarity_scores_file)){
    if(!endsWith(tolower(similarity_scores_file), "rds")){
      message("'similarity_scores_file' does not end with '.rds'. Appending '.RDS' suffix.")
      similarity_scores_file = paste0(similarity_scores_file, ".RDS")
      message(sprintf("Similarity scores will be saved in and/or read from '%s'.", similarity_scores_file))
    }
    if(file.exists(similarity_scores_file)){
      pairs = data.table::as.data.table(readRDS(file = similarity_scores_file))
      generate_scores = F
    }
  }
  
  if(generate_scores){
    # iterate through the rows
    for(i in 1:nrow(pairs)){
      v1_pw = pairs[i, V1]
      v2_pw = pairs[i, V2]
      # get string1
      v1_members = clust1sig_collapsed[term_id == v1_pw, intersection_original]
      # get string2
      v2_members = clust1sig_collapsed[term_id == v2_pw, intersection_original]
      if(v1_members == "NA" | v2_members == "NA"){
        if(include_metab_singletons){
          if(v1_members == "NA" & v2_members == "NA"){
            # both are METAB PWs
            # draw an edge if at least one word overlaps between name and parents
            if(!is.null(parent_pathways)){
              v1_parent = gsub(".*; ","",parent_pathways[[v1_pw]])
              v2_parent = gsub(".*; ","",parent_pathways[[v2_pw]])
            }else{
              v1_parent = ""
              v2_parent = ""
            }
            v1_words = cleanup(paste(clust1sig_collapsed[term_id == v1_pw, term_name], v1_parent))
            v2_words = cleanup(paste(clust1sig_collapsed[term_id == v2_pw, term_name], v2_parent))
            if(length(intersect(v1_words, v2_words))>0){
              s = similarity_cutoff
            }else{
              s = 0
            }
          }else{
            s = 0
          }
        }else{
          s = 0
        }
      }else{
        s = calc_similarity_metric(v1_members, v2_members)
      }
      pairs[i,similarity_score := s]
    }
    
    if(!is.null(similarity_scores_file)){
      saveRDS(pairs, file = similarity_scores_file)
    }
  }
  
  if(verbose) message("Constructing graph...")
  points = clust1sig_collapsed
  data.table::setnames(points, 'term_id', 'Var')
  # points[,symbols := unname(unlist(sapply(intersection, replace_ensembl_with_symbol, feature_to_gene)))]
  # points[,n_genes := as.numeric(gsub(":.*","",symbols))]
  # points[,genes := gsub(".*:","",symbols)]
  
  #if(verbose) message("Defining edges...")
  edges = pairs[similarity_score >= similarity_cutoff]
  
  # skip if there are no edges
  if(nrow(edges) == 0){
    message(sprintf("No similarity scores > %s.", similarity_cutoff))
    return()
  } 
  
  # merge back with points
  edges_v1 = merge(edges, points, by.x='V1', by.y='Var', all.x=T)
  edges = merge(edges_v1, points, by.x='V2', by.y='Var', all.x=T, suffixes = c("_Var1", "_Var2"))
  
  # add intersection (gene symbols)
  edges[,intersection := edge_intersection(genes_Var1, genes_Var2), by=1:nrow(edges)]
  
  # remove nodes with only one supporting gene 
  if(include_metab_singletons){
    points = points[grepl(",", intersection_formatted) | intersection_formatted == "unknown (METAB)"]
  }else{
    points = points[grepl(",", intersection_formatted)]
  }
  
  if(nrow(points) == 0){
    message("No significant pathway enrichments driven by more than one gene.")
    return()
  }
  
  #if(verbose) message("Final touches...")
  
  #if(verbose) message("Labelling groups by most frequent pathway subclass...")
  
  # add class and subclass to edges 
  points[,pathway_class := parent_pathways[Var]]
  points[,pathway_subclass := gsub(".*; ","",pathway_class)]
  
  # if gene lists are longer than 50 characters, add line breaks 
  # make sure to skip METAB nodes 
  points[,enriched_br := format_gene_lists(enriched), by = 1:nrow(points)]
  
  # if pathway name + parent is longer than 50 characters, add line breaks 
  points[,title := sprintf("<b>%s</b> (%s)", term_name, pathway_subclass)]
  points[,title_br := add_line_breaks(title, sep=' '), by = 1:nrow(points)]
  
  # we will color by group later
  viznetwork_nodes = data.table::data.table(id = points[,Var],
                                            value = 100*(points[,n_datasets]^2),
                                            title = sprintf("%s<br><b>P-val:</b> %s<br>%s", # tooltip
                                                            points[,title_br],
                                                            signif(points[,sumlog_p], 2),
                                                            points[,enriched_br]),
                                            color.highlight.background = 'yellow',
                                            color.border = 'black',
                                            shadow = F,
                                            physics = F,
                                            borderWidth = 1,
                                            borderWidthSelected = 3,
                                            pathway_subclass = points[,pathway_subclass])
  
  # add line breaks
  # allow for METAB
  edges[,intersection_label := sprintf("<b>Intersection:</b> %s",intersection)]
  edges[,intersection_br := add_line_breaks(intersection_label), by = 1:nrow(edges)]
  
  viznetwork_edges = data.table::data.table(from = edges[,V1],
                                            to = edges[,V2],
                                            value = (edges[,similarity_score]*50),
                                            title = sprintf("<b>Similarity score:</b> %s<br>%s", # tooltip
                                                            round(edges[,similarity_score],2),
                                                            edges[,intersection_br]),
                                            color.color = "#848484", # visNetwork default
                                            physics = F,
                                            color.highlight = "#545454",
                                            labelHighlightBold = T)
  
  # remove edges without any nodes
  viznetwork_edges = viznetwork_edges[from %in% viznetwork_nodes[,id] & to %in% viznetwork_nodes[,id]]
  
  # define group by igraph cluster 
  net = igraph::graph_from_data_frame(d=viznetwork_edges, vertices=viznetwork_nodes, directed=F)
  clust_to_pw = igraph::cluster_louvain(net)
  viznetwork_nodes[,group := as.character(clust_to_pw$membership)]
  
  #if(verbose) message("Clean up nodes...")
  
  # make cluster singletons all their own cluster 0
  clust_sizes = as.list(table(clust_to_pw$membership))
  singletons = names(clust_sizes)[clust_sizes==1]
  viznetwork_nodes[group %in% singletons, group := '0']
  
  # use most common subclass for group name
  # use these in the legend 
  viznetwork_nodes[,group_number := group]
  if(label_nodes){
    group_labels = sapply(unique(viznetwork_nodes[,group_number]), function(x){
      # get corresponding subclasses
      subclass = viznetwork_nodes[group_number == x, pathway_subclass]
      if(all(subclass=="NULL")){
        return(sprintf("%s: NULL", x))
      }
      subclass = subclass[!is.na(subclass)]
      subclass = subclass[!is.null(subclass)]
      subclass = subclass[subclass!="NULL"]
      subclass_counts = unlist(table(subclass))
      subclass_counts = subclass_counts[order(subclass_counts, decreasing=T)]
      maxclass = names(subclass_counts)[1]
      return(sprintf("%s: %s (%s/%s)", x, maxclass, subclass_counts[[maxclass]], sum(subclass_counts)))
    })
    # add to node data.table
    viznetwork_nodes[,group := group_labels[group_number]]
  }
  
  # help by initializing 
  vn0 = visNetwork::visNetwork(viznetwork_nodes, viznetwork_edges) %>%
    visNetwork::visGroups(groupname = viznetwork_nodes[1,group], color = "red", shape = "triangle") 
  
  # define colors 
  # use Set2, then Set3, then Set1 + Set 3
  n_clusters = length(unique(viznetwork_nodes[,group]))
  if(n_clusters > 12 & n_clusters < 21){
    hex = c(RColorBrewer::brewer.pal(12, "Set3"), 
            RColorBrewer::brewer.pal(9, "Set1"))
  }else if(n_clusters > 8){
    hex = RColorBrewer::brewer.pal(12, "Set3")
  }else if(n_clusters <= 8){
    hex = RColorBrewer::brewer.pal(8, "Set2")
  }else{
    warning("More than 21 clusters? I didn't think that would happen.")
    hex = c(RColorBrewer::brewer.pal(12, "Set3"), 
            RColorBrewer::brewer.pal(9, "Set1"), 
            RColorBrewer::brewer.pal(7, "Set2"))
  }
  hex = hex[1:n_clusters]
  names(hex) = unique(viznetwork_nodes[,group])
  
  #if(verbose) message("Formatting legend...")
  
  # set color by group 
  viznetwork_nodes[,color.background := hex[group]]
  viznetwork_nodes[,color.highlight.border := hex[group]]
  
  if(return_graph_for_cytoscape){
    
    ## add a few more things to viznetwork_nodes
    # clean pathway name
    # sumlog p-values 
    vn = merge(viznetwork_nodes, points[,.(term_name, Var, sumlog_p, omes, tissues, n_datasets)], by.x='id', by.y='Var')
    vn[,title := NULL]
    
    ## add similarity score back to edges
    ve = merge(viznetwork_edges, edges[,.(V2, V1, similarity_score)], by.x=c('to','from'), by.y=c('V2','V1'))
    
    # binary column for each ome and tissue
    omes = unique(unname(unlist(strsplit(vn[,omes],", "))))
    tissues = unique(unname(unlist(strsplit(vn[,tissues],", "))))
    for(c in c(omes, tissues)){
      vn[,(c) := 0]
    }
    for(i in 1:nrow(vn)){
      curr_omes = unname(unlist(strsplit(vn[i, omes], ", ")))
      curr_tissues = unname(unlist(strsplit(vn[i, tissues], ", ")))
      for(c in c(curr_omes, curr_tissues)){
        vn[i,(c) := 1]
      }
    }
    
    # add columns for other omes and tissues 
    all_tissues = MotrpacRatTraining6moData::TISSUE_ABBREV
    all_tissues = all_tissues[!all_tissues %in% c("VENACV","OVARY","TESTES")]
    all_omes = MotrpacRatTraining6moData::ASSAY_ABBREV
    
    for(col in c(all_tissues, all_omes)){
      if(! col %in% colnames(vn)){
        vn[,(col) := 0]
      }
    }
    
    # add -log10 sumlog p for node size
    vn[,sumlog_p_log10 := round(-log10(sumlog_p),2)]
    
    # clean up nodes and edges
    data.table::setnames(vn, 
                         c('group_number', 'group'), 
                         c('group', 'group_label'))
    
    ve[,title := NULL]
    
    cytoscape_igraph = igraph::graph_from_data_frame(ve, directed = F, vertices = vn)
    return(cytoscape_igraph)
  }
  
  # add group label nodes
  if(label_nodes){
    if(add_group_label_nodes){
      # create additional nodes with a strong weight to a random member of that group 
      group_nodes = data.table::data.table(id = names(hex),
                                           value = max(viznetwork_nodes[,value])*2, 
                                           shape = "box",
                                           label = names(hex),
                                           title = names(hex), 
                                           group = names(hex),
                                           physics = F,
                                           labelHighlightBold = F, 
                                           color.background = unname(hex),
                                           color.highlight.background = unname(hex),
                                           color.highlight.border = unname(hex),
                                           color.border = unname(hex),
                                           borderWidthSelected = 0, 
                                           borderWidth = 1, 
                                           shadow = T)
      group_nodes = group_nodes[!grepl("^0: ", id)] # remove label for singleton clusters 
      
      viznetwork_nodes = data.table::rbindlist(list(viznetwork_nodes, group_nodes), fill=T)
      
      # add edges
      # for each group, pick node (smallest pval) where the category is also the label
      rowselect = viznetwork_nodes[,grepl(pathway_subclass, group, fixed = TRUE), by=1:nrow(viznetwork_nodes)][,V1]
      
      sub = viznetwork_nodes[rowselect]
      sub = sub[order(value, decreasing=T)]
      sub = sub[!duplicated(group)]
      
      group_edges = data.table(from = sub[,id],
                               to = sub[,group],
                               labelHighlightBold = F, 
                               color = sub[,color.background],
                               value = max(viznetwork_edges[,value]),
                               labelHighlightBold = F,  
                               title = NULL,
                               physics = F,
                               shadow = T)
      group_edges = group_edges[!grepl("^0: ", to)]
      viznetwork_edges = data.table::rbindlist(list(viznetwork_edges, group_edges), fill=T)
      
      v1 = visNetwork::visNetwork(viznetwork_nodes, viznetwork_edges, main=title) %>%
        visNetwork::visInteraction(tooltipDelay = 10,
                                   tooltipStay = Inf,
                                   hideEdgesOnZoom = T) %>%
        visNetwork::visIgraphLayout(layout = "layout_nicely") 
    }else{
      # nodes for legend
      # order by cluster number
      nums = as.numeric(gsub(":.*","",names(hex)))
      names(nums) = names(hex)
      nums = nums[order(nums, decreasing=F)]
      hex = hex[names(nums)]
      lnodes = data.frame(label = names(hex),
                          shape = c("box"), 
                          color = unname(hex),
                          physics = F)
      
      v1 = visNetwork::visNetwork(viznetwork_nodes, viznetwork_edges, main=title) %>%
        visNetwork::visLegend(zoom = F,
                              addNodes = lnodes,
                              useGroups = F,
                              width = 0.3,
                              stepY = 40) %>%
        visNetwork::visInteraction(tooltipDelay = 10,
                                   tooltipStay = Inf,
                                   hideEdgesOnZoom = T) %>%
        visNetwork::visIgraphLayout(layout = "layout_nicely") 
    }
  }else{
    v1 = visNetwork::visNetwork(viznetwork_nodes, viznetwork_edges, main=title) %>%
      visNetwork::visInteraction(tooltipDelay = 10,
                                 tooltipStay = Inf,
                                 hideEdgesOnZoom = T) %>%
      visNetwork::visIgraphLayout(layout = "layout_nicely") 
  }
    
  if(out_html != "/dev/null"){ 
    if( (file.exists(out_html) & overwrite_html) | !file.exists(out_html)){
      visNetwork::visSave(v1, file = out_html)
      message(sprintf("Network saved in %s", out_html))
      # adjust browser window size 
      cmd = sprintf('sed -i \'s|"browser":{"width":960,"height":500,"padding":40,"fill":false}|"browser":{"width":960,"height":960,"padding":0,"fill":false}|\' %s', out_html)
      system(cmd)
    }
  }
  
  if(verbose){
    message(sprintf("Elapsed time: %ss", round((proc.time() - ptm)[3], 3)))
    cat("\n")
  }
  
  if(return_html){
    return(out_html)
  }else{
    return(v1)
  }
}


#' Calculate pathway similarity metric 
#' 
#' Calculate EnrichmentMap's similarity metric for pathways, which is 50% Jaccard
#' index and 50% overlap score.
#' Function is used internally in [enrichment_network_vis()]. 
#' 
#' @param string1 character, comma-separated list of intersection with pathway1
#' @param string2 character, comma-separated list of intersection with pathway2
#' 
#' @return numeric similarity score
#' 
#' @seealso [enrichment_network_vis()]
#' 
#' @keywords internal
#' 
calc_similarity_metric = function(string1, string2){
  # use Cytoscape Enrichment Map similarity score:
  # jaccard = [size of (A intersect B)] / [size of (A union B)]
  # overlap = [size of (A intersect B)] / [size of (minimum( A , B))]
  # similarity_score = 0.5*jaccard + 0.5*overlap
  set1 = unname(unlist(strsplit(string1, ',')))
  set2 = unname(unlist(strsplit(string2, ',')))
  
  jaccard = length(intersect(set1, set2)) / length(union(set1, set2))
  overlap = length(intersect(set1, set2)) / min(length(set1), length(set2))
  similarity_score = 0.5*jaccard + 0.5*overlap
  
  return(similarity_score)
}


#' Replace Ensembl ID with gene symbol
#' 
#' Function is used internally in [enrichment_network_vis()].
#' 
#' @param x character, comma-separated Ensembl IDs
#' @param map data table, mapping between Ensembl IDs and gene symbols.
#'   Must include the columns "ensembl_gene" and "gene_symbol". 
#' @param return_N boolean, whether to prepend the concatenated 
#'     gene symbols with the number of unique gene symbols.
#'     \code{TRUE} by default. 
#' @param collapse boolean, whether to collapse the gene symbols into a 
#'     comma-separated string. If not, return a vector of unique gene symbols.
#'     \code{TRUE} by default. 
#'     
#' @return either a string of comma-separated gene symbols, 
#'   a vector of unique gene symbols,
#'   or a named list with two values ("N" and "genes")
#'   depending on the values of \code{return_N} and \code{collapse}. 
#'   
#' @seealso [enrichment_network_vis()]
#' 
#' @keywords internal
#' 
replace_ensembl_with_symbol = function(x, 
                                       map, 
                                       return_N = TRUE, 
                                       collapse = TRUE){
  if(!data.table::is.data.table(map)){
    map = data.table::data.table(map)
    data.table::setkey(map, ensembl_gene)
  }

  ensembls = unlist(unname(strsplit(x, ',')))
  if(all(is.na(ensembls))){
    if(return_N & collapse){
      return("0:")
    }else if(!return_N & collapse){
      return("")
    }else if(!return_N & !collapse){
      return("")
    }else{
      return(list(N=0, genes=""))
    }
  } 
  symbols = unique(map[ensembls, on = "ensembl_gene", gene_symbol])
  #symbols = unique(map[ensembl_gene %in% ensembls, gene_symbol])
  symbols = symbols[!is.na(symbols)]
  n_genes = length(symbols)
  # remove LOC## cuz who on earth knows what that is 
  symbols = symbols[!grepl("LOC[0-9]+$",symbols)]  
  # same with RGD# 
  symbols = symbols[!grepl("RGD[0-9]+$",symbols)]
  symbols = symbols[order(symbols)]
  if(return_N & collapse){
    return(paste0(n_genes, ":", paste(symbols, collapse=', ')))
  }else if(!return_N & collapse){
    return(paste(symbols, collapse=', '))
  }else if(!return_N & !collapse){
    return(symbols)
  }else{
    return(list(N=n_genes, genes=symbols))
  }
}


#' Format gene symbols
#' 
#' Format a string of comma-separated gene symbols for use in [enrichment_network_vis()].
#'
#' @param x character, comma-separated gene symbols
#' @param return_N boolean, whether to prepend the concatenated 
#'     gene symbols with the number of unique gene symbols.
#'     \code{TRUE} by default. 
#' @param collapse boolean, whether to collapse the gene symbols into a 
#'     comma-separated string. If not, return a vector of unique gene symbols.
#'     \code{TRUE} by default. 
#'     
#' @return either a string of comma-separated gene symbols, 
#'   a vector of unique gene symbols,
#'   or a named list with two values ("N" and "genes")
#'   depending on the values of \code{return_N} and \code{collapse}. 
#'   
#' @seealso [enrichment_network_vis()]
#' 
#' @keywords internal
#'
format_gene_symbols = function(x, return_N = TRUE, collapse = TRUE){
  symbols = unlist(unname(strsplit(x, ',')))
  if(all(is.na(symbols))){
    if(return_N & collapse){
      return("0:")
    }else if(!return_N & collapse){
      return("")
    }else if(!return_N & !collapse){
      return("")
    }else{
      return(list(N=0, genes=""))
    }
  } 
  symbols = symbols[!is.na(symbols)]
  n_genes = length(symbols)
  # remove LOC## cuz who on earth knows what that is 
  symbols = symbols[!grepl("LOC[0-9]+$",symbols)]  
  # same with RGD# 
  symbols = symbols[!grepl("RGD[0-9]+$",symbols)]
  symbols = symbols[order(symbols)]
  if(return_N & collapse){
    return(paste0(n_genes, ":", paste(symbols, collapse=', ')))
  }else if(!return_N & collapse){
    return(paste(symbols, collapse=', '))
  }else if(!return_N & !collapse){
    return(symbols)
  }else{
    return(list(N=n_genes, genes=symbols))
  }
}


#' Return intersection of gene symbols
#' 
#' Function used internally in [enrichment_network_vis()].
#' 
#' @param symbol1 character, comma-separated gene symbols of intersection from start node
#' @param symbol2 character, comma-separated gene symbols of intersection from end node
#' 
#' @return string of comma-separated gene symbols in intersection 
#' 
#' @seealso [enrichment_network_vis()]
#' 
#' @keywords internal
#' 
edge_intersection = function(symbol1, symbol2){
  symbol1 = unique(unname(unlist(strsplit(symbol1, ', '))))
  symbol2 = unique(unname(unlist(strsplit(symbol2, ', '))))
  gene_intersection = intersect(symbol1, symbol2)
  return(paste0(gene_intersection, collapse=', '))
}


#' Collapse p-values 
#' 
#' Collapse p-values using [metap::sumlog()] when multiple are present.
#' Function is used internally in [enrichment_network_vis()].
#' 
#' @param ps vector of nominal p-values 
#' 
#' @return merged p-value 
#' 
#' @importFrom metap sumlog
#' 
#' @seealso [enrichment_network_vis()]
#' 
#' @keywords internal
#' 
collapse_p = function(ps){
  if(length(ps) == 1){
    return(ps)
  }else{
    return(metap::sumlog(ps)$p)
  }
}


#' Add line breaks to a string 
#' 
#' Given a long string, add line breaks "\<br\>" at specified intervals. 
#' Breaks are always added at the preceding instance of the separator. 
#' Function is used internally in [enrichment_network_vis()]. 
#' 
#' @param x character string
#' @param max_char integer, add line breaks at least every \code{max_char} characters 
#' @param sep character after which to add line break
#' 
#' @return string with additional \<br\> if necessary 
#' 
#' @seealso [enrichment_network_vis()]
#' 
#' @keywords internal
#' 
add_line_breaks = function(x, max_char=50, sep=','){
  if(nchar(x) > max_char){
    charsplits = unname(unlist(strsplit(x, "")))
    newbreaks = seq(max_char, nchar(x), by=max_char)
    for(b in newbreaks){
      # find closest preceding comma
      comma_index = max(which(charsplits[1:b] == sep))
      # replace comma 
      charsplits[comma_index] = sprintf("%s<br>", sep)
    }
    x = paste0(charsplits, collapse="")
  }
  return(x)
}


#' Format gene lists
#' 
#' Apply [add_line_breaks()] to the \code{points$enriched} column. 
#' Function is used internally in [enrichment_network_vis()]. 
#' 
#' @param x value in \code{points$enriched} column
#' 
#' @return character, new value with additional line breaks 
#' 
#' @seealso [enrichment_network_vis()], [add_line_breaks()]
#' 
#' @keywords internal
#' 
format_gene_lists = function(x){
  if(grepl("(METAB)", x)) return(x)
  splits = unname(unlist(strsplit(x, "<br>")))
  newsplits = unname(unlist(sapply(splits, add_line_breaks)))
  return(paste0(newsplits, collapse="<br>"))
}


#' Clean up pathway names
#' 
#' Remove punctuation, numbers, and common words to determine qualitative 
#' overlap between pathways enriched by metabolites versus other pathways.
#' Function used internally in [enrichment_network_vis()].
#'
#' @param x character, pathway names 
#'
#' @return vector of candidate words to determine overlap 
#' 
#' @seealso [enrichment_network_vis()]
#' 
#' @keywords internal
#' 
cleanup = function(x){
  x1 = gsub(":|-|;|'|[0-9]|,","",x)
  x2 = unique(tolower(unlist(strsplit(x1, " "))))
  x2 = x2[x2!=""]
  x2 = x2[!x2 %in% MotrpacRatTraining6mo::STOPWORDS]
  # exclude other common words
  common_pw_words = c("regulation", "cell", "diseases", "pathway", "response", "system", "human", "disease", "systems")
  x2 = x2[!x2 %in% common_pw_words]
  return(x2)
}

