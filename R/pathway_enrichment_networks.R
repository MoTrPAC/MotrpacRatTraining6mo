#' TODO
#' Calculate EnrichmentMap's similarity metric for pathways
#' Function used internally in `enrichment_network_vis()`
#' 
#' @param string1 comma-separated list of intersection with pathway1
#' @param string2 comma-separated list of intersection with pathway2
#' 
#' @return numeric: similarity score (50% Jaccard, 50% Overlap)
.calc_similarity_metric = function(string1, string2){
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


#' TODO
#' Replace Ensembl IDs with gene symbols
#' Function used internally in `enrichment_network_vis()`
#' 
#' @param x string: comma-separated Ensembl IDs
#' @param map data.table master_feature_to_gene
#' @param return_N boolean, whether or not to prepend the concatenated 
#'     gene symbols with the number of unique gene symbols
#' @param collapse boolean, whether or not to collapse the gene symbols into a 
#'     comma-separated string. If not, return a vector of unique gene symbols
#'     
#' @return string: comma-separated gene symbols
.replace_ensembl_with_symbol = function(x, map, return_N = T, collapse = T){
  require(data.table)
  if(!is.data.table(map)){
    map = data.table(map)
    # map = unique(feature_to_gene[,.(ensembl_gene, gene_symbol)])
    # map = map[!is.na(ensembl_gene)]
    setkey(map, ensembl_gene)
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


# TODO
.format_gene_symbols = function(x, return_N = T, collapse = T){
  require(data.table)
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

# TODO
#' Return intersection of gene symbols
#' Function used internally in `enrichment_network_vis()`
#' 
#' @param ensembl1 string: comma-separated gene symbols of intersection from start node
#' @param ensembl2 string: comma-separated gene symbols of intersection from end node
#' 
#' @return string: comma-separated gene symbols in intersection 
.edge_intersection = function(symbol1, symbol2){
  symbol1 = unique(unname(unlist(strsplit(symbol1, ', '))))
  symbol2 = unique(unname(unlist(strsplit(symbol2, ', '))))
  gene_intersection = intersect(symbol1, symbol2)
  return(paste0(gene_intersection, collapse=', '))
}


# TODO
#' Assign hex colors to continuous values 
#' Function used internally in `enrichment_network_vis()`
#' 
#' @param hex_vec e.g. `brewer.pal(9,"YlOrRd")[3:9]`
#' @param values e.g. `-log10(points[,sumlog_p])`
#' 
#' @return vector of hex codes corresponding to values, in the order they were supplied 
.get_feature_colors = function(hex_vec, values){
  require(RColorBrewer)
  my.col = colorRampPalette(hex_vec)(100)
  # get quantile from values; get corresponding quantile from my.col
  values_quant = ecdf(values)(values)
  values_percentile = ceiling(values_quant*100)
  values_colors = my.col[values_percentile]
  return(values_colors)
}


# TODO
#' Collapse p-values using `metap::sumlog()` when multiple are present
#' Function used internally in `enrichment_network_vis()`
#' 
#' @param ps vector of nominal p-values 
#' 
#' @return merged p-value 
.collapse_p = function(ps){
  require(metap)
  if(length(ps) == 1){
    return(ps)
  }else{
    return(sumlog(ps)$p)
  }
}


# TODO
#' Add line breaks to a string 
#' Function used internally in `enrichment_network_vis()`
#' 
#' @param x string
#' @param max_char add line breaks every `max_char` characters 
#' @param sep character after which to add line break
#' 
#' @return string with additional <br> if necessary 
.add_line_breaks = function(x, max_char=50, sep=','){
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


# TODO
#' Apply `.add_line_breaks()` to `points[,enriched]` column
#' #' Function used internally in `enrichment_network_vis()`
#' 
#' @param x value in `points[,enriched]` column
#' 
#' @return string, new value with additional line breaks 
.format_gene_lists = function(x){
  if(grepl("(METAB)", x)) return(x)
  splits = unname(unlist(strsplit(x, "<br>")))
  newsplits = unname(unlist(sapply(splits, .add_line_breaks)))
  return(paste0(newsplits, collapse="<br>"))
}


# TODO
#' Plot a network view of pathway enrichments for a cluster, node, or trajectory
#' 
#' @param sub_enrich data.frame of enrichment results returned from `cluster_pathway_enrichment()`, 
#'     subset to the cluster of interest. Columns must include "adj_p_value", "ome", "tissue", 
#'     "intersection", "computed_p_value", "term_size", "query_size", "intersection_size", "term_name", "term_id". 
#' @param feature_to_gene data.frame: gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/resources/motrpac-mappings-master_feature_to_gene.txt
#' @param classlist list, map of KEGG pathway term ID to class
#' @param corr_thresh similarity metric between pathways below which edges are not drawn between pathways 
#' @param adj_pval_cutoff adj_p_value threshold to define significant enrichments
#' @param verbose boolean, whether or not to print verbose output
#' @param title plot title
#' @param add_group_label_nodes boolean, whether or not to label groups with nodes. If FALSE, use a standard legend instead. 
#' @param out_html output file for HTML
#' @param overwrite_html boolean, whether or not to overwrite `out_html` if it already exists 
#' @param return_html boolean, if TRUE return the path to `out_html`
#' @param save_similarity_scores character or NULL, provide RDS path in which to save/read pairwise similarity scores. 
#'     This is a big time saver if you run more than one iteration.
#' @param multitissue_pathways_only boolean, if TRUE, only include pathways in the input enriched in more than one tissue
#' @param include_metab boolean, if TRUE include pathways enriched only by METAB as singleton nodes
#' @param return_graph_for_cytoscape boolean, if TRUE return igraph::graph_from_data_frame
#' @param intersection_id_type character, type of gene identifier used to define the intersection column in the input, one of c('ensembl_gene', 'gene_symbol')
#' 
#' @result visNetwork graph with `visIgraphLayout(v1,layout = "layout_nicely")`
enrichment_network_vis = function(sub_enrich, 
                                  feature_to_gene,
                                  classlist,
                                  corr_thresh = 0.375,
                                  adj_pval_cutoff = 0.1,
                                  verbose = TRUE,
                                  title = NULL,
                                  add_group_label_nodes = FALSE,
                                  out_html = "/dev/null",
                                  overwrite_html = TRUE,
                                  return_html = FALSE,
                                  save_similarity_scores = NULL,
                                  multitissue_pathways_only = FALSE,
                                  include_metab = TRUE,
                                  return_graph_for_cytoscape = FALSE,
                                  intersection_id_type='ensembl_gene'){
  
  
  ptm = proc.time()
  set.seed(123)
  
  if(!overwrite_html & file.exists(out_html) & out_html != "/dev/null" & return_html){
    message(sprintf("HTML file %s already exists. Returning path", out_html))
    return(out_html)
  }
  
  if(return_graph_for_cytoscape & return_html){
    stop("Only one of 'return_html' and 'return_graph_for_cytoscape' can be set to TRUE.")
  }
  
  require(igraph)
  require(visNetwork)
  require(data.table)
  require(RColorBrewer)
  require(testit)
  require(tm)
  require(MotrpacBicQC)
  
  if(verbose) message("Checking input formats...")
  
  # check format of enrich_res 
  sub_enrich = data.table(sub_enrich)
  req_cols = c("adj_p_value", "ome", "tissue", "intersection", "computed_p_value", 
               "term_size", "query_size", "intersection_size", "term_name", "term_id")
  if(!all(req_cols %in% colnames(sub_enrich))){
    stop("The input is missing required columns. Did you generate it using 'cluster_pathway_enrichment()'?")
  }
  
  # check format of feature_to_gene
  feature_to_gene = data.table(feature_to_gene)
  req_cols = c("feature_ID", "gene_symbol", "ensembl_gene", "kegg_id")
  if(!all(req_cols %in% colnames(feature_to_gene))){
    stop("The feature-to-gene map is missing required columns. Did you read it from 'gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/resources'?")
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
    setkey(map, ensembl_gene)
    clust1sig[,symbols := .replace_ensembl_with_symbol(intersection, map), by = 1:nrow(clust1sig)]
  }else if(intersection_id_type=="gene_symbol"){
    clust1sig[,symbols := .format_gene_symbols(intersection), by = 1:nrow(clust1sig)] 
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
                                        sumlog_p = .collapse_p(computed_p_value),
                                        term_size = sum(term_size),
                                        query_size = sum(query_size),
                                        omes = paste0(unique(ome), collapse=', '),
                                        tissues = paste0(unique(tissue), collapse=', '),
                                        datasets = paste0(dataset, collapse='; '),
                                        n_datasets = sum(n_dataset),
                                        enriched = paste0(enriched, collapse="<br>")), 
                                  by=.(term_name, term_id)]
  
  assert(!any(duplicated(clust1sig_collapsed[,term_id])))
  
  if(!include_metab){
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
  pairs = data.table(t(combn(allpw, 2, simplify = T)))
  pairs = pairs[V1 != V2]
  pairs[,similarity_score := NA_real_]
  
  generate_scores = T
  if(!is.null(save_similarity_scores)){
    if(!endsWith(tolower(save_similarity_scores), "rds")){
      message("'save_similarity_scores' does not end with '.rds'. Appending '.RDS' suffix.")
      save_similarity_scores = paste0(save_similarity_scores, ".RDS")
      message(sprintf("Similarity scores will be saved in and/or read from '%s'.", save_similarity_scores))
    }
    if(file.exists(save_similarity_scores)){
      pairs = as.data.table(readRDS(file = save_similarity_scores))
      generate_scores = F
    }
  }
  
  # remove punctuation and numbers
  .cleanup = function(x){
    x1 = gsub(":|-|;|'|[0-9]|,","",x)
    x2 = unique(tolower(unlist(strsplit(x1, " "))))
    x2 = x2[x2!=""]
    x2 = x2[!x2 %in% tm::stopwords()]
    # exclude other common words
    common_pw_words = c("regulation", "cell", "diseases", "pathway", "response", "system", "human", "disease", "systems")
    x2 = x2[!x2 %in% common_pw_words]
    return(x2)
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
        if(include_metab){
          if(v1_members == "NA" & v2_members == "NA"){
            # both are METAB PWs
            # draw an edge if at least one word overlaps between name and parents
            v1_parent = gsub(".*; ","",classlist[[v1_pw]])
            v2_parent = gsub(".*; ","",classlist[[v2_pw]])
            v1_words = .cleanup(paste(clust1sig_collapsed[term_id == v1_pw, term_name], v1_parent))
            v2_words = .cleanup(paste(clust1sig_collapsed[term_id == v2_pw, term_name], v2_parent))
            if(length(intersect(v1_words, v2_words))>0){
              s = corr_thresh
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
        s = .calc_similarity_metric(v1_members, v2_members)
      }
      pairs[i,similarity_score := s]
    }
    
    if(!is.null(save_similarity_scores)){
      saveRDS(pairs, file = save_similarity_scores)
    }
  }
  
  if(verbose) message("Defining nodes...")
  points = clust1sig_collapsed
  setnames(points, 'term_id', 'Var')
  # points[,symbols := unname(unlist(sapply(intersection, .replace_ensembl_with_symbol, feature_to_gene)))]
  # points[,n_genes := as.numeric(gsub(":.*","",symbols))]
  # points[,genes := gsub(".*:","",symbols)]
  
  if(verbose) message("Defining edges...")
  edges = pairs[similarity_score >= corr_thresh]
  
  # skip if there are no edges
  if(nrow(edges) == 0){
    message(sprintf("No similarity scores > %s.", corr_thresh))
    return()
  } 
  
  # merge back with points
  edges_v1 = merge(edges, points, by.x='V1', by.y='Var', all.x=T)
  edges = merge(edges_v1, points, by.x='V2', by.y='Var', all.x=T, suffixes = c("_Var1", "_Var2"))
  
  # add intersection (gene symbols)
  edges[,intersection := .edge_intersection(genes_Var1, genes_Var2), by=1:nrow(edges)]
  
  # remove nodes with only one supporting gene 
  if(include_metab){
    points = points[grepl(",", intersection_formatted) | intersection_formatted == "unknown (METAB)"]
  }else{
    points = points[grepl(",", intersection_formatted)]
  }
  
  if(nrow(points) == 0){
    message("No significant pathway enrichments driven by more than one gene.")
    return()
  }
  
  if(verbose) message("Final touches...")
  
  #if(verbose) message("Labelling groups by most frequent pathway subclass...")
  
  # add class and subclass to edges 
  points[,pathway_class := classlist[Var]]
  points[,pathway_subclass := gsub(".*; ","",pathway_class)]
  
  # if gene lists are longer than 50 characters, add line breaks 
  # make sure to skip METAB nodes 
  points[,enriched_br := .format_gene_lists(enriched), by = 1:nrow(points)]
  
  # if pathway name + parent is longer than 50 characters, add line breaks 
  points[,title := sprintf("<b>%s</b> (%s)", term_name, pathway_subclass)]
  points[,title_br := .add_line_breaks(title, sep=' '), by = 1:nrow(points)]
  
  # we will color by group later
  viznetwork_nodes = data.table(id = points[,Var],
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
  
  # edge color
  #edge_colors = .get_feature_colors(brewer.pal(9,"YlGnBu")[4:8], edges[,similarity_score])
  
  # add line breaks
  # allow for METAB
  edges[,intersection_label := sprintf("<b>Intersection:</b> %s",intersection)]
  edges[,intersection_br := .add_line_breaks(intersection_label), by = 1:nrow(edges)]
  
  viznetwork_edges = data.table(from = edges[,V1],
                                to = edges[,V2],
                                value = (edges[,similarity_score]*50),
                                title = sprintf("<b>Similarity score:</b> %s<br>%s", # tooltip
                                                round(edges[,similarity_score],2),
                                                edges[,intersection_br]),
                                color.color = "#848484", # visNetwork default
                                physics = F,
                                color.highlight = "#545454",
                                labelHighlightBold = T)
  
  # # remove nodes without an edge before clustering
  # viznetwork_nodes = viznetwork_nodes[id %in% c(viznetwork_edges[,from], viznetwork_edges[,to])]
  
  #if(verbose) message("Identify clusters in the graph...")
  
  # remove edges without any nodes
  viznetwork_edges = viznetwork_edges[from %in% viznetwork_nodes[,id] & to %in% viznetwork_nodes[,id]]
  
  # define group by igraph cluster 
  net = graph_from_data_frame(d=viznetwork_edges, vertices=viznetwork_nodes, directed=F)
  clust_to_pw = igraph::cluster_louvain(net)
  viznetwork_nodes[,group := as.character(clust_to_pw$membership)]
  
  #if(verbose) message("Clean up nodes...")
  
  # remove pathways when they're cluster singletons 
  clust_sizes = as.list(table(clust_to_pw$membership))
  singletons = names(clust_sizes)[clust_sizes==1]
  #viznetwork_nodes = viznetwork_nodes[!group %in% clust_to_remove]
  # put singletons in group 0
  viznetwork_nodes[group %in% singletons, group := '0']
  
  # if(nrow(viznetwork_nodes)==0){
  #   message("No non-singleton pathways.")
  #   return()
  # }
  
  # use most common subclass for group name
  # use these in the legend 
  viznetwork_nodes[,group_number := group]
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
  
  # if(include_metab){
  #   # do we have METAB-only enrichments?
  #   if('METAB' %in% points[,omes]){
  #     # find the group with all of the orphan METAB enrichments 
  #     metab_only_ids = points[omes=="METAB", Var]
  #     metab_group = unique(viznetwork_nodes[id %in% metab_only_ids, group])
  #     if(length(metab_group) > 1){
  #       message(sprintf("I'm surprised there's more than one METAB cluster: %s", paste0(metab_group, collapse=', ')))
  #     }
  #     viznetwork_nodes[group == metab_group, group := sprintf('%s: METAB-specific enrichments', group_number)]
  #   }
  # }
  
  #if(verbose) message("Assigning colors to groups...")
  
  # help by initializing 
  vn0 = visNetwork(viznetwork_nodes, viznetwork_edges) %>%
    visGroups(groupname = viznetwork_nodes[1,group], color = "red", shape = "triangle") 
  
  # define colors 
  # use Set2, then Set3, then Set1 + Set 3
  n_clusters = length(unique(viznetwork_nodes[,group]))
  if(n_clusters > 12 & n_clusters < 21){
    hex = c(brewer.pal(12, "Set3"), brewer.pal(9, "Set1"))
  }else if(n_clusters > 8){
    hex = brewer.pal(12, "Set3")
  }else if(n_clusters <= 8){
    hex = brewer.pal(8, "Set2")
  }else{
    warning("More than 21 clusters? I didn't think that would happen.")
    hex = c(brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), brewer.pal(7, "Set2"))
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
    all_tissues = unique(unname(MotrpacBicQC::tissue_abbr))
    all_tissues = all_tissues[!all_tissues %in% c("VENACV","OVARY","TESTES")]
    all_omes = unique(unname(MotrpacBicQC::assay_abbr))
    all_omes = all_omes[!all_omes %in% c("U-METAB","N-METAB")]
    
    for(col in c(all_tissues, all_omes)){
      if(! col %in% colnames(vn)){
        vn[,(col) := 0]
      }
    }
    
    # add -log10 sumlog p for node size
    vn[,sumlog_p_log10 := round(-log10(sumlog_p),2)]
    
    # clean up nodes and edges
    setnames(vn, 
             c('group_number', 'group'), 
             c('group', 'group_label'))
    
    ve[,title := NULL]
    
    cytoscape_igraph = igraph::graph_from_data_frame(ve, directed = F, vertices = vn)
    return(cytoscape_igraph)
  }
  
  # add group label nodes
  if(add_group_label_nodes){
    # create additional nodes with a strong weight to a random member of that group 
    group_nodes = data.table(id = names(hex),
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
    
    viznetwork_nodes = rbindlist(list(viznetwork_nodes, group_nodes), fill=T)
    
    # add edges
    # for each group, pick node (smallest pval) where the category is also the label
    rowselect = viznetwork_nodes[,grepl(pathway_subclass, group, fixed = TRUE), by=1:nrow(viznetwork_nodes)][,V1]
    # # add row for "METAB-specific enrichments"
    # if(any(grepl("METAB-specific enrichments", viznetwork_nodes[,group]))){
    #   # add a row for this group. connect to any pathway 
    #   # change first METAB row to T
    #   rowselect[grep("METAB-specific enrichments", viznetwork_nodes[,group])[1]] = TRUE
    # }
    
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
    viznetwork_edges = rbindlist(list(viznetwork_edges, group_edges), fill=T)
    
    v1 = visNetwork(viznetwork_nodes, viznetwork_edges, main=title) %>%
      visInteraction(tooltipDelay = 10,
                     tooltipStay = Inf,
                     hideEdgesOnZoom = T) %>%
      visIgraphLayout(layout = "layout_nicely") 
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
    
    v1 = visNetwork(viznetwork_nodes, viznetwork_edges, main=title) %>%
      visLegend(zoom = F,
                addNodes = lnodes,
                useGroups = F,
                width = 0.3,
                stepY = 40) %>%
      visInteraction(tooltipDelay = 10,
                     tooltipStay = Inf,
                     hideEdgesOnZoom = T) %>%
      visIgraphLayout(layout = "layout_nicely") 
  }
  
  if(out_html != "/dev/null"){ 
    if( (file.exists(out_html) & overwrite_html) | !file.exists(out_html)){
      visSave(v1, file = out_html)
      message(sprintf("Network saved in %s", out_html))
      # adjust browser window size 
      cmd = sprintf('sed -i \'s|"browser":{"width":960,"height":500,"padding":40,"fill":false}|"browser":{"width":960,"height":960,"padding":0,"fill":false}|\' %s', out_html)
      system(cmd)
    }
  }
  
  if(verbose){
    message(sprintf("Elapsed time: %s", round((proc.time() - ptm)[3], 3)))
    cat("\n")
  }
  
  if(return_html){
    return(out_html)
  }else{
    return(v1)
  }
}
