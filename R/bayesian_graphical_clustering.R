#' Perform Bayesian graphical clustering of the differential analysis results using repfdr
#' 
#' @import repfdr
#' 
#' @param zscores A numeric matrix. Rows are analytes, columns are conditions (e.g., male week 1). Entries are z-scores.
#' @param min_analyte_posterior_thr A number. The minimal value of sum of posteriors for including an analyte.
#' @param min_prior_for_config A number. The minimal prior probability on the configuration priors.
#' @param naive_edge_sets A logical. TRUE: edge sets are extracted by taking simple intersection of nodes. FALSE: edge sets are extracted from the repfdr solution.
#' 
#' @details 
#' This function is highly specialized for the motrpac data. It assumes that the input zscore matrix has eight columns.
#' Four columns for male z-scores (weeks 1,2,4,8), and four for females. No NA values in the input matrix.
#' The function first runs the Bayesian clustering using repfdr (see repfdr_wrapper) and then extracts the node and edge sets of the graphical solution.
#' 
#' A node is a time-point specific state. It can represent null features (z-score around zero, no effect), up-regulated features (positive z-score), or down regulated features (negative z-scores).
#' An edge represents interaction among two adjacent time point. For example, it can represent features that are up-regulated in week 1 in both males and females, and null in both male and females by week 2.
#' min_analyte_posterior_thr is used to make sure that features without any reasonable fit in the clustering solution will be discarded, whereas min_prior_for_config removes clusters with extremely low prior (i.e., do not have any associated features).
#' 
#' @return A named list with two items: (1) node sets, and (2) edge sets.
bayesian_graphical_clustering<-function(zscores,
      min_analyte_posterior_thr=0.5,min_prior_for_config=0.001,
      naive_edge_sets=T){
  
  if(is.null(dim(zscores))||nrow(zscores)<5 || ncol(zscores)!=8){
    stop("Input zscore matrix does not fit the required input. See documentation for details.")
  }
  
  repfdr_results = repfdr_wrapper(zscores,min_prior_for_config=min_prior_for_config)
  repfdr_clusters = repfdr_results$repfdr_clusters
  repfdr_cluster_posteriors = repfdr_results$repfdr_cluster_posteriors
  
  node_sets = list()
  post_rowsums = rowSums(repfdr_cluster_posteriors,na.rm=T)
  # make sure that features without any reasonable fit in the selected
  # configs will not be included in any solution
  # we do this by simply putting a very large constant value here,
  # this will make sure that these features will not have any
  # high weights in the process below and will therefore be 
  # discarded when we create the graph objects
  post_rowsums[post_rowsums < min_analyte_posterior_thr] = 10^10
  # specify the weeks of each column in zs_smoothed
  study_weeks = c(rep(c("1w","2w","4w","8w"),2))
  for (j in 1:4){
    female_ind = j
    male_ind = j+4 
    week = study_weeks[j]
    for(o1 in c(-1,0,1)){
      for(o2 in c(-1,0,1)){
        curr_configs = repfdr_clusters[female_ind,]==o1 & repfdr_clusters[male_ind,]==o2
        if(sum(curr_configs)==0){next}
        if(sum(curr_configs)==1){
          curr_config_sums = repfdr_cluster_posteriors[,curr_configs]
        }
        else{
          curr_config_sums = rowSums(repfdr_cluster_posteriors[,curr_configs])
        }
        curr_features = rownames(repfdr_cluster_posteriors)[curr_config_sums / post_rowsums > 0.5]
        node_sets[[paste0(week,"_","F",o1,"_","M",o2)]] = curr_features
      }
    }
  }
  # sanity check, we should not see NAs in the node sets
  if(!all(0==sapply(node_sets,function(x)sum(is.na(x))))){
    stop("NA features appear in the node sets")
  }
  # Add week 0 as a dummy node with all features
  node_sets[["0w"]] = rownames(repfdr_cluster_posteriors)
  edge_sets = list()
  if(naive_edge_sets){
    # compute all sets going from a week to its following week using simple 
    # intersection of node sets.
    # below we also provide code to extract edge sets directly from repfdr's output.
    # in practice we observed that the results are almost identical, so we decided to
    # keep this simpler approach.
    weeks = c("0w","1w","2w","4w","8w")
    for(j in 2:length(weeks)){
      this_week_sets = node_sets[grepl(weeks[j],names(node_sets))]
      prev_week_sets = node_sets[grepl(weeks[j-1],names(node_sets))]
      for(n1 in names(this_week_sets)){
        for(n2 in names(prev_week_sets)){
          edge_sets[[paste(n2,n1,sep="---")]] = intersect(node_sets[[n1]],node_sets[[n2]])
        }
      }
    }
  }else{
    edge_sets2 = list()
    for (j in 1:3){
      female_ind = j
      male_ind = j+4
      next_female_ind = female_ind+1
      next_male_ind = male_ind + 1
      week = study_weeks[j]
      next_week = study_weeks[j+1]
      for(o1 in c(-1,0,1)){
        for(o2 in c(-1,0,1)){
          for(o3 in c(-1,0,1)){
            for(o4 in c(-1,0,1)){
              curr_configs = repfdr_clusters[female_ind,]==o1 & repfdr_clusters[male_ind,]==o2 &
                repfdr_clusters[next_female_ind,]==o3 & repfdr_clusters[next_male_ind,]==o4
              if(sum(curr_configs)==0){next}
              if(sum(curr_configs)==1){
                curr_config_sums = repfdr_cluster_posteriors[,curr_configs]
              }
              else{
                curr_config_sums = rowSums(repfdr_cluster_posteriors[,curr_configs])
              }
              curr_features = rownames(repfdr_cluster_posteriors)[curr_config_sums / post_rowsums > 0.5]
              edge_name = paste0(week,"_","F",o1,"_","M",o2,"---",next_week,"_","F",o3,"_","M",o4)
              edge_sets2[[edge_name]] = curr_features
            }
          }
        }
      }
    }
    edge_sets = edge_sets2
  }
  return(list(
    node_sets=node_sets,edge_sets=edge_sets
  ))
}


#' A general wrapper for running repfdr on a matrix of z-scores
#' 
#' @import repfdr,MotrpacRatTraining6moData
#' 
#' @param zscores A numeric matrix. Rows are analytes, columns are conditions (e.g., male week 1). Entries are z-scores.
#' @param min_prior_for_config A number. The minimal prior probability on the configuration priors.
#' 
#' @return A named list with two objects: (1) a binary matrix representation of the clusters, and (2) a matrix with the cluster posteriors.
repfdr_wrapper<-function(zscores,min_prior_for_config = 0.001){
  # Step 1 in repfdr: discretize the z-scores and estimate the probabilities 
  # in each bin.
  # In our work we examined the diagnostic plots manually. The results look reasonable,
  # especially for the qqplots, so we decided to use the current estimation with df=20
  ztobins_res = ztobins(zscores,df=20,type=1,n.bins=150,central.prop = 0.25,
                        plot.diagnostics = F)
  # Step 2 in repfdr: estimate the repfdr model using the EM algorithm
  repfdr_res = repfdr(ztobins_res$pdf.binned.z,
                      ztobins_res$binned.z.mat,non.null = 'replication',
                      control = em.control(max.iter = 500,tol=1e-06,nr.threads = 4),
  )
  # sanity check
  configs = hconfigs(ncol(zscores))
  if(!all(configs == repfdr_res$Pi[,1:ncol(zscores)])){
    stop("Internal error in repfdr, try running the repfdr package using its tutorial")
  }

  # Reorganize the repfdr output for downstream analyses and add the 
  # posterior estimates for every analyte in our dataset
  rownames(repfdr_res$Pi) = apply(repfdr_res$Pi[,1:ncol(zscores)],1,paste,collapse="")
  # We next extract posterior probabilities for clusters with pi > min_prior_for_config
  # This limits the analysis to the configurations that we believe to be
  # well represented in the dataset
  hvec_inds = which(repfdr_res$Pi[,"Pi"] >= min_prior_for_config)
  hvecs = cbind(hvec_inds,repfdr_res$Pi[hvec_inds,"Pi"])
  # the following call returns the posterior probabilities Pr(h|v)
  # for every configuration h in our selected set of configs, and every z-scores
  # vector v in the dataset
  repfdr_posteriors = ldr(
    ztobins_res$pdf.binned.z,ztobins_res$binned.z.mat,
    repfdr_res$Pi,h.vecs = hvec_inds
  )
  # transform the output to objects that will be used below for 
  # the graphical analysis
  # this object keeps the configurations h in a matrix
  repfdr_clusters = t(repfdr_posteriors[,1:ncol(zscores)])
  # this is a simple string representation of the configs
  repfdr_clusters_str = apply(repfdr_clusters,2,paste,collapse="")
  # here we keep P(h) the prior distribution of the configs
  repfdr_clusters_pi = repfdr_res$Pi[repfdr_clusters_str,"Pi"]
  # this matrix keeps the posteriors P(h|v)
  repfdr_cluster_posteriors = t(repfdr_posteriors[,-c(1:ncol(zscores))])
  # to make things easier and explicit we keep the string names of the configs here
  colnames(repfdr_cluster_posteriors) = repfdr_clusters_str
  
  return(list(
    repfdr_clusters = repfdr_clusters,
    repfdr_cluster_posteriors = repfdr_cluster_posteriors
  ))
}


# set our own layout using ggraph
# see http://blog.schochastics.net/post/ggraph-tricks-for-common-problems/

# when using source, this will fail if one of the libraries was not installed
library(ggraph)
library(scatterpie)
library(grid)
library(networkD3)
library(repfdr)
library(data.table)
library(igraph)

#' An auxiliary function for reducing a list of sets by a regex.
#' Useful for filtering DEA analyte sets by tissues or omes.
#' 
#' @params l a named list of character vectors
#' @params regs a character vector of regular expressions
limit_sets_by_regex<-function(sets,regs){
  if(is.null(regs) || length(regs)==0){return(sets)}
  l = list()
  for(nn in names(sets)){
    v = sets[[nn]]
    newv = c()
    for(r in regs){
      r = paste0(r,";")
      newv = union(newv,v[grepl(r,v)])
    }
    l[[nn]] = newv
  }
  return(l)
}

#' Auxiliary function to get the largest paths in a graph.
#' 
#' This is implemented using a dynamic programming approach where
#' we iteratively add the data of the next edge.
#' 
#' When analyzing an edge (x,y) with a set of analytes s, we extend all
#' trajectories that end with x using the new edge, but also all the 
#' trajectories that start with y. 
#' 
#' At the end, because we examine all edges and all extensions we are
#' guaranteed to have covered all full paths.
#' 
#' The min_size parameter is important, since we are interested in paths
#' and edges with at least this number of analytes, then we know that if a
#' current trajectory does not have at least min_size analytes, then since the 
#' set of any extension can only be the same or smaller then we can ignore
#' such paths moving forward. 
#' 
#' However, this parameter has to be considered with care, as specifying a number
#' that is too high will result in no paths in the output.
#' 
#' @param edge_sets a named list of string vectors. The name of an edge is node_id---node_id
#'        edges with no analytes have a NULL set (a set of size zero, but are still represented),
#'        node ids are time_points_Fx_My where x and y represent the up/down state in each sex.
#' @param min_size a number, specifying the minimal path size to be considered
#' @return NULL if no paths of size of at least min_size were found, otherwise
#'         return a data frame that represents all paths of size min_size or greater,
#'         ranked from the largest path to the smallest one.
get_trajectory_sizes_from_edge_sets<-function(edge_sets,min_size=10){
  node_names = unique(unlist(strsplit(names(edge_sets),split="---")))
  node_weeks = sapply(node_names,function(x)strsplit(x,split="_")[[1]][1])
  full_path_size = length(unique(node_weeks))
  arrs = strsplit(names(edge_sets),split="---")
  prevs = sapply(arrs,function(x)x[1])
  nexts = sapply(arrs,function(x)x[2])
  l = list()
  l_sets = c()
  for(j in 1:length(edge_sets)){
    if(length(edge_sets[[j]])<min_size){next}
    ind = length(l)+1
    curr_prev = prevs[j]
    curr_next = nexts[j]
    l[[ind]] = c(curr_prev,curr_next)
    l_sets[[ind]] = edge_sets[[j]]
    for(j2 in 1:length(l)){
      j2_arr = l[[j2]]
      # here the path of j2 starts with the current edge next
      # so the edge represented by j can extend the path j2
      # "on the left" as it represent a previous time interval
      if(j2_arr[1] == curr_next){
        new_ind = length(l)+1
        l[[new_ind]] = c(curr_prev,j2_arr)
        l_sets[[new_ind]] = intersect(l_sets[[j2]],edge_sets[[j]])
      }
      # here the path of j2 end with the current edge prev/start
      # so the edge represented by j can extend the path j2
      # "on the right" as it represent a subsequent time interval
      if(j2_arr[length(j2_arr)]==curr_prev){
        new_ind = length(l)+1
        l[[new_ind]] = c(j2_arr,curr_next)
        l_sets[[new_ind]] = intersect(l_sets[[j2]],edge_sets[[j]])
      }
    }
    keep = sapply(l_sets,length) >= min_size
    l = l[keep]
    l_sets = l_sets[keep]
  }
  # limit the results to the "full" trajectories
  traj_sizes = sapply(l,length)
  # here we assume that at least one full pathway survived
  # the min size filter above, otherwise, return NULL
  traj_inds = traj_sizes==full_path_size
  if(sum(traj_inds)==0){return(NULL)}
  l = l[traj_inds]
  l_sets = l_sets[traj_inds]
  trajectories = cbind(
    t(sapply(l,function(x)x)),
    sapply(l_sets,length)
  )
  trajectories = data.frame(trajectories,stringsAsFactors = F)
  trajectories[[6]] = as.numeric(trajectories[[6]])
  trajectories = trajectories[order(trajectories[[6]],decreasing = T),]
  return(trajectories)
}

#' Keep the edges of the top trajectories of an edge set
filter_edge_sets_by_trajectories<-function(edge_sets,topk=5,min_path_size=5){
  
  traj = get_trajectory_sizes_from_edge_sets(edge_sets,min_size = min_path_size)
  edges_to_keep = c()
  if(!is.null(traj)){
    topk = min(topk,nrow(traj))
    traj = traj[1:topk,]
    for(j in 1:(ncol(traj)-1)){
      curr_edges = paste(traj[,j],traj[,j+1],sep="---")
      edges_to_keep = union(edges_to_keep,curr_edges)
    }
  }
  
  e_copy = copy(edge_sets)
  for(e in names(e_copy)){
    if(!(e %in% edges_to_keep)){
      e_copy[[e]] = character(0) 
    }
  }
  return(e_copy)
}

#' The main function for obtaining a graphical (tree) representation of the differential
#' analysis results.
#' 
#' @param tissues acharacter vector. The set of tissues to take for the analysis. If null take all.
#' @param omes acharacter vector. The set of omes to take for the analysis. If null take all.
#' @param node_sets a named list with the node (state) sets of analyte, see details for analyte name convention.
#' @param edge_sets a named list with the edge (state) sets of analyte, see details for analyte name convention.
#' @param min_size a numeric, a threshold on the set sizes to be considered
#' @param parallel_edges_by_ome a logical. TRUE means that we want to added parallel edges for the different omes.
#' @param parallel_edges_by_tissue a logical. TRUE means that we want to added parallel edges for the different tissues.
#' @param edge_width_range a numeric vector of size 2, a parameter for ggraph
#' @param edge_alpha_range a numeric vector of size 2, a parameter for ggraph
#' @param color_nodes_by_states a logical. If TRUE, nodes are colored by states. 
#'        Red for up-reg, blue for down-reg, green for a discrepancy between the sexes.
#' @param max_trajectories a numeric or NULL. If not NULL then it specifies the number of pathways
#'        to keep when looking into the edge sets after filtering by omes and tissues. If 
#'        parallel_edges_by_tissue = TRUE then take the top trajectories in each tissue.
#' @param highlight_subset a character string or NULL. If not NULL then it specifies the name of a 
#'        node, edge, or path to highlight in the tree. 
#' 
#' 
#' @details 
#'         The function filters the input set to include analytes from the given tissues
#'         and omes (if tissues/omes are not null).
#'         If parallel edges are requested then the relevant edge sizes are computed internally
#'         and are used within ggraph for the output plot.
#'         Analyte names are in the ome;tissue;feature_id format.
#'         Node set names are are in the week_number(1,2,4,8)'w'_F{-1,0,1}_M{-1,0,1} format
#'         for example: 1w_F-1_M0
#'         Edge set names are in the node_a---node_b format, e.g., 4w_F0_M0---8w_F0_M1 
get_tree_plot_for_tissue<-function(
  tissues,
  omes = NULL,
  node_sets,
  edge_sets,
  min_size = 20,
  parallel_edges_by_ome = F,
  parallel_edges_by_tissue = F,
  edge_width_range = c(0,10),
  edge_alpha_range = c(0,1),
  color_nodes_by_states = T,
  max_trajectories = NULL,
  highlight_subset = NULL,
  curvature = 0.1
){
  
  tissues = unique(tissues)
  omes = unique(omes)
  
  tissue_node_sets = limit_sets_by_regex(node_sets,tissues)
  tissue_node_sets = limit_sets_by_regex(tissue_node_sets,omes)
  tissue_edge_sets = limit_sets_by_regex(edge_sets,tissues)
  tissue_edge_sets = limit_sets_by_regex(tissue_edge_sets,omes)
  
  # fix edge width range
  # should not plot 0-width edges
  # not sure this does anything?
  if(edge_width_range[1] == 0){
    edge_width_range[1] = 1e-5
  }
  
  # Set things up to highlight a subset of the tree
  highlight_edge = F
  highlight_node = F
  highlight_path = F
  if(!is.null(highlight_subset)){
    if(grepl(":", highlight_subset)){
      # Check if it's this tissue 
      if(grepl(":", highlight_subset)){
        if(length(tissues)>1){
          warning("'highlight_subset' cannot handle a subset of tissues. Setting 'highlight_subset' to NULL Specify a single tissue or remove the tissue prefix from 'highlight_subset'.")
          highlight_subset = NULL
        }else{
          highlight_tissue = gsub(":.*","",highlight_subset)
          highlight_cluster = gsub(".*:","",highlight_subset)
          if(highlight_tissue != tissues){
            warning("'highlight_subset' tissue prefix does not match the tissue supplied. Setting 'highlight_subset' to NULL.")
            highlight_subset = NULL
          }
        }
      }else{
        highlight_cluster = highlight_subset
      }
    }
    if(!is.null(highlight_subset)){
      if(grepl("---", highlight_subset)){
        highlight_edge = T
      }else if(grepl("->",highlight_subset)){
        highlight_path = T
      }else{
        highlight_node = T
      }
    }
  }
  
  if(!is.null(highlight_subset) & any(parallel_edges_by_ome, parallel_edges_by_tissue, color_nodes_by_states)){
    warning("Setting 'parallel_edges_by_ome', 'parallel_edges_by_tissue', and 'color_nodes_by_states' to FALSE because 'highlight_subset' is TRUE.")
    parallel_edges_by_ome = F
    parallel_edges_by_tissue = F
    color_nodes_by_states = F
  }
  
  # Check the input omes and tissues sets
  # If either ome or tissues is NULL, take them from the data 
  all_analytes = unique(unlist(tissue_node_sets)) # just this subset of the data 
  all_analytes_info = strsplit(all_analytes,split=";")
  data_tissue_set = unique(sapply(all_analytes_info,function(x)x[2]))
  data_ome_set = unique(sapply(all_analytes_info,function(x)x[1]))
  if(!all(tissues %in% data_tissue_set)){
    warning("Some input tissue names are not in the dataset and will be removed")
    tissues = intersect(tissues,data_tissue_set)
  }
  if(!all(omes %in% data_ome_set)){
    warning("Some input ome names are not in the dataset and will be removed")
    omes = intersect(omes,data_ome_set)
  }
  if(is.null(tissues)){
    tissues = data_tissue_set
  }
  if(is.null(omes)){
    omes = data_ome_set
  }
  
  # if we are not asked to get parallel edges and are asked to take 
  # the top trajectories then we can simply filter the edge sets and move on
  if(!is.null(max_trajectories) && !parallel_edges_by_tissue){
    tissue_edge_sets = filter_edge_sets_by_trajectories(tissue_edge_sets,max_trajectories)
  }
  
  # Create the initial matrices/data frame
  # with some additional info (depends on the user input), these will be used as input
  # for generating the graph
  node_info = cbind(
    names(tissue_node_sets),
    sapply(names(tissue_node_sets),function(x)strsplit(x,split="_")[[1]][1]),
    sapply(tissue_node_sets,length)
  )
  colnames(node_info) = c("node","week","size")
  edge_info = cbind(
    sapply(names(tissue_edge_sets),function(x)strsplit(x,split="---")[[1]][1]),
    sapply(names(tissue_edge_sets),function(x)strsplit(x,split="---")[[1]][2]),
    sapply(tissue_edge_sets,length)
  )
  colnames(edge_info) = c("prev","next","size")
  
  # If we are asked to add parallel edges by omes then we need to 
  # compute ome edge sizes and add them
  if(parallel_edges_by_ome){
    for(o in omes){
      v = c()
      for(e in rownames(edge_info)){
        v[e] = sum(grepl(paste0(o,";"),tissue_edge_sets[[e]]))
      }
      edge_info = cbind(edge_info,v)
      colnames(edge_info)[ncol(edge_info)] = o
    }
    # always use the same assay colors
    edge_colors = define_assay_cols() # defined in cluster_viz_fx.R
    # require(RColorBrewer)
    # o = c("METAB","TRNSCRPT","PROT","ACETYL","PHOSPHO","UBIQ","ATAC","METHYL","IMMUNO")
    # edge_colors = brewer.pal(length(o), 'Set1')
    # names(edge_colors) = o
  }
  # if parallel_edges_by_tissue = TRUE then:
  # if we do not need to filter the top trajectories then we can simply count
  # the number of features per edge for a tissue using regular expressions
  if(is.null(max_trajectories) && parallel_edges_by_tissue){
    for(tissue in tissues){
      v = c()
      for(e in rownames(edge_info)){
        v[e] = sum(grepl(paste0(tissue,";"),tissue_edge_sets[[e]]))
      }
      edge_info = cbind(edge_info,v)
      colnames(edge_info)[ncol(edge_info)] = tissue
    }
    # always use the same assay colors
    edge_colors = MotrpacBicQC::tissue_cols
  }
  # if parallel_edges_by_tissue = TRUE then:
  # if we are asked to filter by the top trajectories then we need to 
  # recompute the edge sizes per tissue by limiting the tissue edge sets
  if(!is.null(max_trajectories) && parallel_edges_by_tissue){
    for(tissue in tissues){
      curr_tissue_edge_sets = limit_sets_by_regex(tissue_edge_sets,tissue)
      curr_tissue_edge_sets = filter_edge_sets_by_trajectories(
        curr_tissue_edge_sets,max_trajectories)
      v = c()
      for(e in rownames(edge_info)){
        v[e] = length(curr_tissue_edge_sets[[e]])
      }
      edge_info = cbind(edge_info,v)
      colnames(edge_info)[ncol(edge_info)] = tissue
    }
    # always use the same assay colors
    edge_colors = MotrpacBicQC::tissue_cols
  }
  
  # transform the initial node/edge info matrices to data frames
  # filter edge sizes by the specified user input
  d = as.data.frame(edge_info)
  for(j in 3:ncol(d)){
    d[[j]] = as.numeric(d[[j]])
    d[d[,j]<min_size,j] = 0 # filter by min edge size
  }
  d_nodes = data.frame(node_info,check.names = F)
  d_nodes[[3]] = as.numeric(d_nodes[[3]])
  d_nodes$inds = 0:(nrow(d_nodes)-1)
  
  # If we need to plot parallel edges then we need to transform the edge-based
  # data frame into a long format (to make it work with ggraph)
  if(parallel_edges_by_ome || parallel_edges_by_tissue){
    added_cols = names(d)[-c(1:3)]
    if(!parallel_edges_by_ome){added_cols = setdiff(added_cols,omes)}
    if(!parallel_edges_by_tissue){added_cols = setdiff(added_cols,tissues)}
    newd = d[,1:2]
    newd$size = 0
    newd$type = NA
    for(o in added_cols){
      currd = d[,c(1:2,which(colnames(d)==o))]
      names(currd)[3] = "size"
      currd$type = o
      newd = rbind(newd,currd)
    }
    newd = newd[!is.na(newd$type),]
    d_g = igraph::graph_from_data_frame(newd)
    if(all(E(d_g)$size == 0)){
      message("No non-0 edges. Skipping.")
      return()
    }
    E(d_g)$edge_size = edge_attr(d_g,"size")
    d_g_auto_layout <- create_layout(d_g, layout = 'auto')
  }else{
    d_g = igraph::graph_from_data_frame(d)
    if(all(E(d_g)$size == 0)){
      message("No non-0 edges. Skipping.")
      return()
    }
    E(d_g)$edge_size = edge_attr(d_g,"size")
    d_g_auto_layout <- create_layout(d_g, layout = 'auto')
  }
  
  ################
  # From this point we start generating ggraph info
  # (e.g., layout, labels etc.)
  d_g_ordered_nodes = c("F-1_M1","F1_M-1","F0_M-1","F-1_M-1","F-1_M0",
                        "F0_M0","F0_M1","F1_M1","F1_M0")
  d_g_ordered_cols = c("lightgreen","lightgreen","blue","blue","blue",
                       "gray","red3","red3","red3")
  #d_g_ordered_cols_alt = c("black","black","black","black","black",
  #                     "gray","black","black","black")
  d_g_ordered_cols_alt = c("gray4","gray4","gray4","gray4","gray4",
                           "gray","gray4","gray4","gray4")
  d_g_ordered_shapes = c(18,18,19,19,19,19,19,19,19)
  layer_plot_names = c(
    "F down, M up","F up, M down",
    "M only down","Both down","F only down",
    "No\nchange",
    "M only up","Both up","F only up"
  )
  
  names(d_g_ordered_cols) = d_g_ordered_nodes
  names(d_g_ordered_shapes) = d_g_ordered_nodes
  names(d_g_ordered_cols_alt) = d_g_ordered_nodes
  l_x_lim = c(min(d_g_auto_layout$x),max(d_g_auto_layout$x))
  l_y_lim = c(min(d_g_auto_layout$y),max(d_g_auto_layout$y))
  d_g_our_layout  = copy(d_g_auto_layout)
  # set 5 horiz layers over time
  xjump = (l_x_lim[2]-l_x_lim[1]) / 4
  # set 9 vertical layers over ordered_nodes
  yjump = (l_y_lim[2]-l_y_lim[1]) / 8
  d_g_our_layout[d_g_our_layout$name=="0w","x"] = l_x_lim[1]
  d_g_our_layout[d_g_our_layout$name=="0w","y"] = mean(l_y_lim)
  curr_weeks = c("1w","2w","4w","8w")
  for(j in 1:length(curr_weeks)){
    w = curr_weeks[j]
    d_g_our_layout[grepl(w,d_g_our_layout$name),"x"] = l_x_lim[1] + j*xjump
  }
  for(j in 1:length(d_g_ordered_nodes)){
    n = d_g_ordered_nodes[j]
    d_g_our_layout[grepl(n,d_g_our_layout$name),"y"] = l_y_lim[1] + (j-1)*yjump
  }
  # make sure that the 0w node is in the same line as of the no response
  d_g_our_layout[1,"y"] = d_g_our_layout[grepl("F0_M0",d_g_our_layout$name),"y"][1]
  # set node sizes and other features
  V(d_g)$setsize = d_nodes[V(d_g)$name,"size"]
  V(d_g)$setsize[V(d_g)$name == "0w"] = median(V(d_g)$setsize)
  V(d_g)$label = sapply(V(d_g)$name,
                        function(x){a=strsplit(x,split="w_")[[1]];a[length(a)]})
  V(d_g)$label = gsub("_","\n",V(d_g)$label)
  V(d_g)$col = "gray"
  V(d_g)$alt_col = "gray"
  V(d_g)$shape = 15
  for(j in 1:length(d_g_ordered_nodes)){
    n = d_g_ordered_nodes[j]
    V(d_g)$col[grepl(n,d_g_our_layout$name)] = d_g_ordered_cols[n]
    V(d_g)$alt_col[grepl(n,d_g_our_layout$name)] = d_g_ordered_cols_alt[n]
    V(d_g)$shape[grepl(n,d_g_our_layout$name)] = d_g_ordered_shapes[n]
  }
  
  # Add the set size as a field in the layout data frame
  # This will be used for controlling the node size correctly
  d_g_our_layout["size"] = V(d_g)$setsize[d_g_our_layout$.ggraph.orig_index]
  
  # for(ome in colnames(tissue_ome_data)){
  #   d_g = set_vertex_attr(d_g,ome,V(d_g),tissue_ome_data[V(d_g)$name,ome])
  # }
  grid_group_annotation_y = c()
  for(j in 1:length(d_g_ordered_nodes)){
    n = d_g_ordered_nodes[j]
    grid_group_annotation_y[n] = l_y_lim[1] + (j-1)*yjump
  }
  weeks_annotation_x = c()
  for(j in 1:length(curr_weeks)){
    w = curr_weeks[j]
    weeks_annotation_x[w] = l_x_lim[1] + j*xjump
  }
  
  # Set colors to highlight edge, node, path
  # Red to highlight; gray for everything else 
  if(highlight_edge){
    V(d_g)$col = "gray"
    V(d_g)$alt_col = "gray"
    # Make edge red
    ecols = rep("gray", nrow(d))
    ecols[which(rownames(d)==highlight_cluster)] = "red"
    E(d_g)$col = ecols
  }else if(highlight_node){
    E(d_g)$col = "gray"
    vcols = rep("gray", length(names(V(d_g))))
    vcols[which(names(V(d_g))==highlight_cluster)] = "red"
    V(d_g)$col = vcols
    V(d_g)$alt_col = vcols
  }else if(highlight_path){
    # extract nodes
    curr_nodes = c("0w",unname(unlist(strsplit(highlight_cluster, "->"))))
    # extract edges
    edge_vector = c()
    for(e in 2:length(curr_nodes)){
      edge_vector = c(edge_vector, sprintf("%s---%s", curr_nodes[e-1], curr_nodes[e]))
    }
    # if 8w is not specified in the path, include all edges to 8w
    if(!any(grepl("8w", curr_nodes))){
      edge_vector = c(edge_vector, rownames(d)[grepl(sprintf("%s---8w", curr_nodes[grepl("4w", curr_nodes)]),rownames(d))])
    }
    
    # set colors for edges
    ecols = rep("gray", nrow(d))
    ecols[rownames(d)%in%edge_vector] = "red"
    E(d_g)$col = ecols
    # set colors for nodes
    vcols = rep("gray", length(names(V(d_g))))
    vcols[names(V(d_g))%in%curr_nodes] = "red"
    V(d_g)$col = vcols
    V(d_g)$alt_col = vcols
  }
  
  if(parallel_edges_by_ome || parallel_edges_by_tissue){
    # Assign NA color to 0-weight edges
    E(d_g)$type[E(d_g)$edge_size == 0] = NA
    p = ggraph(d_g,layout = d_g_our_layout,) +  
      geom_edge_fan0(
        aes(width = E(d_g)$edge_size,
            alpha=E(d_g)$edge_size,
            colour = E(d_g)$type),
        lineend = "round",strength = curvature) +
      scale_edge_color_manual(values=edge_colors, name="", limits=names(edge_colors)[names(edge_colors) %in% c(tissues, omes)]) +
      guides(edge_color=guide_legend(override.aes = list(edge_width=3)))
  }else if(!is.null(highlight_subset)){
    # Assign NA color to 0-weight edges
    E(d_g)$col[E(d_g)$edge_size == 0] = NA
    p = ggraph(d_g,layout = d_g_our_layout,) +  
      #geom_edge_link(aes(width = E(d_g)$edge_size,alpha=E(d_g)$edge_size)) + 
      geom_edge_arc(aes(width = E(d_g)$edge_size,
                        alpha = E(d_g)$edge_size, 
                        color = E(d_g)$col),
                    lineend = "round",strength = curvature) +
      scale_edge_colour_identity(guide = "none")
  }else{
    # Assign NA color to 0-weight edges
    E(d_g)$col = "black"
    E(d_g)$col[E(d_g)$edge_size == 0] = NA
    p = ggraph(d_g,layout = d_g_our_layout,) +  
      #geom_edge_link(aes(width = E(d_g)$edge_size,alpha=E(d_g)$edge_size)) + 
      geom_edge_arc(aes(width = E(d_g)$edge_size,
                        alpha = E(d_g)$edge_size,
                        color = E(d_g)$col),
                    lineend = "round",strength = curvature) +
      scale_edge_color_identity(guide = "none")
  }
  
  p = p +
    scale_edge_width(range=edge_width_range,name="Intersect size") + 
    scale_edge_alpha(range=edge_alpha_range,name="Intersect size") 
  
  if(color_nodes_by_states){
    p = p + geom_node_point(aes(size = size),alpha=1,
                            color=V(d_g)$col,shape=V(d_g)$shape)
  }else{
    p = p + geom_node_point(aes(size = size),alpha=1,
                            color=V(d_g)$alt_col,shape=V(d_g)$shape)
  }
  p = p +
    #geom_node_text(aes(label = V(d_g)$label), repel=F) + 
    scale_size(range = c(2,20),name="Number of analytes") +
    theme(panel.background = element_rect(fill="white"),
          legend.key.size = unit(0.4, 'cm'))
  
  names(grid_group_annotation_y) = layer_plot_names
  for(g in names(grid_group_annotation_y)){
    p = p + ggplot2::annotate(geom="text",x=l_x_lim[1] + 0.2, y=grid_group_annotation_y[g],
                              label=g,fontface = "bold")
  }
  names(weeks_annotation_x) = paste("week",c(1,2,4,8))
  for(w in names(weeks_annotation_x)){
    p = p + ggplot2::annotate(geom="text",y=l_y_lim[2]+0.5, x=weeks_annotation_x[w],
                              label=w,fontface = "bold")
  }
  
  # added above
  # # change colors and legend edge size
  # if(parallel_edges_by_ome | parallel_edges_by_tissue){
  #   p = p +
  #     scale_edge_color_manual(values=edge_colors, name="") +
  #     guides(edge_color=guide_legend(override.aes = list(edge_width=3)))
  # }
  
  return(p)
}


#' An auxiliary function for getting the top analyte sets for a set of tissues
extract_tissue_sets<-function(tissues,node_sets,edge_sets,k=3,
                              min_size=20,add_week8=T,omes=NULL){
  
  tissue_node_sets = limit_sets_by_regex(node_sets,tissues)
  tissue_node_sets = limit_sets_by_regex(tissue_node_sets,omes)
  tissue_edge_sets = limit_sets_by_regex(edge_sets,tissues)
  tissue_edge_sets = limit_sets_by_regex(tissue_edge_sets,omes)
  l = list()
  
  # node sets
  node_set_sizes = sapply(tissue_node_sets,length)
  selected_node_sets = names(sort(node_set_sizes,decreasing = T))[1:(k+1)]
  if(add_week8){
    selected_node_sets = union(selected_node_sets,
                               names(tissue_node_sets)[grepl("8w",names(tissue_node_sets))])
  }
  selected_node_sets = setdiff(selected_node_sets,"0w")
  l[selected_node_sets] = tissue_node_sets[selected_node_sets]
  
  # edge sets
  edge_set_sizes = sapply(tissue_edge_sets,length)
  edge_set_sizes = edge_set_sizes[!grepl("0w",names(edge_set_sizes))]
  selected_edge_sets = names(sort(edge_set_sizes,decreasing = T))[1:k]
  l[selected_edge_sets] = tissue_edge_sets[selected_edge_sets]
  
  # path sets: get the top trajectories first
  tissue_top_trajs = get_trajectory_sizes_from_edge_sets(tissue_edge_sets,min_size = min_size)
  if(!is.null(tissue_top_trajs) && nrow(tissue_top_trajs)>0){
    k_paths = min(k,nrow(tissue_top_trajs))
    tissue_top_trajs = tissue_top_trajs[1:k_paths,]
    for(path_i in 1:nrow(tissue_top_trajs)){
      curr_set = tissue_node_sets[[tissue_top_trajs[path_i,2]]]
      for(path_j in 3:5){
        curr_set = intersect(curr_set,
                             tissue_node_sets[[tissue_top_trajs[path_i,path_j]]])
      }
      if(length(curr_set) != tissue_top_trajs[path_i,6]){
        warning("Computed trajectory set size is different from the precomputed size, 
                this happns if edge sets are not simple intersections of node sets")
      }
      l[[paste(tissue_top_trajs[path_i,2:5],collapse="->")]] = curr_set
    }
  }
  l = l[sapply(l,length)>min_size]
  return(l)
}


#' Pull out all non-empty trajectories 
#' 
#' @param node_sets use `load_graph_vis_data()$node_sets`
#' @param edge_sets use `load_graph_vis_data()$edge_sets`
#' @param tissues string vector, optional. tissue subset. all tissues by default
#' @param omes string vector, optional. ome subset. all omes by default
#' 
#' @return named list with one element per trajectories. members are features in the path 
get_all_trajectories = function(edge_sets, 
                                node_sets, 
                                tissues = unique(unname(MotrpacBicQC::tissue_abbr)),
                                omes = unique(unname(MotrpacBicQC::assay_abbr))){
  require("MotrpacBicQC")
  
  l = list()
  tissue_node_sets = limit_sets_by_regex(node_sets,tissues)
  tissue_node_sets = limit_sets_by_regex(tissue_node_sets,omes)
  tissue_edge_sets = limit_sets_by_regex(edge_sets,tissues)
  tissue_edge_sets = limit_sets_by_regex(tissue_edge_sets,omes)
  tissue_top_trajs = get_trajectory_sizes_from_edge_sets(tissue_edge_sets,min_size = 1)
  if(length(tissue_top_trajs) == 0){
    return()
  }
  for(path_i in 1:nrow(tissue_top_trajs)){
    curr_set = tissue_node_sets[[tissue_top_trajs[path_i,2]]]
    for(path_j in 3:5){
      curr_set = intersect(curr_set,
                           tissue_node_sets[[tissue_top_trajs[path_i,path_j]]])
    }
    if(length(curr_set) != tissue_top_trajs[path_i,6]){
      warning("Computed trajectory set size is different from the precomputed size, 
                  this happns if edge sets are not simple intersections of node sets")
    }
    l[[paste(tissue_top_trajs[path_i,2:5],collapse="->")]] = curr_set
  }
  return(l)
}

