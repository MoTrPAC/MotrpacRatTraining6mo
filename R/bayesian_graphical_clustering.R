#' Bayesian graphical clustering
#' 
#' Perform Bayesian graphical clustering of the differential analysis results using \code{repfdr}
#' 
#' @export
#' 
#' @param zscores A numeric matrix. Rows are analytes, columns are conditions (e.g., male week 1). Entries are z-scores.
#' @param min_analyte_posterior_thr A number. The minimal value of sum of posteriors for including an analyte.
#' @param min_prior_for_config A number. The minimal prior probability on the configuration priors.
#' @param df Integer. Parameter for [repfdr::ztobins()]. Degrees of freedom for fitting the estimated density f(z). 
#' @param naive_edge_sets A logical. TRUE: edge sets are extracted by taking simple intersection of nodes. 
#'   FALSE: edge sets are extracted from the repfdr solution.
#' 
#' @details 
#' This function is highly specialized for the MoTrPAC data. It assumes that the input zscore matrix has eight columns.
#' Four columns for male z-scores (weeks 1, 2, 4, 8), and four for females. No NA values in the input matrix.
#' The function first runs the Bayesian clustering using \code{repfdr} (see [repfdr_wrapper()]) and then extracts the node and 
#' edge sets of the graphical solution.
#' 
#' A node is a time-point specific state. It can represent null features (z-score around zero, no effect), up-regulated 
#' features (positive z-score), or down regulated features (negative z-scores).
#' An edge represents interaction among two adjacent time point. For example, it can represent features that are up-regulated 
#' in week 1 in both males and females, and null in both male and females by week 2.
#' \code{min_analyte_posterior_thr} is used to make sure that features without any reasonable fit in the clustering solution will be discarded, 
#' whereas \code{min_prior_for_config} removes clusters with extremely low prior (i.e., do not have any associated features).
#' 
#' Formally, let zscores be the timewise nxm matrix of z-scores over n analytes and m conditions, where m is even 
#' and has the same time points for males and females. 
#' Each z-score corresponds to the effect of training compared to untrained baseline. 
#' We assume that z-scores have a latent configuration in {-1,0,1}, with -1 denotes down-regulation, 
#' 0 denotes null (no change), and 1 denotes up-regulation. 
#' Let h denote the latent configuration of an analyte over all m conditions. 
#' We use the expectation-maximization (EM) process of the repfdr algorithm (Heller and Yekutieli, 2014) 
#' to estimate the prior π(h) and posterior Pr(h|z_i), for every analyte i. 
#' Once the algorithm converges, we discard configurations h with π(h) < min_analyte_posterior_thr and normalize Pr(h|zi) 
#' to sum to 1 (i.e., all posteriors given the same zi). The new posteriors that can be interpreted as a soft clustering solution, 
#' where the greater the value is, the more likely it is for analyte i to participate in cluster h. 
#' 
#' Given the posteriors Pr(h|zi), we assign analytes to nine possible states in each time point (male or females x time points 1w, 2w,4w, and 8w). 
#' Analyte i belongs to a state if the sum of posteriors that are consistent with that state is > min_analyte_posterior_thr.
#' For every pair of states from adjacent time points j and j+1 we define their edge set as the set of analytes whose 
#' differential expression pattern is consistent with the two nodes.. 
#' Thus, the node and edge sets explained above together define a tree structure that represent different differential patterns over sex and time.
#' 
#' Naming format for the output:
#' Node set names: \code{{week}_F{differential regulation status}_M{differential regulation status}}.
#' Example 1: \code{1w_F1_M0} means up-regulation in females week 1 and null (zero effect) in males in week 1.
#' Example 2: \code{1w_F1_M-1} means up-regulation in females week 1 and down-regulation in males in week 1.
#' The edge set object is a named list of string vectors. The name of an edge is \code{[node_id]---[node_id]}.
#' 
#' @return A named list with three items: 
#' \describe{
#'   \item{\code{node_sets}}{named list of node sets}
#'   \item{\code{edge_sets}}{named list of edge sets}
#'   \item{\code{repfdr_results}}{repfdr EM solution}
#' }
#' 
#' @examples 
#' ### Example 1: Simulate data with a single cluster
#' zcolnames = c(
#'   paste("female",c("1w","2w","4w","8w"),sep="_"),
#'   paste("male",c("1w","2w","4w","8w"),sep="_")
#' )
#' zscores = matrix(rnorm(80000), ncol=8, dimnames = list(1:10000,zcolnames))
#' # now add a cluster with a strong signal 
#' zscores[1:500,1:4] = zscores[1:500,1:4] + 5 
#' 
#' # Run the clustering solution wrapper.
#' # When the data are "clean" (e.g., a mixture of two gaussians), we do not
#' # need a high df in the two-groups model estimation
#' # (default is 20, consider at least 10 when analyzing real data)
#' clustering_sol = bayesian_graphical_clustering(zscores, df=5)
#' # check if the clustering solution correctly assigns the first 500 rows 
#' # (with high prob) to the right nodes
#' length(intersect(1:500,clustering_sol$node_sets$`1w_F1_M0`))/500 > 0.95
#' length(intersect(1:500,clustering_sol$node_sets$`2w_F1_M0`))/500 > 0.95
#' length(intersect(1:500,clustering_sol$node_sets$`4w_F1_M0`))/500 > 0.95
#' length(intersect(1:500,clustering_sol$node_sets$`8w_F1_M0`))/500 > 0.95
#' 
#' # examine the node set sizes
#' sapply(clustering_sol$node_sets,length)
#' # examine the edge set sizes
#' sapply(clustering_sol$edge_sets,length)
#' 
#' # extract the top full trajectories in the data
#' # these should be the clusters with at least 10 features
#' min_cluster_size = 10
#' get_trajectory_sizes_from_edge_sets(clustering_sol$edge_sets,min_size = min_cluster_size)
#' 
#' \dontrun{
#' ### Example 2: real data
#' data(REPFDR_INPUTS, package="MotrpacRatTraining6moData")
#' zscores = REPFDR_INPUTS$zs_smoothed
#' rat_data_clustering_sol = bayesian_graphical_clustering(zscores)
#' # extract the largest trajectories
#' get_trajectory_sizes_from_edge_sets(rat_data_clustering_sol$edge_sets, min_size = 50)
#' # plot the top trajectories of the muscle tissues, color edges by tissue
#' get_tree_plot_for_tissue(tissues = c("SKM-GN","HEART","SKM-VL"), 
#'                          omes = "TRNSCRPT",
#'                          node_sets = rat_data_clustering_sol$node_sets,
#'                          edge_sets = rat_data_clustering_sol$edge_sets,
#'                          min_size = 20,
#'                          parallel_edges_by_tissue = TRUE,
#'                          max_trajectories = 3)
#' }
bayesian_graphical_clustering <- function(zscores, 
                                          min_analyte_posterior_thr=0.5,
                                          min_prior_for_config=0.001,
                                          df=20,
                                          naive_edge_sets=TRUE){
  
  if (!requireNamespace("repfdr", quietly = TRUE)){
    stop(
      "Package 'repfdr' must be installed to use this function.",
      call. = FALSE
    )
  }
  
  if(is.null(dim(zscores))||nrow(zscores)<5 || ncol(zscores)!=8){
    stop("Input zscore matrix does not fit the required input. See documentation for details.")
  }
  
  repfdr_results = repfdr_wrapper(zscores,min_prior_for_config=min_prior_for_config,df=df)
  repfdr_clusters = repfdr_results$repfdr_clusters
  repfdr_cluster_posteriors = repfdr_results$repfdr_cluster_posteriors
  
  node_sets = list()
  post_rowsums = rowSums(repfdr_cluster_posteriors,na.rm=TRUE)
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
    node_sets=node_sets,
    edge_sets=edge_sets,
    repfdr_results=repfdr_results
  ))
}


#' \code{repfdr} wrapper
#' 
#' A general wrapper for running \code{repfdr} on a matrix of z-scores.
#' 
#' @export
#' 
#' @param zscores A numeric matrix. Rows are analytes, columns are conditions (e.g., male week 1). Entries are z-scores.
#' @param min_prior_for_config A number. The minimal prior probability on the configuration priors.
#' @param df Integer. Degrees of freedom for fitting the estimated density f(z)
#' 
#' @return A named list with two objects:
#' \describe{
#'   \item{\code{repfdr_clusters}}{a binary matrix representation of the clusters}
#'   \item{\code{repfdr_cluster_posteriors}}{a matrix with the cluster posteriors}
#' }
#' 
#' @details 
#' This wrapper runs inference for the two groups model in each column of zscores, and then run repfdr's EM process.
#' To extract the fuzzy clustering solution, we exclude configurations whose prior probability is lower than min_prior_for_config.
#' 
#' @examples 
#' # Simulate data with a single cluster
#' zcolnames = c(
#'   paste("male",c("1w","2w","4w","8w"),sep="_"),
#'   paste("female",c("1w","2w","4w","8w"),sep="_")
#' )
#' zscores = matrix(rnorm(80000),ncol=8,dimnames = list(1:10000,zcolnames))
#' repfdr_results = repfdr_wrapper(zscores, df=5)
#' # in this example all configurations are null, 
#' # thus the  posteriors of the null cluster (all zeroes) are very high:
#' quantile(repfdr_results$repfdr_cluster_posteriors[,"00000000"])
#' # now add a cluster with a strong signal and rerun
#' zscores[1:500,1:4] = zscores[1:500,1:4] + 5
#' # When the data are "clean" (e.g., a mixture of two gaussians), we do not
#' # need a high df in the two-groups model estimation
#' # (default is 20, consider at least 10 when analyzing real data)
#' repfdr_results = repfdr_wrapper(zscores, df=5)
#' # look at the null cluster after adding the signal above
#' quantile(repfdr_results$repfdr_cluster_posteriors[,"00000000"],probs=c(0.05,0.1,0.5))
#' # now the posteriors of the first 500 rows, 
#' # with respect to the "planted" cluster should have high posteriors:
#' quantile(repfdr_results$repfdr_cluster_posteriors[1:500,"11110000"])
repfdr_wrapper <- function(zscores, min_prior_for_config = 0.001, df = 20){
  
  if (!requireNamespace("repfdr", quietly = TRUE)){
    stop(
      "Package 'repfdr' must be installed to use this function.",
      call. = FALSE
    )
  }
  
  # Step 1 in repfdr: discretize the z-scores and estimate the probabilities 
  # in each bin.
  # In our work we examined the diagnostic plots manually. The results look reasonable,
  # especially for the qqplots, so we decided to use the current estimation with df=20
  nbins = round(min(150,sqrt(nrow(zscores))-1))
  # In some cases, if the df is too high then the EM will not run properly.
  # When this happens, the EM will complete after a single iteration and the 
  # priors will not change. In this loop we check if this has occurred and if so
  # then we lower the df until the EM runs.
  while(df > 2){
    if(nrow(zscores) > 20000){
      # use the paper config only in very large datasets
      ztobins_res = repfdr::ztobins(zscores,
                                    df=df,
                                    type=1,
                                    n.bins=nbins,
                                    central.prop = 0.25,
                                    plot.diagnostics = FALSE,
                                    force.bin.number = TRUE)
    }
    else{
      # for med to small data sizes use the default with increased df
      ztobins_res = repfdr::ztobins(zscores,df=df,n.bins=nbins)
    }
    
    # Step 2 in repfdr: estimate the repfdr model using the EM algorithm
    repfdr_res = repfdr::repfdr(ztobins_res$pdf.binned.z,
                                ztobins_res$binned.z.mat,
                                non.null = 'replication',
                                control = repfdr::em.control(max.iter = 500,tol=1e-06))
    
    # check if the EM process ran
    curr_Pi_non_null = repfdr_res$Pi[-1,]
    num_unique_non_null_priors = length(unique(curr_Pi_non_null[,ncol(curr_Pi_non_null)]))
    if(num_unique_non_null_priors > 3){break}
    df = df-1
  }
  
  # sanity check
  configs = repfdr::hconfigs(ncol(zscores))
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
  repfdr_posteriors = repfdr::ldr(
    ztobins_res$pdf.binned.z,ztobins_res$binned.z.mat,
    repfdr_res$Pi,h.vecs = hvec_inds
  )
  # transform the output to objects that will be used below for 
  # the graphical analysis
  # this object keeps the configurations h in a matrix
  repfdr_clusters = t(repfdr_posteriors[,1:ncol(zscores)])
  if(length(hvec_inds)<2){
    repfdr_cluster_posteriors = rep(1,nrow(zscores))
    repfdr_cluster_posteriors = matrix(repfdr_cluster_posteriors,ncol=1)
    colnames(repfdr_cluster_posteriors) = names(hvec_inds)
    rownames(repfdr_cluster_posteriors) = rownames(zscores)
    return(list(repfdr_clusters=repfdr_clusters,
                repfdr_cluster_posteriors=repfdr_cluster_posteriors))
  }
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


#' Reduce a list of sets by a regex
#' 
#' An auxiliary function useful for filtering differential analyte sets by tissues or omes.
#' 
#' @param sets A named list of character vectors.
#' @param regs A character vector of regular expressions
#' @param append_semicol A logical. If TRUE (the default): append ';' as the suffix of each regex.
#' 
#' @details 
#' This function was added for helping with managing the output of the Bayesian clustering in real data.
#' It takes a list of analyte (e.g., gene ids) sets and removes items that do not match the regular expressions in regs.
#' 
#' @export
#' 
#' @return named list of filtered sets
#' 
#' @examples
#' sets = list(
#'   "cluster1" = c("muscle;g1","heart;g1","muscle;g2"),
#'   "cluster2" = c("muscle;g11","heart;g11","muscle;g12")
#' )
#' # remove non muscle analytes from the clustering solution above:
#' limit_sets_by_regex(sets,"muscle")
limit_sets_by_regex <- function(sets, regs, append_semicol = TRUE){
  if(is.null(regs) || length(regs)==0){return(sets)}
  l = list()
  for(nn in names(sets)){
    v = sets[[nn]]
    newv = c()
    for(r in regs){
      if(append_semicol){r = paste0(r,";")}
      newv = union(newv,v[grepl(r,v)])
    }
    l[[nn]] = newv
  }
  return(l)
}

#' Get trajectory sizes
#' 
#' Auxiliary function for getting the largest paths in a graphical solution.
#' 
#' @details
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
#' See bayesian_graphical_clustering for more details about the graphical analysis.
#' 
#' @param edge_sets A named list of string vectors. The name of an edge is \code{[node_id]---[node_id]}
#'        edges with no analytes have a NULL set (a set of size zero, but are still represented),
#'        node ids are \code{[timepoint]_F[x]_M[y]} where \code{x} and \code{y} represent the up/down state in each sex.
#'        \code{\link[MotrpacRatTraining6moData]{GRAPH_COMPONENTS}$edge_sets} by default.
#' @param min_size An integer specifying the minimal path size to be considered.
#' 
#' @return NULL if no paths of size of at least \code{min_size} were found, otherwise
#'         return a data frame that represents all paths of size min_size or greater,
#'         ranked from the largest path to the smallest one.
#'         
#' @export
#' 
#' @examples 
#' \dontrun{
#' ### Example 1: Simulate data with a single cluster
#' zcolnames = c(
#'   paste("female",c("1w","2w","4w","8w"),sep="_"),
#'   paste("male",c("1w","2w","4w","8w"),sep="_")
#' )
#' zscores = matrix(rnorm(80000),ncol=8,dimnames = list(1:10000,zcolnames))
#' # now add a cluster with a strong signal and rerun
#' zscores[1:500,1:4] = zscores[1:500,1:4] + 5
#' 
#' # run the clustering solution wrapper
#' clustering_sol = bayesian_graphical_clustering(zscores)
#' 
#' # extract the top full trajectories in the data
#' # these should be the clusters with at least 10 features
#' min_cluster_size=10
#' get_trajectory_sizes_from_edge_sets(clustering_sol$edge_sets, 
#'                                     min_size = min_cluster_size)
#' 
#' # extract the edges of the top two full trajectories
#' # this step "cleans" the edge sets by removing edges of trajectories with very few features
#' top2traj_edge_sets = filter_edge_sets_by_trajectories(
#'   clustering_sol$edge_sets, 
#'   topk = 2, 
#'   min_path_size = 10
#' )
#' # examine the new edge set sizes, excluded edges should have zero size
#' sapply(top2traj_edge_sets,length)
#' }
#' 
#' ### Example 2: Use published data 
#' get_trajectory_sizes_from_edge_sets()
get_trajectory_sizes_from_edge_sets <- function(edge_sets = MotrpacRatTraining6moData::GRAPH_COMPONENTS$edge_sets, 
                                                min_size = 10){
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
  trajectories = data.frame(trajectories,stringsAsFactors = FALSE)
  trajectories[[6]] = as.numeric(trajectories[[6]])
  trajectories = trajectories[order(trajectories[[6]],decreasing = TRUE),]
  return(trajectories)
}


#' Filter edge sets to largest trajectories
#' 
#' Keep the edges of the top trajectories of an edge set of a graphical solution.
#' 
#' @param edge_sets A named list of string vectors. The name of an edge is node_id---node_id
#'        edges with no analytes have a NULL set (a set of size zero, but are still represented),
#'        node ids are time_points_Fx_My where x and y represent the up/down state in each sex.
#'        \code{\link[MotrpacRatTraining6moData]{GRAPH_COMPONENTS}$edge_sets} by default.
#' @param topk A number. The maximal number of full trajectories to include in the new solution.
#' @param min_path_size An integer specifying the minimal path size to be considered.
#' 
#' @export
#' 
#' @return 
#' A named list of edge sets. All possible edges in our 9x4 grid will appear in the solution.
#' Edges that are removed will have no features/analytes in their entry.
#' 
#' @examples
#' \dontrun{
#' ### Example 1: Simulate data with a single cluster
#' zcolnames = c(
#'   paste("female",c("1w","2w","4w","8w"),sep="_"),
#'   paste("male",c("1w","2w","4w","8w"),sep="_")
#' )
#' zscores = matrix(rnorm(80000),ncol=8,dimnames = list(1:10000,zcolnames))
#' # now add a cluster with a strong signal and rerun
#' zscores[1:500,1:4] = zscores[1:500,1:4] + 5
#' 
#' # run the clustering solution wrapper
#' clustering_sol = bayesian_graphical_clustering(zscores)
#' 
#' # extract the top full trajectories in the data
#' # these should be the clusters with at least 10 features
#' min_cluster_size=10
#' get_trajectory_sizes_from_edge_sets(clustering_sol$edge_sets, min_size = min_cluster_size)
#' 
#' # extract the edges of the top two full trjectories
#' # this step "cleans" the edge sets by removing edges of trajectories with very few features
#' top2traj_edge_sets = filter_edge_sets_by_trajectories(clustering_sol$edge_sets,
#'                                                       topk = 2,
#'                                                       min_path_size = 10)
#' # examine the new edge set sizes, excluded edges should have zero size
#' sapply(top2traj_edge_sets,length)
#' # for comparison examine the edge sets of the Bayesian clustering solution:
#' sapply(clustering_sol$edge_sets,length)
#' }
#' 
#' ### Example 2: Use published data
#' # Get edges corresponding to 5 largest trajectories in the liver
#' tissue_edge_sets = limit_sets_by_regex(MotrpacRatTraining6moData::GRAPH_COMPONENTS$edge_sets,
#'                                        "LIVER")
#' filter_edge_sets_by_trajectories(tissue_edge_sets)
filter_edge_sets_by_trajectories <- function(edge_sets = MotrpacRatTraining6moData::GRAPH_COMPONENTS$edge_sets, 
                                             topk = 5, 
                                             min_path_size = 5){
  
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
  
  e_copy = data.table::copy(edge_sets)
  for(e in names(e_copy)){
    if(!(e %in% edges_to_keep)){
      e_copy[[e]] = character(0) 
    }
  }
  return(e_copy)
}


#' Graph representation of feature trajectories 
#' 
#' The main function for obtaining a graphical (tree) representation of the differential
#' analysis results.
#' 
#' @param tissues A character vector where values are in [MotrpacRatTraining6moData::TISSUE_ABBREV]. 
#'   The set of tissues to take for the analysis. If NULL take all.
#' @param omes A character vector where values are in [MotrpacRatTraining6moData::ASSAY_ABBREV].
#'   The set of omes to take for the analysis. If NULL take all.
#' @param node_sets A named list with the node (state) sets of analytes/features, see details for analyte name convention.
#'   \code{\link[MotrpacRatTraining6moData]{GRAPH_COMPONENTS}$node_sets} by default.
#' @param edge_sets A named list with the edge (state) sets of analytes/features, see details for analyte name convention.
#'   \code{\link[MotrpacRatTraining6moData]{GRAPH_COMPONENTS}$edge_sets} by default.
#' @param min_size A numeric. The threshold on the set sizes to be considered.
#' @param parallel_edges_by_ome A logical. TRUE means that we want to added parallel edges for the different omes.
#' @param parallel_edges_by_tissue A logical. TRUE means that we want to added parallel edges for the different tissues.
#' @param edge_width_range A numeric vector of size 2, a parameter for ggraph
#' @param edge_alpha_range A numeric vector of size 2, a parameter for ggraph
#' @param color_nodes_by_states A logical. If TRUE, nodes are colored by states. Red for up-reg, blue for down-reg, green for a discrepancy between the sexes.
#' @param max_trajectories A numeric or NULL. If not NULL then it specifies the number of pathways to keep when looking into the edge sets after filtering by omes and tissues. If  parallel_edges_by_tissue = TRUE then take the top trajectories in each tissue.
#' @param highlight_subset A character string or NULL. If not NULL then it specifies the name of a node, edge, or path to highlight in the tree. 
#' @param curvature A number between 0 and 1, where larger values yield more curved edges. 
#' 
#' @return a graph object
#' 
#' @importFrom igraph graph_from_data_frame edge_attr V E
#' @importFrom ggraph create_layout ggraph geom_edge_fan0 scale_edge_color_manual geom_edge_arc scale_edge_colour_identity scale_edge_width scale_edge_alpha geom_node_point
#' @importFrom grid unit
#' 
#' @export
#' 
#' @details 
#' The function filters the input set to include analytes from the given tissues and omes (if tissues/omes are not null).
#' If parallel edges are requested then the relevant edge sizes are computed internally and are used within ggraph for the output plot.
#' Analyte names are in the \code{\link[MotrpacRatTraining6moData]{ASSAY_ABBREV};\link[MotrpacRatTraining6moData]{TISSUE_ABBREV};[feature_ID]} format.
#' Node set names are are in the \code{[{1,2,4,8}]w_F[{-1,0,1}]_M[{-1,0,1}]} format for example: \code{1w_F-1_M0}.
#' Edge set names are in the \code{[node_a]---[node_b]} format, e.g., \code{4w_F0_M0---8w_F0_M1}.
#' 
#' @seealso [bayesian_graphical_clustering()] for more details about the graphical analysis.
#' 
#' @examples
#' \dontrun{
#' ### Example 1: redo the analysis using the rat data differential analysis results (z-scores)
#' data(REPFDR_INPUTS, package="MotrpacRatTraining6moData")
#' zscores = REPFDR_INPUTS$zs_smoothed
#' rat_data_clustering_sol = bayesian_graphical_clustering(zscores)
#' # extract the largest trajectories
#' get_trajectory_sizes_from_edge_sets(rat_data_clustering_sol$edge_sets,min_size = 50)
#' # plot the top trajectories of the muscle tissues, color edges by tissue
#' get_tree_plot_for_tissue(tissues = c("SKM-GN","HEART","SKM-VL"),
#'                          omes = "TRNSCRPT",
#'                          node_sets = rat_data_clustering_sol$node_sets, 
#'                          edge_sets = rat_data_clustering_sol$edge_sets, 
#'                          min_size = 20,
#'                          parallel_edges_by_tissue = TRUE,
#'                          max_trajectories = 3)
#' }
#' 
#' ### Example 2: load the graphical solutions from MotrpacRatTraining6moData 
#' ### and plot without rerunning the algorithm
#' get_tree_plot_for_tissue(
#'   tissues=c("SKM-GN","HEART","SKM-VL"),
#'   omes="TRNSCRPT",
#'   min_size = 20,
#'   parallel_edges_by_tissue = TRUE,
#'   max_trajectories = 3
#' )
#'
get_tree_plot_for_tissue <- function(
  tissues,
  omes = NULL,
  node_sets = MotrpacRatTraining6moData::GRAPH_COMPONENTS$node_sets,
  edge_sets = MotrpacRatTraining6moData::GRAPH_COMPONENTS$edge_sets,
  min_size = 20,
  parallel_edges_by_ome = FALSE,
  parallel_edges_by_tissue = FALSE,
  edge_width_range = c(0,10),
  edge_alpha_range = c(0,1),
  color_nodes_by_states = TRUE,
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
  highlight_edge = FALSE
  highlight_node = FALSE
  highlight_path = FALSE
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
        highlight_edge = TRUE
      }else if(grepl("->",highlight_subset)){
        highlight_path = TRUE
      }else{
        highlight_node = TRUE
      }
    }
  }
  
  if(!is.null(highlight_subset) & any(parallel_edges_by_ome, parallel_edges_by_tissue, color_nodes_by_states)){
    warning("Setting 'parallel_edges_by_ome', 'parallel_edges_by_tissue', and 'color_nodes_by_states' to FALSE because 'highlight_subset' is TRUE.")
    parallel_edges_by_ome = FALSE
    parallel_edges_by_tissue = FALSE
    color_nodes_by_states = FALSE
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
    edge_colors = MotrpacRatTraining6moData::ASSAY_COLORS
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
    edge_colors = MotrpacRatTraining6moData::TISSUE_COLORS
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
    edge_colors = MotrpacRatTraining6moData::TISSUE_COLORS
  }
  
  # transform the initial node/edge info matrices to data frames
  # filter edge sizes by the specified user input
  d = as.data.frame(edge_info)
  for(j in 3:ncol(d)){
    d[[j]] = as.numeric(d[[j]])
    d[d[,j]<min_size,j] = 0 # filter by min edge size
  }
  d_nodes = data.frame(node_info,check.names = FALSE)
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
    if(all(igraph::E(d_g)$size == 0)){
      message("No non-0 edges. Skipping.")
      return()
    }
    igraph::E(d_g)$edge_size = edge_attr(d_g,"size")
    d_g_auto_layout <- ggraph::create_layout(d_g, layout = 'auto')
  }else{
    d_g = igraph::graph_from_data_frame(d)
    if(all(igraph::E(d_g)$size == 0)){
      message("No non-0 edges. Skipping.")
      return()
    }
    igraph::E(d_g)$edge_size = edge_attr(d_g,"size")
    d_g_auto_layout <- ggraph::create_layout(d_g, layout = 'auto')
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
  d_g_our_layout  = data.table::copy(d_g_auto_layout)
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
  igraph::V(d_g)$setsize = d_nodes[V(d_g)$name,"size"]
  igraph::V(d_g)$setsize[V(d_g)$name == "0w"] = stats::median(igraph::V(d_g)$setsize)
  igraph::V(d_g)$label = sapply(V(d_g)$name,
                        function(x){a=strsplit(x,split="w_")[[1]];a[length(a)]})
  igraph::V(d_g)$label = gsub("_","\n",igraph::V(d_g)$label)
  igraph::V(d_g)$col = "gray"
  igraph::V(d_g)$alt_col = "gray"
  igraph::V(d_g)$shape = 15
  for(j in 1:length(d_g_ordered_nodes)){
    n = d_g_ordered_nodes[j]
    igraph::V(d_g)$col[grepl(n,d_g_our_layout$name)] = d_g_ordered_cols[n]
    igraph::V(d_g)$alt_col[grepl(n,d_g_our_layout$name)] = d_g_ordered_cols_alt[n]
    igraph::V(d_g)$shape[grepl(n,d_g_our_layout$name)] = d_g_ordered_shapes[n]
  }
  
  # Add the set size as a field in the layout data frame
  # This will be used for controlling the node size correctly
  d_g_our_layout["size"] = igraph::V(d_g)$setsize[d_g_our_layout$.ggraph.orig_index]

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
    igraph::V(d_g)$col = "gray"
    igraph::V(d_g)$alt_col = "gray"
    # Make edge red
    ecols = rep("gray", nrow(d))
    ecols[which(rownames(d)==highlight_cluster)] = "red"
    igraph::E(d_g)$col = ecols
  }else if(highlight_node){
    igraph::E(d_g)$col = "gray"
    vcols = rep("gray", length(names(igraph::V(d_g))))
    vcols[which(names(igraph::V(d_g))==highlight_cluster)] = "red"
    igraph::V(d_g)$col = vcols
    igraph::V(d_g)$alt_col = vcols
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
    igraph::E(d_g)$col = ecols
    # set colors for nodes
    vcols = rep("gray", length(names(igraph::V(d_g))))
    vcols[names(igraph::V(d_g))%in%curr_nodes] = "red"
    igraph::V(d_g)$col = vcols
    igraph::V(d_g)$alt_col = vcols
  }
  
  if(parallel_edges_by_ome || parallel_edges_by_tissue){
    # Assign NA color to 0-weight edges
    igraph::E(d_g)$type[igraph::E(d_g)$edge_size == 0] = NA
    p = ggraph::ggraph(d_g,layout = d_g_our_layout,) +  
      ggraph::geom_edge_fan0(
        ggplot2::aes(width = E(d_g)$edge_size,
            alpha=E(d_g)$edge_size,
            colour = E(d_g)$type),
        lineend = "round",
        strength = curvature) +
      ggraph::scale_edge_color_manual(values=edge_colors, name="", limits=names(edge_colors)[names(edge_colors) %in% c(tissues, omes)]) +
      ggplot2::guides(edge_color=ggplot2::guide_legend(override.aes = list(edge_width=3)))
  }else if(!is.null(highlight_subset)){
    # Assign NA color to 0-weight edges
    igraph::E(d_g)$col[igraph::E(d_g)$edge_size == 0] = NA
    p = ggraph::ggraph(d_g,layout = d_g_our_layout,) +  
      ggraph::geom_edge_arc(ggplot2::aes(width = E(d_g)$edge_size,
                        alpha = igraph::E(d_g)$edge_size, 
                        color = igraph::E(d_g)$col),
                    lineend = "round",strength = curvature) +
      ggraph::scale_edge_colour_identity(guide = "none")
  }else{
    # Assign NA color to 0-weight edges
    igraph::E(d_g)$col = "black"
    igraph::E(d_g)$col[igraph::E(d_g)$edge_size == 0] = NA
    p = ggraph::ggraph(d_g,layout = d_g_our_layout,) +  
      #geom_edge_link(aes(width = E(d_g)$edge_size,alpha=E(d_g)$edge_size)) + 
      ggraph::geom_edge_arc(ggplot2::aes(width = igraph::E(d_g)$edge_size,
                        alpha = igraph::E(d_g)$edge_size,
                        color = igraph::E(d_g)$col),
                    lineend = "round",strength = curvature) +
      ggraph::scale_edge_color_identity(guide = "none")
  }
  
  p = p +
    ggraph::scale_edge_width(range=edge_width_range,name="Intersect size") + 
    ggraph::scale_edge_alpha(range=edge_alpha_range,name="Intersect size") 
  
  if(color_nodes_by_states){
    p = p + ggraph::geom_node_point(ggplot2::aes(size = size),alpha=1,
                            color=igraph::V(d_g)$col,shape=igraph::V(d_g)$shape)
  }else{
    p = p + ggraph::geom_node_point(ggplot2::aes(size = size),alpha=1,
                            color=igraph::V(d_g)$alt_col,shape=igraph::V(d_g)$shape)
  }
  p = p +
    ggplot2::scale_size(range = c(2,20),name="Number of analytes") +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill="white"),
          legend.key.size = grid::unit(0.4, 'cm'))
  
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
  
  return(p)
}


#' An auxiliary function for getting the top analyte sets for a set of tissues.
#' 
#' @param tissues A character vector. The names of the tissues (one or more) to be considered.
#' @param node_sets A named list with the node (state) sets of analytes/features.
#'   \code{\link[MotrpacRatTraining6moData]{GRAPH_COMPONENTS}$node_sets} by default.
#' @param edge_sets A named list with the edge (state) sets of analytes/features.
#'   \code{\link[MotrpacRatTraining6moData]{GRAPH_COMPONENTS}$edge_sets} by default.
#' @param min_size An integer specifying the minimal path size to be considered.
#' @param k An integer. How many node/edge/trajectory sets to extract. See details.
#' @param add_week8 A logical. TRUE (default): add all week 8 nodes to the node set.
#' @param omes A character vector (optional). The names of the omes (one or more) to be considered.
#' 
#' @return A named list. Names correspond to set names (can correspond to a node, edge, or a trajectory).
#' 
#' @export
#' 
#' @details 
#' This function is useful for extracting the largest sets for a specific set of tissues and omes. 
#' Thus, it is a useful step before running enrichment analysis on sets of analytes identified by the graphical clustering analysis.
#' 
#' By specifying k the user can control how many sets to include. The default is three, 
#' which means that the top largest node sets, edge sets, and full trajectories will be extracted (three each).
#' 
#' Naming format:
#' Analyte names are in the ome;tissue;feature_id format.
#' Node set names: {week}_F{differential regulation status}_M{differential regulation status}.
#' Example 1: 1w_F1_M0 means up-regulation in females week 1 and null (zero effect) in males in week 1.
#' Example 2: 1w_F1_M-1 means up-regulation in females week 1 and down-regulation in males in week 1.
#' The edge set object is a named list of string vectors. The name of an edge is node_id---node_id.
#' A full trajectory is annotated by node_id->node_id->node_id->node_id. For example: 1w_F-1_M-1->2w_F-1_M-1->4w_F-1_M-1->8w_F-1_M-1.
#' 
#' See bayesian_graphical_clustering for more details about the graphical analysis.
#' 
#' @examples
#' filtered_solution = extract_tissue_sets(
#'    tissues = c("HEART","SKM-VL","SKM-GN")
#' )
extract_tissue_sets<-function(tissues,
                              node_sets = MotrpacRatTraining6moData::GRAPH_COMPONENTS$node_sets,
                              edge_sets = MotrpacRatTraining6moData::GRAPH_COMPONENTS$edge_sets,
                              k=3,
                              min_size=20,
                              add_week8=TRUE,
                              omes=NULL){
  
  tissue_node_sets = limit_sets_by_regex(node_sets,tissues)
  tissue_node_sets = limit_sets_by_regex(tissue_node_sets,omes)
  tissue_edge_sets = limit_sets_by_regex(edge_sets,tissues)
  tissue_edge_sets = limit_sets_by_regex(tissue_edge_sets,omes)
  l = list()
  
  # node sets
  node_set_sizes = sapply(tissue_node_sets,length)
  selected_node_sets = names(sort(node_set_sizes,decreasing = TRUE))[1:(k+1)]
  if(add_week8){
    selected_node_sets = union(selected_node_sets,
                               names(tissue_node_sets)[grepl("8w",names(tissue_node_sets))])
  }
  selected_node_sets = setdiff(selected_node_sets,"0w")
  l[selected_node_sets] = tissue_node_sets[selected_node_sets]
  
  # edge sets
  edge_set_sizes = sapply(tissue_edge_sets,length)
  edge_set_sizes = edge_set_sizes[!grepl("0w",names(edge_set_sizes))]
  selected_edge_sets = names(sort(edge_set_sizes,decreasing = TRUE))[1:k]
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


#' Remove all non-empty trajectories 
#' 
#' @param node_sets A named list with the node (state) sets of analytes/features.
#'   \code{\link[MotrpacRatTraining6moData]{GRAPH_COMPONENTS}$node_sets} by default.
#' @param edge_sets A named list with the edge (state) sets of analytes/features.
#'   \code{\link[MotrpacRatTraining6moData]{GRAPH_COMPONENTS}$edge_sets} by default.
#' @param tissues string vector, optional. tissue subset. all tissues by default
#' @param omes string vector, optional. ome subset. all omes by default
#' 
#' @export
#' 
#' @return named list with one element per trajectories. members are features in the path 
#' 
#' @examples
#' # Get lists of features belonging to all trajectories in the liver
#' get_all_trajectories(tissues = "LIVER")
#'
get_all_trajectories = function(edge_sets = MotrpacRatTraining6moData::GRAPH_COMPONENTS$edge_sets, 
                                node_sets = MotrpacRatTraining6moData::GRAPH_COMPONENTS$node_sets, 
                                tissues = MotrpacRatTraining6moData::TISSUE_ABBREV,
                                omes = MotrpacRatTraining6moData::ASSAY_ABBREV){
  
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
                  this happens if edge sets are not simple intersections of node sets")
    }
    l[[paste(tissue_top_trajs[path_i,2:5],collapse="->")]] = curr_set
  }
  return(l)
}


#' Extract main graphical clusters
#' 
#' Return a data frame with features from the 2 largest nodes, 2 largest edges, 10 largest non-null paths,
#' and all 8-week nodes from the graphical representation of training-regulated
#' features in each tissue. This code replicates the graphical clusters for which
#' pathway enrichment was performed for the landscape manuscript. 
#' 
#' @return a data frame with 5 columns and one row per combination of 
#'   feature ID and cluster:
#' \describe{
#'   \item{\code{feature}}{`r feature()`}
#'   \item{\code{cluster}}{character, cluster label}
#'   \item{\code{ome}}{`r assay()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{feature_ID}}{`r feature_ID()`} 
#' }
#' 
#' @export
#'
#' @examples
#' cluster_df = extract_main_clusters()
#' 
#' @details 
#' Notes about cluster labels:  
#' * All clusters are prefixed with the tissue abbreviation and a colon, e.g. "SKM-GN:"  
#' * Nodes are defined by the time point and state in each sex, where state is 1 for up,  
#' 0 for null, and -1 for down. For example, "1w_F-1_M-1" is a node that characterizes molecules 
#' at the "1w" time point that are down-regulated in females ("F-1") and down-regulated in males ("M-1"). 
#' These three pieces of information (time point, female state, male state) are separated by underscores ("_")   
#' * Edges contain "---" and connect a pair of nodes  
#' * Paths contain "->" and connect four nodes  
#' 
extract_main_clusters = function(){
  edge_sets = MotrpacRatTraining6moData::GRAPH_COMPONENTS$edge_sets
  node_sets = MotrpacRatTraining6moData::GRAPH_COMPONENTS$node_sets
  tissues = MotrpacRatTraining6moData::TISSUE_ABBREV
  tissues = tissues[!tissues %in% c("OVARY","TESTES","VENACV")]
  tree_analysis_selected_sets = c()
  
  # Add top 2 nodes, edges, paths from each cluster 
  # Note that null nodes and edges are currently INCLUDED 
  # All 8-week nodes are also included 
  for(tissue in tissues){
    curr_sets = extract_tissue_sets(tissue,
                                    k = 2,
                                    min_size = 20,
                                    node_sets,
                                    edge_sets)
    if(length(curr_sets)==0) next 
    names(curr_sets) = paste(tissue,names(curr_sets),sep=":")
    tree_analysis_selected_sets = c(tree_analysis_selected_sets,curr_sets)
  }
  
  # Add top 10 paths from each cluster 
  tree_analysis_selected_paths = list()
  for(tissue in tissues){
    curr_sets = extract_tissue_sets(tissue,k=10,min_size = 20,node_sets,edge_sets)
    if(length(curr_sets)==0) next
    names(curr_sets) = paste(tissue,names(curr_sets),sep=":")
    tree_analysis_selected_paths = c(tree_analysis_selected_paths,curr_sets)
  }
  # Keep paths only 
  tree_analysis_selected_paths = tree_analysis_selected_paths[grepl("->", names(tree_analysis_selected_paths))]
  # Remove all null paths
  tree_analysis_selected_paths = tree_analysis_selected_paths[!grepl("1w_F0_M0->2w_F0_M0->4w_F0_M0->8w_F0_M0", names(tree_analysis_selected_paths))]
  # Add them only if not already in `tree_analysis_selected_sets`
  tree_analysis_selected_paths = tree_analysis_selected_paths[!names(tree_analysis_selected_paths) %in% names(tree_analysis_selected_sets)]
  
  selected_sets = c(tree_analysis_selected_sets, tree_analysis_selected_paths)
  selected_sets_df  = check_cluster_res_format(selected_sets)
  return(selected_sets_df)
}


#' Check clustering results format 
#' 
#' This function will generate warnings or errors if there are issues with the format of the input. 
#' Used internally in [extract_main_clusters()].
#'
#' @param cluster_res Either a data frame or a list of lists. 
#'   If a data frame, it needs at least two columns: "feature" and "cluster". 
#'   The "feature" column should be in the format 
#'   '[MotrpacRatTraining6moData::ASSAY_ABBREV];[MotrpacRatTraining6moData::TISSUE_ABBREV];feature_ID'. 
#'   If a list of lists, each sublist must be named with the cluster name (character string), 
#'   and the values must be features in the format 
#'   '[MotrpacRatTraining6moData::ASSAY_ABBREV];[MotrpacRatTraining6moData::TISSUE_ABBREV];feature_ID'. 
#'
#' @return a data frame with 5 columns:
#' \describe{
#'   \item{\code{feature}}{`r feature()`}
#'   \item{\code{cluster}}{character, cluster label}
#'   \item{\code{ome}}{`r assay()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{feature_ID}}{`r feature_ID()`} 
#' }
#' 
#' @seealso [extract_main_clusters()]
#' 
#' @export 
#' 
#' @examples 
#' # Note this example is redundant because check_cluster_res_format()
#' # is called within extract_main_clusters(), but it provides an 
#' # example input and output
#' cluster_res = extract_main_clusters()
#' cluster_res_checked = check_cluster_res_format(cluster_res)
#' 
check_cluster_res_format = function(cluster_res){
  
  if(is.list(cluster_res) & !is.data.frame(cluster_res) & !data.table::is.data.table(cluster_res)){
    # convert to data.table
    dtlist = list()
    for(c in names(cluster_res)){
      dtlist[[c]] = data.table::data.table(feature=unlist(cluster_res[[c]]), cluster=c)
    }
    cluster_res = data.table::rbindlist(dtlist)
  }
  
  cluster_res = as.data.frame(cluster_res, check.names=F, stringsAsFactors=F)
  
  ## check colnames
  errors = c()
  if(!"feature" %in% colnames(cluster_res)){
    errors = c(errors, "\n'feature' not found in colnames of clustering results. Please reformat.")
  }
  if(!"cluster" %in% colnames(cluster_res)){
    errors = c(errors, "\n'cluster' not found in colnames of clustering results. Please reformat.")
  }
  if(length(errors)>0){
    stop(errors)
  }
  cluster_res$feature = as.character(cluster_res$feature)
  
  ## check that "feature" is in the right format
  ## should be "OME;TISSUE;feature_ID"
  
  # first check that there are two semicolons in each cell 
  if(!all(grepl("*;*;*", cluster_res$feature))){
    stop("Are all features in the 'feature' column of the clustering results in the format 'OME;TISSUE;FEATURE_ID'? At least one 'feature' does not have 2 semicolons.")
  }
  
  # then check that there are only 3 ;-separated strings
  if(!all(unlist(lapply(strsplit(cluster_res$feature, ';'), length))==3)){
    stop("There are more than 3 ;-separated strings in the 'feature' column of the clustering results. Make sure the 'feature' column is in the required format 'OME;TISSUE;FEATURE_ID'.")
  }
  
  # then check that the strings are correct 
  cluster_res = annotate_cluster_res(cluster_res)
  
  errors = c()
  assay_abbr = MotrpacRatTraining6moData::ASSAY_ABBREV
  if(!all(cluster_res$ome %in% assay_abbr)){
    not_rec = unique(cluster_res$ome[!cluster_res$ome %in% unname(assay_abbr)])
    errors = c(errors, sprintf("\nUnrecognized ome labels: %s", paste(not_rec, collapse=', ')))
  }
  tissue_abbr = MotrpacRatTraining6moData::TISSUE_ABBREV
  if(!all(cluster_res$tissue %in% tissue_abbr)){
    not_rec = unique(cluster_res$tissue[!cluster_res$tissue %in% unname(tissue_abbr)])
    errors = c(errors, sprintf("\nUnrecognized tissue labels: %s", paste(not_rec, collapse=', ')))
  } 
  if(length(errors)>0){
    stop(errors)
  }
  
  # ## check that "cluster" is in the right format
  # if(length(unique(cluster_res$cluster)) > 50){
  #   warning("There are more than 50 unique values in the 'cluster' column. Is this right?")
  # }
  
  return(cluster_res)
  
}


#' Annotate clustering results
#' 
#' Add "ome", "tissue", and "feature_ID" columns to the input.
#' Run within [check_cluster_res_format()]. 
#'
#' @param cluster_res data frame with at least two columns: "feature" and "cluster". 
#'   The "feature" column must be in the format 
#'   "[MotrpacRatTraining6moData::ASSAY_ABBREV];[MotrpacRatTraining6moData::TISSUE_ABBREV];feature_ID".
#'
#' @return the input as a data frame with additional "ome", "tissue", and "feature_ID" columns
#'   and unfactorized "cluster" and "feature" columns
#'   
#' @seealso [check_cluster_res_format()]
#' 
annotate_cluster_res = function(cluster_res){
  cluster_res = as.data.frame(cluster_res, stringsAsFactors = F, check.names = F)
  cluster_res$cluster = as.character(cluster_res$cluster) # keep it from being a factor
  cluster_res$feature = as.character(cluster_res$feature) # keep it from being a factor
  cluster_res$ome = sapply(cluster_res$feature, function(x) unname(unlist(strsplit(x, ';')))[1])
  cluster_res$tissue = sapply(cluster_res$feature, function(x) unname(unlist(strsplit(x, ';')))[2])
  cluster_res$feature_ID = sapply(cluster_res$feature, function(x) unname(unlist(strsplit(x, ';')))[3])
  return(cluster_res)
}


#' Plot feature composition of clusters 
#' 
#' Plot number of features in each cluster and fractions of features that belong to each ome and tissue. 
#' 
#' @param cluster_res Either a data frame or a list of lists. 
#'   If a data frame, it needs at least two columns: "feature" and "cluster". 
#'   The "feature" column should be in the format 
#'   '[MotrpacRatTraining6moData::ASSAY_ABBREV];[MotrpacRatTraining6moData::TISSUE_ABBREV];feature_ID'. 
#'   If a list of lists, each sublist must be named with the cluster name (character string), 
#'   and the values must be features in the format 
#'   '[MotrpacRatTraining6moData::ASSAY_ABBREV];[MotrpacRatTraining6moData::TISSUE_ABBREV];feature_ID'. 
#' 
#' @export 
#' 
#' @return [ggplot2::ggplot()] object
#' 
#' @examples 
#' # Get top 10 largest paths, nodes, edges in gastrocnemius 
#' # Exclude additional 8-week nodes
#' clusters = extract_tissue_sets("SKM-GN", k=10, add_week8=FALSE)
#' # Select paths only 
#' clusters = clusters[grepl("->", names(clusters))]
#' 
#' # Plot distribution of features
#' plot_features_per_cluster(clusters)
#' 
plot_features_per_cluster = function(cluster_res){
  
  cluster_res = check_cluster_res_format(cluster_res)
  
  # type = total, tissue, ome
  c1 = as.data.frame(table(cluster_res$cluster), stringsAsFactors=F)
  colnames(c1) = c('cluster','N')
  c1$colour = 'none'
  c1$type = 'total'
  
  # tissue
  c2 = as.data.frame(table(cluster_res$tissue, cluster_res$cluster), stringsAsFactors=F)
  colnames(c2) = c('colour','cluster','N')
  c2 = merge(c2, c1[,c('cluster','N')], by='cluster')
  c2$N = c2$N.x/c2$N.y*100
  c2$type = 'tissue'
  
  # ome
  c3 = as.data.frame(table(cluster_res$ome, cluster_res$cluster), stringsAsFactors=F)
  colnames(c3) = c('colour','cluster','N')
  c3$type = 'ome'
  c3 = merge(c3, c1[,c('cluster','N')], by='cluster')
  c3$N = c3$N.x/c3$N.y*100
  c = rbindlist(list(c1,c2,c3), use.names=T, fill=T)
  c$type = factor(c$type, levels=c('total','ome','tissue'))
  
  tissues = unique(cluster_res$tissue)
  tissues = tissues[order(tissues)]
  omes = unique(cluster_res$ome)
  omes = omes[order(omes)]
  
  # always use the same assay colors
  assay_cols = MotrpacRatTraining6moData::ASSAY_COLORS
  assay_cols = assay_cols[omes]
  
  c$colour = factor(c$colour, levels=c(rev(omes), rev(tissues), 'none'))
  c$cluster = factor(as.character(c$cluster))
  
  # order by N
  c = c[order(c$N, decreasing=T),]
  order = unique(c$cluster[c$type=='total'])
  
  # remove "tissue" or "ome" facet if there is only one present  
  fill_vals = c(MotrpacRatTraining6moData::TISSUE_COLORS,
                MotrpacRatTraining6moData::ASSAY_COLORS,
                none='black')
  breaks = c(omes, 'none', tissues)
  labels = c(omes, '----------', tissues)
  if(length(unique(cluster_res$tissue))==1){
    c = c[type!="tissue"]
    fill_vals = c(assay_cols,none='black')
    breaks = c(omes)
    labels = c(omes)
  }
  if(length(unique(cluster_res$ome))==1){
    c = c[type!="ome"]
    fill_vals = c(MotrpacRatTraining6moData::TISSUE_COLORS,none='black')
    breaks = c(tissues)
    labels = c(tissues)
  }
  
  c[type=="total", type := "Total N features"]
  c[type=="tissue", type := "% of feat. by tissue"]
  c[type=="ome", type := "% of feat. by ome"]
  g = ggplot2::ggplot(c, aes(x=cluster, y=N, fill=colour)) +
    ggplot2::geom_bar(stat='identity') +
    ggplot2::theme_classic() +
    ggplot2::facet_wrap(~type, ncol=3, scales='free_x') +
    ggplot2::scale_fill_manual(values = fill_vals,
                               breaks = breaks,
                               labels = labels) +
    ggplot2::coord_flip() +
    ggplot2::labs(x='Cluster') +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank()) +
    ggplot2::scale_x_discrete(limits=rev(order))
  
  # # make nicer data
  # data = data.table::data.table(table(cluster_res$cluster, cluster_res$ome, cluster_res$tissue), stringsAsFactors = F)
  # colnames(data) = c('cluster','ome','tissue','N')
  
  return(g)
}
