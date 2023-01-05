#' Pathway enrichment for graphical clusters 
#' 
#' Wrapper for multi-tissue, multi-omic pathway enrichment of clustering or 
#' graphical results.
#' Pathway enrichment is performed using [gprofiler2::gost()] separately for 
#' each unique combination of tissue, assay/ome, and cluster. 
#' 
#' @param cluster_res Either a data frame or a list of lists. 
#'   If a data frame, it needs at least two columns: "feature" and "cluster". 
#'   The "feature" column should be in the format 
#'   '[MotrpacRatTraining6moData::ASSAY_ABBREV];[MotrpacRatTraining6moData::TISSUE_ABBREV];feature_ID'. 
#'   If a list of lists, each sublist must be named with the cluster name (character string), 
#'   and the values must be features in the format 
#'   '[MotrpacRatTraining6moData::ASSAY_ABBREV];[MotrpacRatTraining6moData::TISSUE_ABBREV];feature_ID'. 
#' @param databases character vector of g:Profiler pathway databases to query. 
#'   "KEGG" and "REAC" (REACTOME) by default. 
#'   Current options include: GO (GO:BP, GO:MF, GO:CC to select a particular GO branch), 
#'   KEGG, REAC, TF, MIRNA, CORUM, HP, HPA, WP. 
#'   See [gprofiler2 documentation](https://rdrr.io/cran/gprofiler2/man/gost.html) for an up-to-date list. 
#' @param feature_to_gene data frame, map between \code{intersection_id_type} and gene symbols. 
#'   Columns must include "feature_ID", "gene_symbol", "ensembl_gene", and "kegg_id".
#'   [MotrpacRatTraining6moData::FEATURE_TO_GENE] by default. 
#' @param gene_identifier_type character, column in \code{feature_to_gene} that 
#'   matches the gene identifier type in \code{universe}.
#'   "ensembl_gene" by default. 
#' @param universe list of lists of character vectors, 
#'   named first by assay (i.e., [MotrpacRatTraining6moData::ASSAY_ABBREV])
#'   and then by tissue (i.e., [MotrpacRatTraining6moData::TISSUE_ABBREV]). 
#'   Vectors provide the full set of gene symbols associated with features tested
#'   during differential analysis. For example, \code{universe$TRNSCRPT$LUNG}
#'   should be a character vector of expressed genes in the lung, where the type
#'   of gene identifier matches \code{gene_identifier_type}. 
#'   \code{[MotrpacRatTraining6moData::GENE_UNIVERSES]$ensembl_gene} by default. 
#' @param kegg_db_destination character, target directory for KEGG 
#'   database used for \code{FELLA} pathway enrichment with metabolites. 
#'   Creates database if it doesn't exist yet.
#' @param fella_method character, enrichment method for [FELLA::enrich()], 
#'   one of "hypergeom" or "diffusion", passed to [run_fella()] 
#' @param min_input_set_size integer, input must have this minimum number of 
#'   unique mappable gene IDs to attempt enrichment with [gprofiler2::gost()] 
#' @param min_pw_set_size integer, pathway must have at least this many members 
#'   to attempt enrichment with [gprofiler2::gost()] 
#' @param max_pw_set_size integer, pathway must have no more than this many members 
#'   to attempt enrichment with [gprofiler2::gost()]  
#' @param adjust_p boolean, whether to adjust nominal p-values for multiple testing (IHW by tissue)
#' @param num_cores optional integer, number of cores to register if parallel computing is desired 
#' @param logfile optional character, path to log of failed iterations
#' @param maxattempt integer, max number of consecutive null results from [gprofiler2::gost()] 
#'   before giving up
#' 
#' @return data frame with enrichment results, or NULL if no enrichment results were returned:
#' \describe{
#'   \item{\code{query}}{character, the name of the input query which by default is the order of query with the prefix "query_" (from [gprofiler2::gost()])}
#'   \item{\code{term_size}}{integer, number of genes that are annotated to the term (from [gprofiler2::gost()])}
#'   \item{\code{query_size}}{integer, number of genes that were included in the query (from [gprofiler2::gost()])}
#'   \item{\code{intersection_size}}{integer, the number of genes in the input query that are annotated to the corresponding term (from [gprofiler2::gost()])}
#'   \item{\code{precision}}{double, the proportion of genes in the input list that are annotated to the function, defined as \code{intersection_size/query_size} (from [gprofiler2::gost()])}
#'   \item{\code{recall}}{double, the proportion of functionally annotated genes that the query recovers, defined as \code{intersection_size/term_size} (from [gprofiler2::gost()])}
#'   \item{\code{term_id}}{character, unique term/pathway identifier (from [gprofiler2::gost()])}
#'   \item{\code{source}}{character, the abbreviation of the data source for the term/pathway (from [gprofiler2::gost()])}
#'   \item{\code{term_name}}{character, term/pathway name (from [gprofiler2::gost()])}
#'   \item{\code{effective_domain_size}}{integer, the total number of genes in the universe used for the hypergeometric test (from [gprofiler2::gost()])}
#'   \item{\code{source_order}}{integer, numeric order for the term within its data source (from [gprofiler2::gost()])}
#'   \item{\code{parents}}{list of term IDs that are hierarchically directly above the term. 
#'     For non-hierarchical data sources this points to an artificial root node (from [gprofiler2::gost()]).}
#'   \item{\code{evidence_codes}}{character, comma-separated evidence codes (from [gprofiler2::gost()])}
#'   \item{\code{intersection}}{character, input gene IDs that intersect with the term/pathway (from [gprofiler2::gost()])}
#'   \item{\code{gost_adj_p_value}}{double, improperly adjusted hypergeometric p-value 
#'     from [gprofiler2::gost()]. For reference only; should not be used to filter results 
#'     unless there was only a single ome/tissue/cluster combination in the input.}
#'   \item{\code{computed_p_value}}{double, nominal hypergeometric p-value, 
#'     computed from the [gprofiler2::gost()] output}
#'   \item{\code{cluster}}{character, cluster label}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{ome}}{`r ome()`}
#'   \item{\code{adj_p_value}}{double, adjusted nominal p-value 
#'     \code{computed_p_value} using IHW with tissue as a covariate}
#' }
#' 
#' @export
#' 
#' @details 
#' [FELLA::enrich()] is used for pathway enrichment of metabolites; [gprofiler2::gost()] is used for 
#' all other omes assuming features have been mapped to genes. 
#' 
#' This function was used to generate [MotrpacRatTraining6moData::GRAPH_PW_ENRICH]. 
#' 
#' @examples 
#' \dontrun{
#' # Use graphical clusters as an example
#' cluster_res = extract_main_clusters()
#' # Pick a single graphical cluster
#' # Gastrocnemius features up-regulated in both males and females at 8 weeks of training
#' cluster_res = cluster_res[cluster_res$cluster == "SKM-GN:8w_F1_M1",]
#' 
#' # Example 1: Run pathway enrichment for this cluster on a single core
#' pw_enrich = cluster_pathway_enrichment(cluster_res)
#' 
#' # Example 2: Run pathway enrichment for this cluster on 4 cores
#' pw_enrich = cluster_pathway_enrichment(cluster_res, num_cores = 4)
#' 
#' # Example 3: Same as above, but include metabolites. 
#' # Use FELLA's hypergeometric method for enrichment. 
#' pw_enrich = cluster_pathway_enrichment(cluster_res, 
#'                                        num_cores = 4,
#'                                        kegg_db_destination = "~/KEGGdb/test",
#'                                        fella_method = "hypergeom")
#' }
cluster_pathway_enrichment = function(cluster_res, 
                                      databases = c('REAC','KEGG'), 
                                      feature_to_gene = MotrpacRatTraining6moData::FEATURE_TO_GENE, 
                                      gene_identifier_type = "ensembl_gene",
                                      universe = MotrpacRatTraining6moData::GENE_UNIVERSES$ensembl_gene, 
                                      kegg_db_destination = NULL,
                                      fella_method = "hypergeom", 
                                      min_input_set_size = 1, 
                                      min_pw_set_size = 10, 
                                      max_pw_set_size = 200,
                                      adjust_p = TRUE,
                                      num_cores = NULL,
                                      logfile = "/dev/null",
                                      maxattempt = 50){
  
  if(!requireNamespace("foreach", quietly = TRUE)){
    stop(
      "Package 'foreach' must be installed to run 'cluster_pathway_enrichment()'.",
      call. = FALSE
    )
  }
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    stop(
      "Package 'gprofiler2' must be installed to perform pathway enrichment.",
      call. = FALSE
    )
  }
  
  kegg_databaseDir = kegg_db_destination
  
  feature_to_gene_map = data.table::data.table(feature_to_gene)
  
  # initialize log file 
  write(sprintf("#%s", as.character(date())), file = logfile, append = F)
  write("cluster,tissue,ome,n_features,reason", file = logfile, append = T)
  
  # load KEGG data if path exists (only need to do this once before querying)
  if(!is.null(kegg_databaseDir)){
    if(!requireNamespace("FELLA", quietly = TRUE) | !requireNamespace("KEGGREST", quietly = TRUE)){
      stop(
        "Packages 'FELLA' and 'KEGGREST' must be installed to perform pathway enrichment with metabolites.",
        call. = FALSE
      )
    }
    make_kegg_db(kegg_databaseDir)
    # load database
    fella.data = FELLA::loadKEGGdata(
      databaseDir = kegg_databaseDir, 
      internalDir = FALSE,
      loadMatrix = c("hypergeom","diffusion"))
  }
  
  cluster_res = check_cluster_res_format(cluster_res)
  
  if(!requireNamespace("IHW", quietly = TRUE)){
    if(length(unique(cluster_res$tissue)) > 1){
      stop(
        "Package 'IHW' must be installed to perform multiple testing adjustment when more than one tissue is present in the input.",
        call. = FALSE
      )
    }
  }
  
  if(is.null(kegg_databaseDir) & 'METAB' %in% cluster_res$ome){
    warning("It looks like you have metabolites in your clustering result, but you did not provide a path for the KEGG database ('kegg_db_destination'). Pathway enrichment will be skipped for metabolites.\n")
  }
  
  # subset feature_to_gene_map to speed up intersection 
  #message("Subsetting feature-to-gene map...")
  feature_to_gene_map = feature_to_gene_map[feature_ID %in% cluster_res$feature_ID]
  #message("Done.")
  
  # get enrichments per tissue per ome
  # add a progress bar 
  #count = 0
  #total = length(unique(cluster_res$cluster))*length(unique(cluster_res$tissue))*length(unique(cluster_res$ome))
  #pb = txtProgressBar(min = 0, max = total, style = 3)
  
  ## set things up for parallelization
  # set up cores
  if(is.null(num_cores)){
    message("Number of cores not defined with 'num_cores' argument. Setting the number of cores to 1. Explicity define 'num_cores' to take advantage of parallelization. If you aren't sure how many cores are available to you, try 'length(parallel::mcaffinity()).'")
    num_cores = 1
  }
  
  if(num_cores > 1){
    if(!requireNamespace("doParallel", quietly = TRUE)){
      stop(
        "Package 'doParallel' must be installed to run 'cluster_pathway_enrichment()' with more than one core.",
        call. = FALSE
      )
    }
    doParallel::registerDoParallel(num_cores)
    message(sprintf("Registered %s core(s).", num_cores))
  }
  
  # 1 row per task 
  iterations = unique(data.table::as.data.table(cluster_res[,c('cluster','tissue','ome')]))
  
  # there is not a simple way to have anything within foreach print on the master node :( 
  message(sprintf("Running enrichments for %s gene sets on %s core(s). This may take a while...", nrow(iterations), num_cores))
  
  # only use dopar if more than one core is registered
  `%myinfix%` = ifelse(num_cores > 1, foreach::`%dopar%`, foreach::`%do%`)
  enrich_list = foreach::foreach(i = 1:nrow(iterations)) %myinfix% {
    if(i %% maxattempt == 0) message(sprintf("%s out of %s...", i, nrow(iterations)))
    res_dt = run_single_enrichment(feature_to_gene_map,
                                   universe,
                                   cluster_res,
                                   gene_identifier_type,
                                   databases,
                                   iterations, i,
                                   min_input_set_size,
                                   fella.data,
                                   fella_method, 
                                   logfile,
                                   maxattempt)
    res_dt
  }
  message("Done.")
  
  
  enrich_res = data.table::rbindlist(enrich_list, fill=T)
  
  # remove tests outside of the parameters
  if(min_pw_set_size>0 | max_pw_set_size<Inf){
    if('term_size' %in% colnames(enrich_res)){
      enrich_res = enrich_res[(term_size<=max_pw_set_size & term_size>=min_pw_set_size) | is.na(term_size)]
    }
  }
  
  # check if there are no results
  if(nrow(enrich_res)==0){
    message("No enrichment results were returned. Are you sure you are using the right universe and 'feature' format in the input?")
    return()
  }
  
  # fill in missing pathway names
  if(!'term_name' %in% colnames(enrich_res)){
    enrich_res[,term_name := NA_character_]
  }
  # pw = unique(enrich_res[,.(term_name, term_id)])
  # id_to_name = pw$term_name
  # names(id_to_name) = pw$term_id
  # id_to_name1 = id_to_name[!is.na(id_to_name)]
  # id_to_name2 = id_to_name[is.na(id_to_name)] # these should only be KEGG
  # id_to_name2 = id_to_name2[grepl('KEGG:\\d+$', names(id_to_name2))]
  # 
  # # fill in the rest with KEGG API
  # keggid2keggname = as.list(KEGGREST::keggList("pathway"))
  # names(keggid2keggname) = gsub("[A-z:]","",names(keggid2keggname))
  # id_to_name3 = unname(lapply(names(id_to_name2), 
  #                             function(x) unname(unlist(keggid2keggname[substring(x, 6)]))))
  # names(id_to_name3) = names(id_to_name2)
  # id_to_name = c(id_to_name1, id_to_name3)
  # enrich_res[is.na(term_name), term_name := as.character(sapply(term_id,
  #                                                               function(x) id_to_name[[x]]))]
  # enrich_res = enrich_res[term_name == 'NULL', term_name := term_id]
  
  # fill in missing KEGG IDs (for METAB)
  if(any(is.na(enrich_res[,term_name]))){
    keggid2keggname = as.list(KEGGREST::keggList("pathway"))
    names(keggid2keggname) = gsub("[A-z:]","",names(keggid2keggname))
    names(keggid2keggname) = paste0("KEGG:", names(keggid2keggname))
    enrich_res[is.na(term_name), term_name := as.character(sapply(term_id, function(x) keggid2keggname[[x]]))]
  }
  
  # make variable for IHW
  enrich_res = enrich_res[!is.na(computed_p_value)]
  if(adjust_p){
    if(length(unique(enrich_res[,tissue])) == 1){
      enrich_res[,adj_p_value := stats::p.adjust(computed_p_value, method="BH")]
    }else{
      ihw_res = IHW::ihw(enrich_res[,computed_p_value],
                         factor(enrich_res[,tissue]),
                         alpha=0.05)
      enrich_res[,adj_p_value := IHW::adj_pvalues(ihw_res)]
    }
  }
  
  # convert to data.frame 
  enrich_res = as.data.frame(enrich_res, stringsAsFactors=F)
  enrich_res$significant = NULL
  
  return(enrich_res)
}


#' Run single iteration of pathway enrichment 
#' 
#' Run single iteration of pathway enrichment (FELLA for METAB or gprofiler for all other omes)
#' NOT intended to be run directly. Run within [cluster_pathway_enrichment()]. 
#' 
#' @param feature_to_gene_map data table, map between \code{intersection_id_type} and gene symbols. 
#'   Columns must include "feature_ID", "gene_symbol", "ensembl_gene", and "kegg_id".
#' @param universe character vector of genes in the universe or background
#' @param cluster_res data frame output by [check_cluster_res_format()]
#' @param gene_identifier_type column in \code{feature_to_gene_map} that 
#'   matches the gene identifier type in \code{universe}
#' @param databases character vector of g:Profiler pathway databases to query. 
#'   "KEGG" and "REAC" (REACTOME) by default. 
#'   Current options include: GO (GO:BP, GO:MF, GO:CC to select a particular GO branch), 
#'   KEGG, REAC, TF, MIRNA, CORUM, HP, HPA, WP. 
#'   See [gprofiler2 documentation](https://rdrr.io/cran/gprofiler2/man/gost.html) for an up-to-date list. 
#' @param iterations data frame passed from [cluster_pathway_enrichment()] 
#'   that defines each iteration of this function
#' @param i integer index passed from [cluster_pathway_enrichment()] that specifies 
#'   which row of \code{iterations} to use
#' @param min_input_set_size integer, input must have this minimum number of 
#'   unique mappable gene IDs to attempt enrichment with [gprofiler2::gost()] 
#' @param fella.data FELLA database returned by [FELLA::loadKEGGdata()] passed to [run_fella()]
#' @param fella_method character, enrichment method for [FELLA::enrich()], 
#'   one of "hypergeom" or "diffusion", passed to [run_fella()] 
#' @param logfile optional character, path to log of failed iterations
#' @param maxattempt integer, max number of consecutive null results from [gprofiler2::gost()] 
#'   before giving up
#'
#' @seealso [cluster_pathway_enrichment()], [run_fella()]
#' 
#' @return data table with enrichment results, or NULL if no enrichment results were returned 
#' 
#' @keywords internal
#' 
run_single_enrichment = function(feature_to_gene_map, 
                                 universe, 
                                 cluster_res, 
                                 gene_identifier_type, 
                                 databases, 
                                 iterations, i, 
                                 min_input_set_size,
                                 fella.data,
                                 fella_method, 
                                 logfile = "/dev/null",
                                 maxattempt = 50){
  
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    stop(
      "Package 'gprofiler2' must be installed to perform pathway enrichment.",
      call. = FALSE
    )
  }
  
  CLUSTER = iterations[i,cluster]
  OME = iterations[i, ome]
  TISSUE = iterations[i, tissue]
  
  # get all features 
  curr_features = unique(cluster_res[cluster_res$cluster==CLUSTER &
                                       cluster_res$tissue==TISSUE &
                                       cluster_res$ome==OME, 'feature_ID'])
  if(OME=='METAB'){
    curr_identifiers = unique(feature_to_gene_map[feature_ID %in% curr_features, kegg_id])
  }else{
    curr_identifiers = unique(feature_to_gene_map[feature_ID %in% curr_features, get(gene_identifier_type)])
  }
  curr_identifiers = curr_identifiers[!is.na(curr_identifiers)]
  
  # define cluster labels
  labels = list(cluster=CLUSTER,
                tissue=TISSUE,
                ome=OME)
  
  curr_universe = universe[[OME]][[TISSUE]]
  if(length(curr_universe)==0){
    # add to log file 
    write(sprintf("%s,%s,%s,%s,no universe", CLUSTER, TISSUE, OME, length(curr_identifiers)), file = logfile, append = TRUE)
    return() 
  } 
  
  # check if all identifiers are in universe
  if(! all(curr_identifiers %in% curr_universe)){
    message(paste0("\n",sprintf("%s out of %s identifiers in %s %s %s are not in the current universe. Removing these features.",
                                length(curr_identifiers[!curr_identifiers %in% curr_universe]),
                                length(curr_identifiers),
                                CLUSTER, OME, TISSUE)))
    print(curr_identifiers[!curr_identifiers %in% curr_universe])
    curr_identifiers = curr_identifiers[curr_identifiers %in% curr_universe]
  }
  
  if(length(curr_identifiers) == 0){
    # add to log file 
    write(sprintf("%s,%s,%s,%s,no input features", CLUSTER, TISSUE, OME, length(curr_identifiers)), file = logfile, append = TRUE)
    return()  
  }
  if(length(curr_identifiers) < min_input_set_size){
    # add to log file 
    write(sprintf("%s,%s,%s,%s,too few input features", CLUSTER, TISSUE, OME, length(curr_identifiers)), file = logfile, append = TRUE)
    return()  
  }
  
  if(OME=='METAB'){ # KEGG IDs
    res_dt = run_fella(curr_identifiers, 
                       curr_universe, 
                       fella.data, 
                       fella_method)
    if(is.null(res_dt)){
      # add to log file 
      write(sprintf("%s,%s,%s,%s,null FELLA results", CLUSTER, TISSUE, OME, length(curr_identifiers)), file = logfile, append = TRUE)
      return()
    } 
  }else{ # genes
    # run gost
    res = NULL
    attempt = 0
    while(is.null(res)){
      attempt = attempt + 1
      res = gprofiler2::gost(curr_identifiers,
                             organism='rnorvegicus',
                             significant=FALSE, # return all results
                             correction_method='fdr', # it will not return unadjusted p-values. calculate those yourself
                             sources=databases,
                             custom_bg=curr_universe,
                             domain_scope='custom',
                             evcodes=TRUE)
      if(attempt == maxattempt){
        message(sprintf("No result returned for %s %s %s after %s attempts. Returning NULL.",
                        TISSUE, CLUSTER, OME, maxattempt))
        # add to log file 
        write(sprintf("%s,%s,%s,%s,no non-null results after %s attempts", CLUSTER, TISSUE, OME, length(curr_identifiers), maxattempt), file = logfile, append = TRUE)
        return()
      }
    }
    
    res_dt = data.table::data.table(res$result)
    
    # calculate raw p-values
    res_dt[,gost_adj_p_value := p_value]
    res_dt[,p_value := NULL]
    res_dt[,computed_p_value := stats::phyper(intersection_size-1, term_size, effective_domain_size-term_size, query_size,lower.tail= FALSE)]
    # res_dt[,my_bh := p.adjust(p_value, method='fdr'), by=source]
    # plot(-log10(res_dt[,my_bh]), -log10(res_dt[,adj_p_value]))
  }
  
  # add cluster labels 
  for(l in 1:length(labels)){
    res_dt[,(names(labels)[l]) := labels[l]]
  }
  
  return(res_dt)
}


#' Metabolomics pathway enrichment 
#' 
#' Run \code{FELLA} for metabolomics enrichment in KEGG pathways. 
#' Return results in a format similar to [gprofiler2::gost()] results.  
#'
#' @param input character vector, input features (KEGG IDs)
#' @param background character vector, background features (KEGG IDs)
#' @param fella.data \code{FELLA} database returned by [FELLA::loadKEGGdata()]
#' @param method character, \code{FELLA} enrichment method, one of "hypergeom", "diffusion"
#' @param niter integer, number of iterations with which to estimate p-values by simulation. 
#'   Only applies for the diffusion method of enrichment. 
#' 
#' @export
#' 
#' @return data table of pathway enrichment results:
#' \describe{
#'   \item{\code{term_size}}{double, number of KEGG IDs that are annotated to the term}
#'   \item{\code{query_size}}{integer, number of KEGG IDs that were included in the query}
#'   \item{\code{intersection_size}}{double, the number of KEGG IDs in the input query 
#'     that are annotated to the corresponding term}
#'   \item{\code{term_id}}{character, unique term/pathway identifier}
#'   \item{\code{source}}{character, database source of term/pathway}
#'   \item{\code{computed_p_value}}{double, nominal enrichment p-value}
#'   \item{\code{kegg_id}}{character, KEGG ID for \code{term_id}} 
#' }
#' 
#' @examples 
#' \dontrun{
#' # Make KEGG database
#' kegg_db = "~/KEGGdb/test"
#' make_kegg_db(kegg_db)
#' 
#' # Get FELLA data 
#' fella.data = FELLA::loadKEGGdata(databaseDir = kegg_db, 
#'                                  internalDir = FALSE,
#'                                  loadMatrix = c("hypergeom","diffusion"))
#' 
#' # Get input features
#' # Metabolites in SKM-GN:8w_F1_M1 that map to a KEGG ID
#' input = c("C02918","C00195","C01967","C00016","C04438","C02294","C00003","C00006",
#'           "C00157","C00350","C04475","C00344","C06125","C00550","C00387","C04230",
#'           "C00073","C00864","C00670","C00836","C00319")
#' 
#' # Get universe/background list 
#' background = MotrpacRatTraining6moData::GENE_UNIVERSES$gene_symbol$METAB$`SKM-GN`
#' 
#' # Example 1: method "hypergeom"
#' res = run_fella(input, background, fella.data, method="hypergeom")
#' 
#' # Example 2: method "diffusion" (more powerful but slower and more difficult to interpret)
#' res = run_fella(input, background, fella.data, method="diffusion")
#' }
run_fella = function(input, background, fella.data, method="hypergeom", niter=1e5){
  if(!requireNamespace("FELLA", quietly = TRUE)){
    stop(
      "Package 'FELLA' must be installed to perform pathway enrichment with metabolites.",
      call. = FALSE
    )
  }
  if(method=="hypergeom"){
    out = tryCatch(
      {
        myAnalysis = FELLA::enrich(
          compounds = input,
          method = "hypergeom",
          data = fella.data,
          compoundsBackground = background,
          p.adjust = "none")
        
        # extract all results myself because generateResultsTable() is too limited
        pvalues = myAnalysis@hypergeom@pvalues
        path_hits = myAnalysis@hypergeom@pathhits
        path_background = myAnalysis@hypergeom@pathbackground
        res_dt = data.table::data.table(
          term_size=path_background,
          query_size=length(myAnalysis@userinput@metabolites),
          intersection_size=path_hits,
          term_id=names(pvalues),
          source='KEGG',
          computed_p_value=unname(pvalues)
          )
        # get term_id in same format as gost
        res_dt[,kegg_id := term_id]
        res_dt[,term_id := gsub('rno','KEGG:',term_id)]
      },
      error=function(cond) {
        return(NULL)
      }
    )    
  }else if(method=="diffusion"){
    out = tryCatch(
      {
        myAnalysis = FELLA::enrich(
          compounds = input,
          method = "diffusion",
          approx = "simulation",
          niter = niter,
          data = fella.data,
          compoundsBackground = background)
        
        res_dt = data.table::data.table(
          query_size = length(myAnalysis@userinput@metabolites),
          term_id = names(myAnalysis@diffusion@pscores),
          source = 'KEGG',
          computed_p_value=unname(myAnalysis@diffusion@pscores)
          )
        
        # get term_id in same format as gost
        res_dt[,kegg_id := term_id]
        res_dt[,term_id := gsub('rno','KEGG:',term_id)]
        
        # only keep pathway results, not compounds, modules, etc.
        res_dt = res_dt[grepl("rno", kegg_id)]
      },
      error=function(cond) {
        return(NULL)
      }
    )    
  }else{
    stop("'method' must be one of 'hypergeom', 'diffusion'.")
  }
  return(out)
}


#' Make KEGG database
#' 
#' Make KEGG database for \code{FELLA} pathway enrichment with metabolites. 
#' Run internally in [cluster_pathway_enrichment()].
#'
#' @param kegg_db_destination character, target directory for KEGG 
#'   database. Parent directories are created if they do not yet exist. 
#' @export
#' @seealso [cluster_pathway_enrichment()], [run_fella()]
#' @examples 
#' \dontrun{
#' make_kegg_db("KEGGDB/20220921")
#' }
make_kegg_db = function(kegg_db_destination){
  if(!requireNamespace("FELLA", quietly = TRUE)){
    stop(
      "Package 'FELLA' must be installed to perform pathway enrichment with metabolites.",
      call. = FALSE
    )
  }
  if(endsWith(kegg_db_destination, "/")){
    kegg_db_destination = gsub("/$","",kegg_db_destination)
  }
  if(!dir.exists(kegg_db_destination)){
    kegg_path = dirname(kegg_db_destination)
    if(!dir.exists(kegg_path)){
      dir.create(kegg_path, recursive = TRUE)
    }
    kegg.graph = FELLA::buildGraphFromKEGGREST(organism='rno')
    FELLA::buildDataFromGraph(
      keggdata.graph = kegg.graph,
      databaseDir = kegg_db_destination, # directory must not exist yet
      internalDir = FALSE,
      matrices = c("hypergeom","diffusion"),
      niter = 100)
    message(sprintf("Created KEGG database:\n %s", kegg_db_destination))
  }else{
    message(sprintf("KEGG database already exists:\n %s", kegg_db_destination))
  }
}


#' Custom pathway enrichment for graphical clusters 
#' 
#' Wrapper for multi-tissue, multi-omic pathway enrichment of clustering or 
#' graphical results with a user-supplied list of pathways.
#' Pathway enrichment is performed using a hypergeometric test separately for 
#' each unique combination of tissue, assay/ome, and cluster. 
#' 
#' @param cluster_res Either a data frame or a list of lists. 
#'   If a data frame, it needs at least two columns: "feature" and "cluster". 
#'   The "feature" column should be in the format 
#'   '[MotrpacRatTraining6moData::ASSAY_ABBREV];[MotrpacRatTraining6moData::TISSUE_ABBREV];feature_ID'. 
#'   If a list of lists, each sublist must be named with the cluster name (character string), 
#'   and the values must be features in the format 
#'   '[MotrpacRatTraining6moData::ASSAY_ABBREV];[MotrpacRatTraining6moData::TISSUE_ABBREV];feature_ID'. 
#' @param pathway_member_list named list of character vectors where names are pathway names 
#'   and values are pathway members. Pathway members must match values in the 
#'   \code{gene_identifer_type} column of \code{feature_to_gene}. 
#' @param source optional character string to define the source of \code{pathway_member_list} 
#' @param feature_to_gene data frame, map between \code{intersection_id_type} and gene symbols. 
#'   Columns must include "feature_ID", "gene_symbol", "ensembl_gene", and "kegg_id".
#'   [MotrpacRatTraining6moData::FEATURE_TO_GENE] by default. 
#' @param gene_identifier_type character, column in \code{feature_to_gene} that 
#'   matches the gene identifier type in \code{universe}.
#'   "gene_symbol" by default. 
#' @param universe list of lists of character vectors, 
#'   named first by assay (i.e., [MotrpacRatTraining6moData::ASSAY_ABBREV])
#'   and then by tissue (i.e., [MotrpacRatTraining6moData::TISSUE_ABBREV]). 
#'   Vectors provide the full set of gene symbols associated with features tested
#'   during differential analysis. For example, \code{universe$TRNSCRPT$LUNG}
#'   should be a character vector of expressed genes in the lung, where the type
#'   of gene identifier matches \code{gene_identifier_type}. 
#'   \code{[MotrpacRatTraining6moData::GENE_UNIVERSES]$gene_symbol} by default. 
#' @param add_ensembl_intersection bool, whether to add a \code{intersection_ensembl}
#'   column, which converts gene IDs in the intersection to Ensembl IDs
#' @param min_input_set_size integer, input must have this minimum number of 
#'   unique mappable gene IDs to attempt enrichment 
#' @param min_pw_set_size integer, pathway must have at least this many members 
#'   to attempt enrichment
#' @param max_pw_set_size integer, pathway must have no more than this many members 
#'   to attempt enrichment
#' @param adjust_p boolean, whether to adjust nominal p-values for multiple testing (IHW by tissue)
#' @param num_cores optional integer, number of cores to register if parallel computing is desired 
#' @param logfile optional character, path to log of failed iterations
#' 
#' @export 
#'
#' @seealso [pathway_hypergeom_test()], [cluster_pathway_enrichment()]
#' 
#' @return data frame with enrichment results, or NULL if no enrichment results were returned:
#' \describe{
#'   \item{\code{term_size}}{integer, number of genes that are annotated to the term}
#'   \item{\code{query_size}}{integer, number of genes that were included in the query}
#'   \item{\code{intersection_size}}{integer, the number of genes in the input query that are annotated to the corresponding term}
#'   \item{\code{term_id}}{character, unique term/pathway identifier}
#'   \item{\code{source}}{character, the abbreviation of the data source for the term/pathway}
#'   \item{\code{term_name}}{character, term/pathway name}
#'   \item{\code{effective_domain_size}}{integer, the total number of genes in the universe used for the hypergeometric test}
#'   \item{\code{intersection}}{character, input gene IDs that intersect with the term/pathway}
#'   \item{\code{computed_p_value}}{double, nominal hypergeometric p-value}
#'   \item{\code{cluster}}{character, cluster label}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{ome}}{`r ome()`}
#'   \item{\code{adj_p_value}}{double, adjusted nominal p-value 
#'     \code{computed_p_value} using IHW with tissue as a covariate}
#' }
#' 
#' @examples 
#' # Use graphical clusters as an example
#' cluster_res = extract_main_clusters()
#' # Pick a single graphical cluster
#' # Gastrocnemius features up-regulated in both males and females at 8 weeks of training
#' cluster_res = cluster_res[cluster_res$cluster == "SKM-GN:8w_F1_M1",]
#' 
#' # Make a toy pathway member list with human gene symbols
#' pathways = list("TCA cycle" = c('ACO2','CS','FH','MDH1','OGDH','PDHA1','PDHA2','SDHC','SUCLG1'))
#' 
#' # Convert human gene symbols to rat gene symbols 
#' data("RAT_TO_HUMAN_GENE", package = "MotrpacRatTraining6moData")
#' for (pw in names(pathways)){
#'   newmembers = c()
#'   for (m in pathways[[pw]]){
#'     # get rat symbol
#'     rat = RAT_TO_HUMAN_GENE$RAT_SYMBOL[RAT_TO_HUMAN_GENE$HUMAN_ORTHOLOG_SYMBOL == m][1]
#'     if(!is.na(rat)){
#'       newmembers = c(newmembers, rat)
#'     }
#'   }
#'   pathways[[pw]] = newmembers
#' }
#' 
#' # Perform pathway enrichment 
#' custom_cluster_pathway_enrichment(cluster_res, 
#'                                   pathway_member_list = pathways,
#'                                   min_pw_set_size = 1)
#' 
custom_cluster_pathway_enrichment = function(cluster_res, 
                                             pathway_member_list, 
                                             source = 'custom', 
                                             feature_to_gene = MotrpacRatTraining6moData::FEATURE_TO_GENE, 
                                             gene_identifier_type = "gene_symbol",
                                             universe = MotrpacRatTraining6moData::GENE_UNIVERSES$gene_symbol, 
                                             add_ensembl_intersection = TRUE, 
                                             min_input_set_size = 1, 
                                             min_pw_set_size = 10, 
                                             max_pw_set_size = 200,
                                             adjust_p = TRUE,
                                             num_cores = NULL,
                                             logfile = "/dev/null"){
  
  if(!requireNamespace("foreach", quietly = TRUE)){
    stop(
      "Package 'foreach' must be installed to run 'custom_cluster_pathway_enrichment()'.",
      call. = FALSE
    )
  }
  
  feature_to_gene_map = data.table::data.table(feature_to_gene)
  
  # initialize log file 
  write(sprintf("#%s", as.character(date())), file = logfile, append = F)
  write("cluster,tissue,ome,n_features,reason", file = logfile, append = T)
  
  cluster_res = check_cluster_res_format(cluster_res)
  
  if(!requireNamespace("IHW", quietly = TRUE)){
    if(length(unique(cluster_res$tissue)) > 1){
      stop(
        "Package 'IHW' must be installed to perform multiple testing adjustment when more than one tissue is present in the input.",
        call. = FALSE
      )
    }
  }
  
  # subset feature_to_gene_map to speed up intersection 
  #message("Subsetting feature-to-gene map...")
  feature_to_gene_map = feature_to_gene_map[feature_ID %in% cluster_res$feature_ID]
  #message("Done.")
  
  ## set things up for parallelization
  # set up cores
  if(is.null(num_cores)){
    message("Number of cores not defined with 'num_cores' argument. Setting the number of cores to 1. Explicity define 'num_cores' to take advantage of parallelization. If you aren't sure how many cores are available to you, try 'length(parallel::mcaffinity()).'")
    num_cores = 1
  }
  
  if(num_cores > 1){
    if(!requireNamespace("doParallel", quietly = TRUE)){
      stop(
        "Package 'doParallel' must be installed to run 'cluster_pathway_enrichment()' with more than one core.",
        call. = FALSE
      )
    }
    doParallel::registerDoParallel(num_cores)
    message(sprintf("Registered %s core(s).", num_cores))
  }
  
  # 1 row per task 
  iterations = unique(data.table::as.data.table(cluster_res[,c('cluster','tissue','ome')]))
  
  # there is not a simple way to have anything within foreach print on the master node :( 
  message(sprintf("Running enrichments for %s gene sets on %s core(s). This may take a while...", nrow(iterations), num_cores))
  
  # only use dopar if more than one core is registered
  `%myinfix%` = ifelse(num_cores > 1, foreach::`%dopar%`, foreach::`%do%`)
  enrich_list = foreach::foreach(i = 1:nrow(iterations)) %myinfix% {
    if(i %% 50 == 0) message(sprintf("%s out of %s...", i, nrow(iterations)))
    res_dt = pathway_hypergeom_test(feature_to_gene_map,
                                    universe,
                                    cluster_res,
                                    pathway_member_list, 
                                    source, 
                                    gene_identifier_type,
                                    iterations, i,
                                    min_input_set_size,
                                    logfile,
                                    add_ensembl_intersection)
    res_dt
  }
  message("Done.")
  
  enrich_res = data.table::rbindlist(enrich_list, fill=T)
  
  # remove tests outside of the parameters
  if(min_pw_set_size>0 | max_pw_set_size<Inf){
    enrich_res = enrich_res[term_size<=max_pw_set_size & term_size>=min_pw_set_size]
  }
  
  # make variable for IHW
  enrich_res = enrich_res[!is.na(computed_p_value)]
  
  # check if there are no results
  if(nrow(enrich_res)==0){
    message("No enrichment results were returned. Are you sure you are using the right universe and 'feature' format in the input?")
    return()
  }
  
  # make variable for IHW
  if(adjust_p){
    if(length(unique(enrich_res[,tissue])) == 1){
      enrich_res[,adj_p_value := stats::p.adjust(computed_p_value, method="BH")]
    }else{
      ihw_res = IHW::ihw(enrich_res[,computed_p_value],
                         factor(enrich_res[,tissue]),
                         alpha=0.05)
      enrich_res[,adj_p_value := IHW::adj_pvalues(ihw_res)]
    }
  }
  
  # convert to data.frame 
  enrich_res = as.data.frame(enrich_res, stringsAsFactors=F)
  
  return(enrich_res)
}


#' Custom pathway enrichment test
#' 
#' Worker function run by [custom_cluster_pathway_enrichment()]. 
#' Not intended to be run independently. 
#' 
#' @param feature_to_gene data frame, map between \code{intersection_id_type} and gene symbols. 
#'   Columns must include "feature_ID", "gene_symbol", "ensembl_gene", and "kegg_id".
#'   [MotrpacRatTraining6moData::FEATURE_TO_GENE] by default. 
#' @param universe list of lists of character vectors, 
#'   named first by assay (i.e., [MotrpacRatTraining6moData::ASSAY_ABBREV])
#'   and then by tissue (i.e., [MotrpacRatTraining6moData::TISSUE_ABBREV]). 
#'   Vectors provide the full set of gene symbols associated with features tested
#'   during differential analysis. For example, \code{universe$TRNSCRPT$LUNG}
#'   should be a character vector of expressed genes in the lung, where the type
#'   of gene identifier matches \code{gene_identifier_type}. 
#'   \code{[MotrpacRatTraining6moData::GENE_UNIVERSES]$gene_symbol} by default. 
#' @param cluster_res Either a data frame or a list of lists. 
#'   If a data frame, it needs at least two columns: "feature" and "cluster". 
#'   The "feature" column should be in the format 
#'   '[MotrpacRatTraining6moData::ASSAY_ABBREV];[MotrpacRatTraining6moData::TISSUE_ABBREV];feature_ID'. 
#'   If a list of lists, each sublist must be named with the cluster name (character string), 
#'   and the values must be features in the format 
#'   '[MotrpacRatTraining6moData::ASSAY_ABBREV];[MotrpacRatTraining6moData::TISSUE_ABBREV];feature_ID'. 
#' @param pathway_member_list named list of character vectors where names are pathway names 
#'   and values are pathway members. Pathway members must match values in the 
#'   \code{gene_identifer_type} column of \code{feature_to_gene}. 
#' @param source optional character string to define the source of \code{pathway_member_list} 
#' @param gene_identifier_type character, column in \code{feature_to_gene} that 
#'   matches the gene identifier type in \code{universe}.
#'   "gene_symbol" by default. 
#' @param iterations data frame passed from [custom_cluster_pathway_enrichment()] 
#'   that defines each iteration of this function
#' @param i integer index passed from [custom_cluster_pathway_enrichment()] that specifies 
#'   which row of \code{iterations} to use
#' @param min_input_set_size integer, input must have this minimum number of 
#'   unique mappable gene IDs to attempt enrichment  
#' @param logfile optional character, path to log of failed iterations
#' @param add_ensembl_intersection bool, whether to add a \code{intersection_ensembl}
#'   column, which converts gene IDs in the intersection to Ensembl IDs
#' 
#' @return data table with enrichment results, or NULL if no enrichment results were returned:
#' \describe{
#'   \item{\code{term_size}}{integer, number of genes that are annotated to the term}
#'   \item{\code{query_size}}{integer, number of genes that were included in the query}
#'   \item{\code{intersection_size}}{integer, the number of genes in the input query that are annotated to the corresponding term}
#'   \item{\code{term_id}}{character, unique term/pathway identifier}
#'   \item{\code{source}}{character, the abbreviation of the data source for the term/pathway}
#'   \item{\code{term_name}}{character, term/pathway name}
#'   \item{\code{effective_domain_size}}{integer, the total number of genes in the universe used for the hypergeometric test}
#'   \item{\code{intersection}}{character, input gene IDs that intersect with the term/pathway}
#'   \item{\code{computed_p_value}}{double, nominal hypergeometric p-value}
#'   \item{\code{cluster}}{character, cluster label}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{ome}}{`r ome()`}
#' }
#' 
#' @seealso [custom_cluster_pathway_enrichment()]
#' 
#' @keywords internal
#' 
pathway_hypergeom_test = function(feature_to_gene,
                                  universe,
                                  cluster_res,
                                  pathway_member_list, 
                                  source, 
                                  gene_identifier_type,
                                  iterations, i,
                                  min_input_set_size,
                                  logfile,
                                  add_ensembl_intersection){
  
  CLUSTER = iterations[i,cluster]
  OME = iterations[i, ome]
  TISSUE = iterations[i, tissue]
  feature_to_gene_map = feature_to_gene
  
  # get all features 
  curr_features = unique(cluster_res[cluster_res$cluster==CLUSTER &
                                       cluster_res$tissue==TISSUE &
                                       cluster_res$ome==OME, 'feature_ID'])
  if(OME=='METAB'){
    curr_identifiers = unique(feature_to_gene_map[feature_ID %in% curr_features, kegg_id])
  }else{
    curr_identifiers = unique(feature_to_gene_map[feature_ID %in% curr_features, get(gene_identifier_type)])
  }
  curr_identifiers = curr_identifiers[!is.na(curr_identifiers)]
  
  # define cluster labels
  labels = list(cluster=CLUSTER,
                tissue=TISSUE,
                ome=OME)
  
  curr_universe = universe[[OME]][[TISSUE]]
  if(length(curr_universe)==0){
    # add to log file 
    write(sprintf("%s,%s,%s,%s,no universe", CLUSTER, TISSUE, OME, length(curr_identifiers)), file = logfile, append = TRUE)
    return() 
  } 
  
  # check if all identifiers are in universe
  if(! all(curr_identifiers %in% curr_universe)){
    message(paste0("\n",sprintf("%s out of %s identifiers in %s %s %s are not in the current universe. Removing these features.",
                                length(curr_identifiers[!curr_identifiers %in% curr_universe]),
                                length(curr_identifiers),
                                CLUSTER, OME, TISSUE)))
    print(curr_identifiers[!curr_identifiers %in% curr_universe])
    curr_identifiers = curr_identifiers[curr_identifiers %in% curr_universe]
  }
  
  if(length(curr_identifiers) == 0){
    # add to log file 
    write(sprintf("%s,%s,%s,%s,no input features", CLUSTER, TISSUE, OME, length(curr_identifiers)), file = logfile, append = TRUE)
    return()  
  }
  if(length(curr_identifiers) < min_input_set_size){
    # add to log file 
    write(sprintf("%s,%s,%s,%s,too few input features", CLUSTER, TISSUE, OME, length(curr_identifiers)), file = logfile, append = TRUE)
    return()  
  }
  
  # check if feature identifiers match up 
  if(!any(unlist(pathway_member_list) %in% curr_universe)){
    # add to log file 
    write(sprintf("%s,%s,%s,%s,pathway members not found in universe", CLUSTER, TISSUE, OME, length(curr_identifiers)), file = logfile, append = TRUE)
    return()  
  }
  
  res_list = list()
  for(curr_pathway in names(pathway_member_list)){
    
    pw_members = pathway_member_list[[curr_pathway]]
    # subset pathway members down to the current universe
    pw_members = pw_members[pw_members %in% curr_universe]
    
    if(length(pw_members) == 0){
      next
    }
    
    intersection = intersect(pw_members, curr_identifiers)
    if(add_ensembl_intersection){
      # convert to ensembl
      intersection2 = unique(feature_to_gene_map[get(gene_identifier_type) %in% intersection, ensembl_gene])
      intersection2 = intersection2[!is.na(intersection2)]
    }
    
    query_size = length(curr_identifiers)
    term_size = length(pw_members)
    effective_domain_size = length(curr_universe)
    intersection_size = length(intersection)
    
    # run hypergeometric test
    hypergeom_p = stats::phyper(intersection_size-1, term_size, effective_domain_size-term_size, query_size,lower.tail= FALSE)
    
    curr_res_dt = data.table::data.table(term_size = term_size,
                                         query_size = query_size,
                                         intersection_size = intersection_size,
                                         term_id = gsub(" |-|,|\\.|;|:","_", curr_pathway),
                                         source = source,
                                         term_name = curr_pathway,
                                         effective_domain_size = effective_domain_size,
                                         intersection = paste0(intersection, collapse=','),
                                         computed_p_value = hypergeom_p)
    if(add_ensembl_intersection){
      curr_res_dt[,intersection_ensembl := paste0(intersection2, collapse=',')]
    }
    
    res_list[[curr_pathway]] = curr_res_dt
  }
  res_dt = data.table::rbindlist(res_list, fill=T)
  
  # add cluster labels 
  for(l in 1:length(labels)){
    res_dt[,(names(labels)[l]) := labels[l]]
  }
  
  # remove rows with 0 intersection, as gprofiler does
  res_dt = res_dt[intersection_size > 0]
  if(nrow(res_dt)==0){
    return()
  }
  
  return(res_dt)
}


# TODO
# low priority 
plot_top_enrichments_per_cluster = function(){}


#' Gene pathway enrichment
#'
#' Perform pathway enrichment for a list of genes. 
#' Wrapper for [gprofiler2::gost()]. 
#'
#' @param input vector of gene identifiers
#' @param background vector of gene identifiers in universe/background. Must include the 
#'   input. Should include all genes that were candidates for \code{input}, e.g., 
#'   expressed genes. If performing pathway enrichment on genes identified through \code{MotrpacRatTraining6moData},
#'   see [MotrpacRatTraining6moData::GENE_UNIVERSES]. See \code{custom_bg} argument for [gprofiler2::gost()]. 
#' @param organism character, species name. Default: "rnorvegicus" See \code{organism}
#'   argument for [gprofiler2::gost()]. 
#' @param databases character, vector of databases in which to search. See \code{sources}
#'   argument for [gprofiler2::gost()]. 
#' @param min_pw_set_size integer, exclude pathways smaller than this. Default: 10
#' @param max_pw_set_size integer, exclude pathways larger than this. Default: 200
#' @param return_gem bool, whether to return a data frame compatible with the
#'   Generic Enrichment Map (GEM) file format, which can be used as an input for the 
#'   [Cytoscape EnrichmentMap application](https://apps.cytoscape.org/apps/enrichmentmap).
#'   Default: FALSE
#'
#' @return a data frame of pathway enrichment results. If \code{return_gem=TRUE},
#'   column names are changed add added to be compatible with the Generic Enrichment Map (GEM) file format.
#'   Otherwise, the columns are defined as follows:
#' \describe{
#'   \item{\code{query}}{character, the name of the input query which by default is the order of query with the prefix "query_"}
#'   \item{\code{significant}}{bool, whether the [gprofiler2::gost()] enrichment p-value is less than 0.05}
#'   \item{\code{term_size}}{integer, number of genes that are annotated to the term}
#'   \item{\code{query_size}}{integer, number of genes that were included in the query (from [gprofiler2::gost()])}
#'   \item{\code{intersection_size}}{integer, the number of genes in the input query that are annotated to the corresponding term}
#'   \item{\code{precision}}{double, the proportion of genes in the input list that are annotated to the function, defined as \code{intersection_size/query_size}}
#'   \item{\code{recall}}{double, the proportion of functionally annotated genes that the query recovers, defined as \code{intersection_size/term_size}}
#'   \item{\code{term_id}}{character, unique term/pathway identifier}
#'   \item{\code{source}}{character, the abbreviation of the data source for the term/pathway}
#'   \item{\code{term_name}}{character, term/pathway name}
#'   \item{\code{effective_domain_size}}{integer, the total number of genes in the universe used for the hypergeometric test}
#'   \item{\code{source_order}}{integer, numeric order for the term within its data source}
#'   \item{\code{parents}}{list of term IDs that are hierarchically directly above the term. 
#'     For non-hierarchical data sources this points to an artificial root node}
#'   \item{\code{evidence_codes}}{character, comma-separated evidence codes}
#'   \item{\code{intersection}}{character, input gene IDs that intersect with the term/pathway}
#'   \item{\code{gost_adj_p_value}}{double, improperly adjusted hypergeometric p-value 
#'     from [gprofiler2::gost()]. For reference only; should not be used to filter results 
#'     unless there was only a single list of genes in the input.}
#'   \item{\code{computed_p_value}}{double, nominal hypergeometric p-value, 
#'     computed from the [gprofiler2::gost()] output}
#'   \item{\code{BH_adj_p_value}}{double, BH-adjusted p-values, calculated on \code{computed_p_value}} 
#'}
#' @export
#'
#' @examples
#' # Perform pathway enrichment for differential transcripts in the liver
#' diff = MotrpacRatTraining6moData::TRAINING_REGULATED_FEATURES
#' input_feat = diff$feature_ID[diff$tissue == "LIVER" & diff$assay == "TRNSCRPT"]
#' map = MotrpacRatTraining6moData::FEATURE_TO_GENE_FILT
#' input = unique(map$gene_symbol[map$feature_ID %in% input_feat])
#' background = MotrpacRatTraining6moData::GENE_UNIVERSES$gene_symbol$TRNSCRPT$LIVER
#' res = gene_pathway_enrichment(input, background)
#' head(res)
gene_pathway_enrichment = function(input, 
                                   background, 
                                   organism = 'rnorvegicus', 
                                   databases = c('REAC','WP','KEGG'),
                                   min_pw_set_size = 10,
                                   max_pw_set_size = 200,
                                   return_gem = FALSE){
  
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    stop(
      "Package 'gprofiler2' must be installed to run 'gene_pathway_enrichment()'.",
      call. = FALSE
    )
  }
  
  res = gprofiler2::gost(input,
                         organism = organism,
                         significant = FALSE, # return all results
                         correction_method = 'fdr', # it will not return unadjusted p-values. calculate those yourself
                         sources = databases,
                         custom_bg = background,
                         domain_scope = 'custom',
                         evcodes = TRUE)
  
  res_dt = data.table::data.table(res$result)
  
  if(nrow(res_dt)==0){
    warning("No results returned by 'gprofiler2::gost()'.\n")
    return()
  }
  
  # calculate raw p-values
  res_dt[,gost_adj_p_value := p_value]
  res_dt[,p_value := NULL]
  res_dt[,computed_p_value := stats::phyper(intersection_size-1, term_size, effective_domain_size-term_size, query_size, lower.tail= FALSE)]
  
  # remove tests outside of the parameters
  if(min_pw_set_size>0 | max_pw_set_size<Inf){
    res_dt = res_dt[term_size<=max_pw_set_size & term_size>=min_pw_set_size]
  }
  
  # fix parents column 
  res_dt[,parents := sapply(parents, function(x) paste0(unlist(x), collapse=','))]
  
  res_dt[,BH_adj_p_value := stats::p.adjust(computed_p_value, method='BH')]
  res_dt = res_dt[order(BH_adj_p_value, decreasing=F)]
  
  if(return_gem){
    res_dt[,Phenotype := 1]
    res_dt = res_dt[,.(term_id, term_name, computed_p_value, BH_adj_p_value, Phenotype, intersection)]
    colnames(res_dt) = c('GO.ID','Description','p.Val','FDR','Phenotype','Genes')
  }
  
  return(as.data.frame(res_dt))
}

