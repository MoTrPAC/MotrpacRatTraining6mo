#' Prepare GCT file for GSEA
#' 
#' Take a set of differential analysis results, either for a specific tissue and ome or
#' from a user-supplied table, convert to GCT format, and write to file for GSEA.  
#'
#' @param tissue `r tissue()`. Must be specified if cusotm \code{input} is not provided.   
#' @param assay character, assay abbreviation, one of "PROT", "PHOSPHO", "TRNSCRPT", "ACETYL", "UBIQ" 
#' @param gene_id_type character, gene identifier type. Must match the gene ID type in the
#'   gene set database you plan on using for GSEA. 
#'   One of "gene_symbol", "entrez", "ensembl", "refseq". Default: "gene_symbol"
#' @param outdir character, output directory for GCT file. The directory is created if it does not already exist. 
#'   Current directory by default. 
#' @param outfile_prefix character, prefix for output GCT file. By default, this prefix
#'   includes the specified tissue and assay and current date. Must be specified for custom input data. 
#' @param input optional data frame if the user wants to perform this analysis for 
#'   a custom set of differential analysis results. 
#'   Required columns are "tscore", "feature_ID", and \code{cast_vars} 
#'   OR "tscore", \code{cast_vars}, and the gene identifier 
#'   indicated by \code{gene_id_type}. If a "feature_ID" column exists but not  
#'   a column corresponding to \code{gene_id_type}, then \code{feature_to_gene_map}
#'   must map between "feature_ID" and \code{gene_id_type}. 
#' @param cast_vars character vector of column names in the differential analysis results 
#'   that are used to convert the table from long to wide format, with t-scores as the value variable. 
#'   See [data.table::dcast()] for more details. Default: "sex", "comparison_group"  
#' @param feature_to_gene_map data frame, map between "feature_ID" and \code{gene_id_type}.
#'   [MotrpacRatTraining6moData::FEATURE_TO_GENE_FILT] if not otherwise specified. 
#'
#' @return character, path of the GCT file  
#' @export
#' @importFrom methods new
#'
#' @examples
#' # Input for GSEA where gene set databases use gene symbols as IDs.
#' # This applies for the default gene set database 
#' # when 'method="gsea"' in 'ssGSEA2_wrapper()'.
#' prepare_gsea_input("HEART","PROT",outdir="/tmp")
#' prepare_gsea_input("LIVER","PHOSPHO",outdir="/tmp")
#' prepare_gsea_input("LIVER","ACETYL",outdir="/tmp")
#' prepare_gsea_input("LIVER","UBIQ",outdir="/tmp")
#' prepare_gsea_input("LIVER","TRNSCRPT",outdir="/tmp")
#' 
#' # Input for GSEA with Mitocarta (i.e., 'method="gsea_mitocarta"' 
#' # in 'ssGSEA2_wrapper()'), which uses RefSeq IDs 
#' prepare_gsea_input("LIVER","TRNSCRPT",outdir="/tmp",gene_id_type="refseq")
#' 
#' # "Custom" input
#' res = combine_da_results(tissues = "KIDNEY", assays = "PROT")
#' # add dummy column
#' res$gene_symbol = res$feature_ID
#' prepare_gsea_input(input=res, outdir="/tmp", outfile_prefix="KIDNEY_PROT")
#' 
#' 
#' @details 
#' T-scores from the timewise differential analysis results are used for scores. 
#' Feature-level data is summarized into gene-level data using the maximum absolute t-score. 
#' 
prepare_gsea_input = function(tissue = NULL, 
                              assay = NULL, 
                              gene_id_type = "gene_symbol", 
                              outdir = ".", 
                              outfile_prefix = NULL, 
                              input = NULL, 
                              cast_vars = c("sex","comparison_group"),
                              feature_to_gene_map = NULL){
  
  if(!requireNamespace("cmapR", quietly = TRUE)){
    stop(
      "Package 'cmapR' must be installed to run 'da_results_to_gct()'.",
      call. = FALSE
    )
  }
  
  if(is.null(input) & any(is.null(c(tissue, assay)))){
    stop("If 'input' is not specified, both 'tissue' and 'assay' must be specified.")
  }

  if(!is.null(feature_to_gene_map)){
    feature_map = data.table::as.data.table(feature_to_gene_map)
    # check colnames
    if(!all(c("feature_ID","gene_id_type") %in% colnames(feature_map))){
      stop(sprintf("When 'feature_to_gene_map' is specified, columns must include 'feature_ID' and '%s'.", gene_id_type))
    }
  }
  
  # define outfile 
  date = gsub("-","",Sys.Date())
  if(is.null(outfile_prefix)){
    if(!is.null(input)){
      stop("If a custom input is provided, 'outfile_prefix' must be defined.")
    }
    outfile_prefix = sprintf("MotrpacRatTraining6mo_gsea_%s_%s_%s", tissue, assay, date)
  }
  if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
  }
  
  # load differential analysis results
  merged = NULL
  if(!is.null(input)){
    da = data.table::as.data.table(input)
    # check for required columns 
    req_cols = c("tscore", gene_id_type, cast_vars)
    missing = c()
    for(c in req_cols){
      if(!c %in% colnames(da)){
        missing = c(missing, c)
      }
    }
    if(length(missing)>0){
      stop(sprintf("The following required columns are missing from the custom input:\n%s",
                   paste0(missing, collapse=", ")))
    }
    merged = da
  }else{
    # get feature metadata
    if(!assay %in% c("PROT","PHOSPHO","TRNSCRPT","ACETYL","UBIQ")){
      stop(sprintf("GSEA/PTM-SEA not compatible with %s differential analysis results.", assay))
    }
    da = data.table::as.data.table(combine_da_results(tissues = tissue, assays = assay))
  }
  
  # handle gene id type 
  if(!gene_id_type %in% colnames(da)){
    if(is.null(feature_to_gene_map)){
      feature_map = data.table::as.data.table(MotrpacRatTraining6moData::FEATURE_TO_GENE_FILT)
      if(gene_id_type == "entrez"){
        feature_map[,entrez := entrez_gene]
      }else if(gene_id_type == "ensembl"){
        feature_map[,ensembl := ensembl_gene]
      }else if(gene_id_type == "refseq"){
        # there isn't a refseq column, but protein feature_ID are RefSeq IDs, so take advantage of those 
        refseq_map = unique(feature_map[grepl("^[A-Z]P_[0-9]", feature_ID) & !is.na(gene_symbol), .(gene_symbol, feature_ID)])
        setnames(refseq_map,"feature_ID","refseq")
        refseq_map[,refseq := gsub("_[A-Z].*","",refseq)]
        refseq_map = unique(refseq_map)
        feature_map = unique(feature_map[!is.na(gene_symbol) & feature_ID %in% da[,feature_ID],.(feature_ID, gene_symbol)])
        feature_map = merge(feature_map, refseq_map, by="gene_symbol")
      }
      feature_map = unique(feature_map[!is.na(get(gene_id_type)), c("feature_ID",gene_id_type), with=FALSE])
    }
  }
  
  if(is.null(merged)){
    # merge with differential analysis results
    merged = merge(da, feature_map, by="feature_ID", allow.cartesian=TRUE)
    merged[,feature_ID := NULL]
    merged = unique(merged)
  }
  
  # handle column names 
  if(!"tscore" %in% colnames(merged)){
    merged[,tscore := NA_real_]
  }
  if(!"zscore" %in% colnames(merged)){
    merged[,zscore := NA_real_]
  }
  
  # make table of scores 
  merged[,gene := get(gene_id_type)]
  merged[is.na(tscore), tscore := zscore]
  merged = merged[!is.na(tscore), .(gene, sex, comparison_group, tscore)]
  # get max t-score per gene ID
  merged_max = merged[!is.na(gene), list(max_tscore = max(tscore)),
                      by=.(gene, sex, comparison_group)]
  # cast wide
  form = sprintf("gene ~ %s", paste0(cast_vars, collapse=" + "))
  merged_wide = data.table::dcast(merged_max, formula = eval(parse(text=form)), value.var='max_tscore')
  
  # convert to matrix
  rn = merged_wide[,gene]
  merged_wide[,gene := NULL]
  mat = as.matrix(merged_wide)
  colnames(mat) = colnames(merged_wide)
  rownames(mat) = rn
  
  # make GCT
  gct = methods::new("GCT", mat = mat)
  
  # write to file s
  outfile = sprintf("%s/%s.gct", outdir, outfile_prefix)
  cmapR::write_gct(gct, ofile = outfile, appenddim = FALSE)
  path = list.files(path=dirname(outfile), pattern=basename(outfile), full.names=TRUE)[1]
  return(path)
}


#' Prepare PTM-SEA input
#' 
#' Write a GCT file with timewise differential analysis t-scores as values and 
#' human flanking sequences as IDs. \code{rat_to_human_ptm_map} is used to map from rat
#' to human phosphosites, and a human protein FASTA file is used to add human 
#' flanking sequences. 
#'
#' @param tissue `r tissue()`. Must be specified if cusotm \code{input} is not provided.   
#' @param fasta [Biostrings::XStringSet] object returned from reading in 
#'   a human protein FASTA file with [Biostrings::readAAStringSet()]. Names of the 
#'   [Biostrings::XStringSet] object should be set to the human protein accession, 
#'   e.g., "Q96QG7". If not specified, the result of [load_uniprot_human_fasta()] is used, which 
#'   returns the version of the UniProt human protein FASTA used in the manuscript. 
#' @param rat_to_human_ptm_map data frame with columns "ptm_id_rat_refseq" and "ptm_id_human_uniprot"
#'   used to map PTMs from rat to human. Default: [MotrpacRatTraining6moData::RAT_TO_HUMAN_PHOSPHO]
#' @param outdir character, output directory for GCT file. The directory is created if it does not already exist. 
#'   Current directory by default. 
#' @param outfile_prefix character, prefix for output GCT file. By default, this prefix
#'   includes the specified tissue and current date. Must be specified for custom input data. 
#' @param input optional data frame if the user wants to perform this analysis for 
#'   a custom set of differential analysis results. Required columns are "feature_ID", "tscore", and \code{cast_vars}. 
#' @param cast_vars character vector of column names in the differential analysis results 
#'   that are used to convert the table from long to wide format, with t-scores as the value variable. 
#'   See [tidyr::pivot_wider()] for more details. Default: "sex", "comparison_group"  
#'
#' @return character, full path to GCT file 
#' @export
#' 
#' @seealso [load_uniprot_human_fasta()], [find_flanks()], [MotrpacRatTraining6moData::RAT_TO_HUMAN_PHOSPHO]
#' 
#' @importFrom purrr map_chr
#' @importFrom methods new
#' @importFrom tibble column_to_rownames remove_rownames
#' @importFrom dplyr filter mutate select group_by summarise_all
#' @importFrom tidyr pivot_wider separate_rows
#' 
#' @examples
#' \dontrun{
#' # Using existing differential analysis results
#' prepare_ptmsea_input("HEART", outdir="/tmp")
#' 
#' # Using a "custom" input
#' res = combine_da_results(tissues = "HEART", assays = "PHOSPHO")
#' # add dummy column
#' prepare_ptmsea_input(input=res, outdir="/tmp", outfile_prefix="HEART_PHOSPHO")
#' }
prepare_ptmsea_input = function(tissue = NULL, 
                                fasta = NULL, 
                                rat_to_human_ptm_map = MotrpacRatTraining6moData::RAT_TO_HUMAN_PHOSPHO,
                                outdir = ".", 
                                outfile_prefix = NULL, 
                                input = NULL, 
                                cast_vars = c("sex","comparison_group")){
  
  missing = c()
  for(pkg in c("stringr", "BiocGenerics", "cmapR")){
    if(!requireNamespace(pkg, quietly = TRUE)){
      missing = c(missing, pkg)
    }
  }
  if(length(missing)>0){
    stop(sprintf("The following packages must be installed to run 'prepare_ptmsea_input()':\n%s",
                 paste(missing, collapse=", ")))
  }
  
  if(is.null(input) & is.null(tissue)){
    stop("If 'input' is not specified, 'tissue' must be specified.")
  }
  
  # define outfile 
  date = gsub("-","",Sys.Date())
  if(is.null(outfile_prefix)){
    if(!is.null(input)){
      stop("If a custom input is provided, 'outfile_prefix' must be defined.")
    }
    outfile_prefix = sprintf("MotrpacRatTraining6mo_ptmsea_%s_%s", tissue, date)
  }
  if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
  }
  
  # load differential analysis results
  if(!is.null(input)){
    dea.results = as.data.frame(input)
    # check for required columns 
    req_cols = c("tscore", "feature_ID", cast_vars)
    missing = c()
    for(c in req_cols){
      if(!c %in% colnames(dea.results)){
        missing = c(missing, c)
      }
    }
    if(length(missing)>0){
      stop(sprintf("The following required columns are missing from the custom input:\n%s",
                   paste0(missing, collapse=", ")))
    }
  }else{
    dea.results = combine_da_results(tissues=tissue, assays="PHOSPHO")
  }
  
  # load data 
  rat2human = as.data.frame(rat_to_human_ptm_map)
  if(is.null(fasta)){
    human.fasta = load_uniprot_human_fasta()
  }else{
    human.fasta = fasta
  }
  
  # annotate with human flanks 
  message("Adding human flanking sequences...")
  rat2human.flanks <- rat2human %>% 
    dplyr::filter(ptm_id_rat_refseq %in% dea.results$feature_ID) %>%
    dplyr::filter(!is.na(ptm_id_human_uniprot)) %>%
    unique %>%
    dplyr::mutate(flanking_sequence = purrr::map_chr(ptm_id_human_uniprot, 
                                                     find_flanks, 
                                                     human.fasta)) %>%
    dplyr::mutate(ptm_id = ptm_id_rat_refseq)
  
  dea.results = merge(dea.results, 
                      rat2human.flanks, 
                      by.x="feature_ID",
                      by.y="ptm_id")
  
  message("Formatting and writing GCT file...")
  tscore.table <- dea.results %>%
    dplyr::mutate(signed_logpval = tscore) %>%
    #Make into wide format
    dplyr::select(feature_ID, flanking_sequence, sex, comparison_group, signed_logpval) %>%
    unique() %>%
    #Expand flanking_sequence that are grouped
    tidyr::pivot_wider(names_from = cast_vars,
                       values_from = signed_logpval) %>%
    #Now we combine rows with the same flanking sequence
    tidyr::separate_rows(flanking_sequence, sep = "\\|") %>%
    dplyr::select(-feature_ID) %>% 
    dplyr::group_by(flanking_sequence) %>% 
    dplyr::summarise_all(mean, na.rm=T) %>%
    #Make modifications uppercase
    dplyr::mutate(flanking_sequence = toupper(flanking_sequence)) %>%
    #Substitute dashes with underscores
    dplyr::mutate(flanking_sequence = gsub("-","_",flanking_sequence)) %>%
    #Add phospho mark to the sequence
    dplyr::mutate(flanking_sequence = paste0(flanking_sequence,"-p")) %>%
    as.data.frame()
  
  tscore.gct <- methods::new("GCT",
                              mat = tscore.table %>% tibble::remove_rownames() %>%
                                tibble::column_to_rownames(var = 'flanking_sequence') %>% 
                                as.matrix(),
                              rdesc = data.frame(id.flanking = as.character(tscore.table$flanking_sequence)))
  tscore.gct@rdesc$id.flanking <- as.character(tscore.gct@rdesc$id.flanking)
  outfile = sprintf("%s/%s.gct", outdir, outfile_prefix)
  cmapR::write_gct(tscore.gct, ofile = outfile, appenddim = F)
  
  path = list.files(path=dirname(outfile), pattern=basename(outfile), full.names=TRUE)[1]
  message("Done.")
  return(path)
}


#' Run ssGSEA2 or PTM-SEA
#' 
#' Run ssGSEA or PTM-SEA with [ssGSEA2.0](https://github.com/broadinstitute/ssGSEA2.0). 
#' GCT files with results are written to \code{outdir}. 
#'
#' @param input_gct_path character, path to GCT file used as input, e.g., from 
#'   [prepare_gsea_input()] or [prepare_ptmsea_input()]
#' @param method character, one of "gsea","gsea_mitocarta","ptmsea". Must be 
#'   "gsea" or "ptmsea" is \code{gene_set_db_path} is specified. 
#' @param gene_set_db_path optional character or character vector of GMT files
#'   that should be used as the gene set database(s). If not specified, the corresponding
#'   database used for manuscript analyses is used (see details). 
#' @param outdir character, output directory for results. The directory is created if it does not already exist. 
#'   Current directory by default. 
#' @param output_prefix character, prefix for results. By default, the prefix is taken
#'   from \code{input_gct_path}. 
#' @param nperm integer, number of permutations. Default 10000 for methods "gsea", "gsea_mitocarta"
#'   and 1000 for method "ptmsea". 
#' @param ... other arguments passed to [ssGSEA2::run_ssGSEA2()], e.g., for parallelization 
#'
#' @return Large named list with full GSEA output. Documentation for this list is
#'   not currently available. Most users will be interested
#'   in the GCT files written to \code{outdir} instead. 
#' @export 
#' 
#' @seealso [prepare_gsea_input()], [prepare_ptmsea_input()], [ssGSEA2::run_ssGSEA2()] 
#' 
#' @details
#' **Gene set databases:**
#' 
#' If \code{gene_set_db_path} is not specified, than the following gene set databases
#' are used for the given method:  
#' 
#' * gsea: \code{c2.cp.v7.0.symbols.RAT.gmt}  
#' * gsea_mitocarta: \code{Rat.Refseq.MitoPathways3.0.gmt}  
#' * ptmsea: \code{ptm.sig.db.all.flanking.human.v1.9.0.gmt}  
#'
#' @examples
#' \dontrun{
#' # For these examples, we use a single permutation for the sake of running
#' # the example quickly. In practice, a large number of permutations should 
#' # be used, e.g., 1000, 10000
#' 
#' # Run GSEA with differential analysis results from liver global proteome
#' input = prepare_gsea_input(tissue = "LIVER",
#'                            assay = "PROT", 
#'                            outdir = "/tmp")
#' res = ssGSEA2_wrapper(input, method = "gsea", outdir = "/tmp", nperm = 1)
#' 
#' # Run GSEA with Mitocarta database,
#' # which uses RefSeq IDs as pathway members
#' input = prepare_gsea_input(tissue = "SKM-GN",
#'                            assay = "TRNSCRPT", 
#'                            outdir = "/tmp",
#'                            gene_id_type = "refseq")
#' res = ssGSEA2_wrapper(input, method = "gsea_mitocarta", outdir = "/tmp", nperm = 1)
#' 
#' # Use example from the Broad
#' # Download example input
#' download.file(url = paste0("https://raw.githubusercontent.com/nicolerg/ssGSEA2.0/",
#'     "master/example/PI3K_pert_logP_n2x23936.gct"),
#'   destfile = "/tmp/PI3K_pert_logP_n2x23936.gct")
#' # Download gene set database
#' download.file(url = paste0("https://raw.githubusercontent.com/nicolerg/ssGSEA2.0/",
#'     "master/example/ptm.sig.db.all.flanking.human.v1.8.1.gmt"),
#'   destfile = "/tmp/ptm.sig.db.all.flanking.human.v1.8.1.gmt")
#' # Using a small number of permutations for the sake of the example. 1000 recommended
#' res = ssGSEA2_wrapper("/tmp/PI3K_pert_logP_n2x23936.gct",
#'                      method = "ptmsea", 
#'                      gene_set_db_path = "/tmp/ptm.sig.db.all.flanking.human.v1.8.1.gmt",
#'                      outdir = "/tmp",
#'                      nperm = 1)
#' }
ssGSEA2_wrapper = function(input_gct_path, 
                          method = c("gsea","gsea_mitocarta","ptmsea"), 
                          gene_set_db_path = NULL,
                          output_prefix = NULL,
                          outdir = ".",
                          nperm = NULL,
                          ...){
  
  if(!requireNamespace("ssGSEA2", quietly = TRUE)){
    stop(
      paste("Package 'ssGSEA2' must be installed to run 'ssGSEA2_wrapper()'.",
            "For R >= 4.0, install with devtools::install_github('nicolerg/ssGSEA2.0').",
            "For R 3.6, see installation instructions: https://github.com/nicolerg/ssGSEA2.0/blob/master/README.md#r-36"),
      call. = FALSE
    )
  }
  
  if(!is.null(gene_set_db_path) & !method %in% c("ptmsea","gsea")){
    stop("If 'gene_set_db_path' is provided, 'method' must be one of 'gsea', 'ptmsea'.")
  }
  
  if(is.null(output_prefix)){
    output_prefix = basename(input_gct_path)
    output_prefix = gsub("\\.[^\\.]*$", "", output_prefix)
  }else{
    if(grepl("/", output_prefix)){
      warning("'output_prefix' should not include paths. Specify the output directory with 'outdir' instead.")
      if(outdir == "."){
        outdir = dirname(output_prefix)
      }
      output_prefix = basename(output_prefix)
    }
  }
  
  if(is.null(gene_set_db_path)){
    if(method == "gsea"){
      gene_set_db_url = "https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/gsea/c2.cp.v7.0.symbols.RAT.gmt"
    }else if(method == "gsea_mitocarta"){
      gene_set_db_url = "https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/gsea/Rat.Refseq.MitoPathways3.0.gmt"
    }else if(method == "ptmsea"){
      gene_set_db_url = "https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/gsea/ptm.sig.db.all.flanking.human.v1.9.0.gmt"
    }else{
      stop("'method' must be one of 'gsea', 'ptmsea', 'gsea_mitocarta'.")
    }
    # download db file 
    gene_set_db_path = sprintf("%s/%s", outdir, basename(gene_set_db_url))
    utils::download.file(url = gene_set_db_url, destfile = gene_set_db_path)
  }
  
  if(is.null(nperm)){
    if(grepl("gsea", method)){
      nperm = 10000
    }else if(method=="ptmsea"){
      nperm = 1000
    }
  }

  #Setting sample.norm.type to none and correl.type to rank uses the actual signed log p-vals
  gsea.out <- ssGSEA2::run_ssGSEA2(
    input.ds = input_gct_path, # input.ds is a path to a file
    output.prefix = output_prefix,
    output.directory = outdir, 
    gene.set.databases = gene_set_db_path, # list of paths to gene set files 
    sample.norm.type = "none",
    weight = 0.75,
    correl.type = "rank",
    statistic = "area.under.RES",
    output.score.type = "NES",
    min.overlap = 5,
    extended.output = TRUE,
    global.fdr = FALSE,
    nperm = nperm,
    log.file = sprintf("%s/run.log", outdir),
    ...)
  
  main_files = list.files(path = outdir, pattern = "scores|pvalues|combined", full.names = TRUE)
  main_files = main_files[grepl(output_prefix, main_files)]
  message(sprintf("Output GCT files written:\n%s", paste(main_files, collapse="\n")))
  
  #Original parameters used for GSEA before 2/16/2021
  # gsea.out <- ssGSEA2(input.ds = input.ds, output.prefix = output.prefix, gene.set.databases = gene.set.databases,
  #                     sample.norm.type = "rank", weight=0.75, correl.type="z.score", statistic="area.under.RES",
  #                     output.score.type="NES", nperm=1000, min.overlap=5, extended.output=T, global.fdr=F)
  # gsea.out <- ssGSEA2(input.ds = input.ds, output.prefix = output.prefix, gene.set.databases = gene.set.databases,
  #                     sample.norm.type = "rank", weight=0, correl.type="z.score", statistic="area.under.RES",
  #                     output.score.type="NES", nperm=1000, min.overlap=5, extended.output=T, global.fdr=F)
  
  #Original parameters used for PTM-SEA before 3/30/2021
  # gsea.out <- ssGSEA2(input.ds = input.ds, output.prefix = output.prefix, gene.set.databases = gene.set.databases,
  #                     sample.norm.type = "rank", weight=0.75, correl.type="z.score", statistic="area.under.RES", 
  #                     output.score.type="NES", nperm=1000, min.overlap=5, extended.output=T, global.fdr=F)
  
  return(gsea.out)
}


#' Load UniProt human canonical protein FASTA file 
#'
#' @param scratchdir character, directory in which the file from Google Cloud
#'   Storage should be downloaded
#'
#' @return [Biostrings::XStringSet] object returned from reading in 
#'   the FASTA file with [Biostrings::readAAStringSet()]
#' @export
#' @seealso [find_flanks()]
#'
#' @examples
#' fasta = load_uniprot_human_fasta()
#' head(fasta)
#' 
#' @details 
#' The human proteome FATSA for canonical proteins was downloaded from UniProt on 
#' 2/3/2021 (UniProtKB query "reviewed:true AND proteome:up000005640"). 
load_uniprot_human_fasta = function(scratchdir="."){
  
  missing = c()
  for(pkg in c("stringr", "BiocGenerics")){
    if(!requireNamespace(pkg, quietly = TRUE)){
      missing = c(missing, pkg)
    }
  }
  if(length(missing)>0){
    stop(sprintf("The following packages must be installed to run 'find_flanks()':\n%s",
                 paste(missing, collapse=", ")))
  }
  
  url = "https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/gsea/uniprot-reviewed-homo_sapiens_20210203.fasta"
  local = sprintf("%s/%s", scratchdir, basename(url))
  download.file(url, destfile = local)
  human.fasta = Biostrings::readAAStringSet(local)
  names(human.fasta) = stringr::str_split_fixed(names(human.fasta),"\\|",n = 3)[,2]
  file.remove(local)
  return(human.fasta)
}


#' Find human flanks
#' 
#' Find flanking sequence in human ortholog of PTM 
#'
#' @param x character, human accession for PTM, e.g., "Q96QG7_S548s"
#' @param fasta [Biostrings::XStringSet] object returned from reading in 
#'   a human protein FASTA file with [Biostrings::readAAStringSet()]. Names of the 
#'   [Biostrings::XStringSet] object should be set to the human protein accession, 
#'   e.g., "Q96QG7". Can be the result of [load_uniprot_human_fasta()], which 
#'   returns the version of the Uniprot human protein FASTA used in the manuscript. 
#' @param type character, type of accession 
#'
#' @return character, pipe-separated string of human flanking sequences for the PTM
#' @export
#' @seealso [load_uniprot_human_fasta()]
#'
#' @examples
#' fasta = load_uniprot_human_fasta(scratchdir="/tmp")
#' x = "Q96QG7_S548s"
#' find_flanks(x, fasta)
find_flanks = function(x, fasta, type = "uniprot"){
  
  missing = c()
  for(pkg in c("stringr", "BiocGenerics", "XVector")){
    if(!requireNamespace(pkg, quietly = TRUE)){
      missing = c(missing, pkg)
    }
  }
  if(length(missing)>0){
    stop(sprintf("The following packages must be installed to run 'find_flanks()':\n%s",
                 paste(missing, collapse=", ")))
  }

  #If missing mapping, return NA
  if(is.na(x)){
    return(NA)
  }
  
  #If this is not a single character object, report error
  if(!is.character(x) | length(x) != 1){
    stop("Input must be a character vector of length 1")
  }
  
  #Extract accession and position or position(s)
  if(type == "uniprot"){
    accession <- stringr::str_split_fixed(x,"_",2)[1]
    position <- stringr::str_split_fixed(x,"_",2)[2] %>%
      stringr::str_extract_all("\\d+") %>% unlist %>%
      as.numeric
  }else{
    accession <- stringr::str_extract(x,"^.._\\d+\\.\\d")
    position <- stringr::str_split_fixed(x,"_",3)[3] %>%
      stringr::str_extract_all("\\d+") %>% unlist %>%
      as.numeric
  }
  
  #Extract sequence of the accession
  current.seq <- fasta[accession]
  
  #Vector with flanking sequences
  output <- c()
  
  #Loop through all aminoacid positions
  for(i in position){
    
    #If the position is less than 8 from the edge, then add underscores
    if(i > BiocGenerics::width(current.seq)){
      #Skip  
      warning(paste("Invalid amino acid location", accession, i))
    }else if (i < 8){ 
      flanks <- stringr::str_c(stringr::str_dup("_",8-i),
                               toString(XVector::subseq(current.seq,1,i+7)))
    }else if(i > BiocGenerics::width(current.seq)-7) {
      flanks <- stringr::str_c(toString(XVector::subseq(current.seq,
                                               i-7,
                                               BiocGenerics::width(current.seq))),
                               stringr::str_dup("_",7-(BiocGenerics::width(current.seq)-i)))
    }else{
      flanks <- XVector::subseq(current.seq,i-7,i+7) %>% as.character
    }
    
    output <- c(output,flanks)
  }
  
  if(is.null(output)){
    return(NA)
  } else {
    return(paste(output,collapse="|"))
  }
}
