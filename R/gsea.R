#' Prepare GCT file for GSEA
#' 
#' Take a set of differential analysis results, either for a specific tissue and ome or
#' from a user-supplied table, convert to GCT format, and write to file for GSEA.  
#'
#' @param tissue `r tissue()` 
#' @param assay character, assay abbreviation, one of "PROT", "PHOSPHO", "TRNSCRPT", "ACETYL", "UBIQ" 
#' @param outdir character, output directory for GCT file. The directory is created if it does not already exist. 
#'   Current directory by default. 
#' @param outfile_prefix character, prefix for output GCT file. By default, this prefix
#'   includes the specified tissue and assay and current date. Must be specified for custom input data. 
#' @param input optional data frame if the user wants to perform this analysis for 
#'   a custom set of results. Required columns are "tscore", "gene_symbol", and \code{cast_vars}. 
#' @param cast_vars character vector of column names in the differential analysis results 
#'   that are used to convert the table from long to wide format, with t-scores as the value variable. 
#'   See [data.table::dcast()] for more details. Default: "sex", "comparison_group"  
#'
#' @return character, path of the GCT file  
#' @export
#' @importFrom methods new
#'
#' @examples
#' prepare_gsea_input("HEART","PROT",outdir="/tmp")
#' prepare_gsea_input("LIVER","PHOSPHO",outdir="/tmp")
#' prepare_gsea_input("LIVER","ACETYL",outdir="/tmp")
#' prepare_gsea_input("LIVER","UBIQ",outdir="/tmp")
#' prepare_gsea_input("LIVER","TRNSCRPT",outdir="/tmp")
#' 
#' # "custom" input
#' res = combine_da_results(tissues = "KIDNEY", assays = "PROT")
#' # add dummy column
#' res$gene_symbol = res$feature_ID
#' prepare_gsea_input(input=res, outdir="/tmp", outfile_prefix="KIDNEY_PROT")
#' 
#' @details 
#' T-scores from the timewise differential analysis results are used for scores. 
#' Feature-level data is summarized into gene-level data using the maximum absolute t-score. 
#' 
prepare_gsea_input = function(tissue = NULL, assay = NULL, outdir = ".", outfile_prefix = NULL, input = NULL, cast_vars = c("sex","comparison_group")){
  
  if(!requireNamespace("cmapR", quietly = TRUE)){
    stop(
      "Package 'cmapR' must be installed to run 'da_results_to_gct()'.",
      call. = FALSE
    )
  }
  
  if(is.null(input) & any(is.null(c(tissue, assay)))){
    stop("If 'input' is not specified, both 'tissue' and 'assay' must be specified.")
  }
  
  # define outfile 
  date = gsub("-","",Sys.Date())
  if(is.null(outfile_prefix)){
    if(!is.null(input)){
      stop("If a custom input is provided, 'outfile_prefix' must be defined.")
    }
    outfile_prefix = sprintf("MotrpacRatTraining6mo_%s_%s_%s", tissue, assay, date)
  }
  if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
  }
  
  # load differential analysis results
  if(!is.null(input)){
    da = data.table::as.data.table(input)
    # check for required columns 
    req_cols = c("tscore", "gene_symbol", cast_vars)
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
  }else{
    # get feature metadata
    feature_map = data.table::data.table(get("FEATURE_TO_GENE_FILT", envir=as.environment("package:MotrpacRatTraining6moData")))
    if(assay %in% c("PROT","PHOSPHO","TRNSCRPT","ACETYL","UBIQ")){
      curr_feature_map = feature_map
    }else{
      stop(sprintf("GSEA/PTM-SEA not compatible with %s differential analysis results.", assay))
    }
    da = data.table::as.data.table(combine_da_results(tissues = tissue, assays = assay))
    
    # merge with differential analysis results
    merged = merge(da, curr_feature_map, by="feature_ID")
  }
  
  # handle column names 
  if(!"tscore" %in% colnames(merged)){
    merged[,tscore := NA_real_]
  }
  if(!"zscore" %in% colnames(merged)){
    merged[,zscore := NA_real_]
  }
  
  # make table of scores 
  merged[is.na(tscore), tscore := zscore]
  merged = merged[!is.na(tscore), .(feature_ID, gene_symbol, sex, comparison_group, tscore)]
  # get max t-score per gene 
  merged_max = merged[!is.na(gene_symbol),list(max_tscore = max(tscore)),
                      by=.(gene_symbol, sex, comparison_group)]
  # cast wide
  form = sprintf("gene_symbol ~ %s", paste0(cast_vars, collapse=" + "))
  merged_wide = data.table::dcast(merged_max, formula = eval(parse(text=form)), value.var='max_tscore')
  
  # convert to matrix
  rn = merged_wide[,gene_symbol]
  merged_wide[,gene_symbol := NULL]
  mat = as.matrix(merged_wide)
  colnames(mat) = colnames(merged_wide)
  rownames(mat) = rn
  
  # make GCT
  gct = methods::new("GCT", mat = mat)
  
  # write to file 
  cmapR::write_gct(gct, ofile = sprintf("%s/%s.gct", outdir, outfile_prefix), appenddim = FALSE)
  return(sprintf("%s/%s.gct", outdir, outfile_prefix))
}


#' Prepare PTM-SEA input
#'
#' @param tissue `r tissue()`
#' @param fasta [Biostrings::XStringSet] object returned from reading in 
#'   a human protein FASTA file with [Biostrings::readAAStringSet()]. Names of the 
#'   [Biostrings::XStringSet] object should be set to the human protein accession, 
#'   e.g., "Q96QG7". If not specified, the result of [load_uniprot_human_fasta()] is used, which 
#'   returns the version of the Uniprot human protein FASTA used in the manuscript. 
#' @param outfile character, optional output file name for GCT file. By default, 
#'   the file name includes the tissue and date and is saved to the current working directory.  
#'
#' @return character, full path to GCT file 
#' @export
#' 
#' @seealso [load_uniprot_human_fasta()], [find_flanks()]
#' 
#' @importFrom purrr map_chr
#' @importFrom methods new
#' @importFrom tibble column_to_rownames remove_rownames
#' @importFrom dplyr filter mutate select group_by summarise_all
#' @importFrom tidyr pivot_wider separate_rows
#' 
#' @examples
#' res = prepare_ptmsea_input("HEART")
prepare_ptmsea_input = function(tissue, fasta=NULL, outfile=NULL){
  
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
  
  if(is.null(outfile)){
    outfile = sprintf("ptmsea_input_tscore_%s_%s.gct",
                      tissue,
                      gsub("-","",Sys.Date()))
  }
  
  # load data 
  rat2human = MotrpacRatTraining6moData::RAT_TO_HUMAN_PHOSPHO
  dea.results = combine_da_results(tissues=tissue, assays="PHOSPHO")
  dea.annot = MotrpacRatTraining6moData::FEATURE_TO_GENE_FILT
  dea.results = merge(dea.results, dea.annot, by="feature_ID")
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
    tidyr::pivot_wider(names_from = c("sex","comparison_group"),
                       values_from = signed_logpval) %>%
    #Expand flanking_sequence that are grouped
    tidyr::separate_rows(flanking_sequence, sep = "\\|") %>%
    #Now we combine rows with the same flanking sequence
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
  cmapR::write_gct(tscore.gct, ofile = outfile, appenddim = F)
  
  path = list.files(path=dirname(outfile), pattern=basename(outfile), full.names=TRUE)[1]
  message("Done.")
  return(path)
}


gsea = function(tissue, assay, gene_set, gene_centric = TRUE, outdir = "."){
  
  if(!requireNamespace("ssGSEA2", quietly = TRUE)){
    stop(
      paste("Package 'ssGSEA2' must be installed to run 'gsea()'.",
            "Install with devtools::install_github('nicolerg/ssGSEA2.0')."),
      call. = FALSE
    )
  }

  # #ssGSEA script - https://github.com/broadinstitute/ssGSEA2.0
  # script.dir <- "~/Projects/ssGSEA2.0/"
  # source("~/Projects/ssGSEA2.0/src/ssGSEA2.0.R")
  # library(tidyverse)
  current.tissue = tissue
  current.ome = assay
  
  input.ds = da_results_to_gct(current.tissue, 
                               current.ome, 
                               outdir = outdir, 
                               outfile_prefix = NULL, 
                               input = NULL, 
                               cast_vars = c("sex","comparison_group"))

  #ssGSEA parameters
  if(gene_centric){
    extension = "_gene-centric.gct"
  } else{
    extension = ".gct"
  }
  # #Input data file in GMT format
  # if(current.ome == "prot-ph"){
  #   input.ds <- paste0("~/Projects/motrpac-pass/mawg-data/pass1b-06/",current.ome,"/GCT_files/pass1b-06_",current.tissue,"_",
  #                      current.ome,"_protein-centric_signed_logp",extension)
  # 
  # } else if(current.ome == "transcript-rna-seq"){
  #   input.ds <- paste0("~/Projects/motrpac-pass/mawg-data/pass1b-06/",current.ome,"/GCT_files/pass1b-06_",current.tissue,"_",
  #                      current.ome,"_signed_logp",extension)
  # 
  # } else if(current.ome == "prot-ac"){
  #   input.ds <- paste0("~/Projects/motrpac-pass/mawg-data/pass1b-06/",current.ome,"/GCT_files/pass1b-06_",current.tissue,"_",
  #                      current.ome,"_protein-centric_signed_logp",extension)
  # } else {
  #   input.ds <- paste0("~/Projects/motrpac-pass/mawg-data/pass1b-06/",current.ome,"/GCT_files/pass1b-06_",current.tissue,"_",
  #                      current.ome,"_signed_logp",extension)
  # }
  # cat("reading file: ", input.ds)
  # #Output prefix
  # output.prefix <- paste0("pass1b-06_",current.tissue,"_",current.ome)
  # #Gene set directory
  # gene.set.dir <- params$gene_set_dir
  # #List of gene sets to evaluate
  # gene.set.file <- params$gene_set
  # #Paste gene set and directories
  # gene.set.databases <- as.list(paste0(gene.set.dir,gene.set.file))
  # #Output directory
  # if(is.null(params$output_dir)){
  #   output.dir <- paste0("~/Projects/motrpac-pass/mawg-data/pass1b-06/",current.ome,"/gsea/",
  #                        gene.set.file,"/",current.tissue,"/")
  # } else {
  #   output.dir <- file.path(params$output_dir,current.tissue)
  # }
  # dir.create(output.dir,recursive=T)

  # ssGSEA analysis
  
  ## recommended to run directly on the command line
  # current.wd <- getwd()
  # setwd(output.dir)
  
  #Setting sample.norm.type to none and correl.type to rank uses the actual signed log p-vals
  gsea.out <- ssGSEA2::run_ssGSEA2(
    input.ds = input.ds, # input.ds is a path to a file
    output.prefix = output.prefix,
    gene.set.databases = gene.set.databases, # list of paths to gene set files 
    sample.norm.type = "none",
    weight = 0.75,
    correl.type = "rank",
    statistic = "area.under.RES",
    output.score.type = "NES",
    min.overlap = 5,
    extended.output = TRUE,
    global.fdr = FALSE,
    nperm = 10000,
    log.file = sprintf("%s/run.log", outdir))

  #Original parameters used before 2/16/2021
  # gsea.out <- ssGSEA2(input.ds = input.ds, output.prefix = output.prefix, gene.set.databases = gene.set.databases,
  #                     sample.norm.type = "rank", weight=0.75, correl.type="z.score", statistic="area.under.RES",
  #                     output.score.type="NES", nperm=1000, min.overlap=5, extended.output=T, global.fdr=F)
  # gsea.out <- ssGSEA2(input.ds = input.ds, output.prefix = output.prefix, gene.set.databases = gene.set.databases,
  #                     sample.norm.type = "rank", weight=0, correl.type="z.score", statistic="area.under.RES",
  #                     output.score.type="NES", nperm=1000, min.overlap=5, extended.output=T, global.fdr=F)
}

ptmsea = function(tissue){
  
  if(!requireNamespace("ssGSEA2", quietly = TRUE)){
    stop(
      paste("Package 'ssGSEA2' must be installed to run 'ptmsea()'.",
            "Install with devtools::install_github('nicolerg/ssGSEA2.0')."),
      call. = FALSE
    )
  }
  
  # #ssGSEA script - https://github.com/broadinstitute/ssGSEA2.0
  # script.dir <- "~/ssGSEA2.0/"
  # source("~/ssGSEA2.0/src/ssGSEA2.0.R")
  # #Helper tools
  # if(file.exists("~/helper-functions/helper.functions.v1.1.R")){
  #   source("~/helper-functions/helper.functions.v1.1.R")
  # }else{
  #   devtools::source_url("https://raw.githubusercontent.com/pierremj/helper-functions/master/helper.functions.v1.1.R")
  # }
  # library(limma)
  # library(gprofiler2)
  # library(tidyverse)
  
  #ssGSEA parameters
  
  #These are the parameters that should be changed to determine input/output

  #Add suffix if one wishes to use human ortholog PTM sites
  if(params$use_human_flanks){
    suffix <- "_human_orthologs"
  } else {
    suffix <- ""
  }
  #Tissue selection
  current.tissue <- params$tissue_code
  #Construct path to input file
  input.ds <- paste0("~/motrpac-pass/mawg-data/pass1b-06/prot-ph/GCT_files/pass1b-06_",current.tissue,
                     "_prot-ph_signed_logp",suffix,".gct")
  cat("Reading file:", input.ds)
  #Output prefix
  output.prefix <- paste0("pass1b-06_",current.tissue,"_prot-ph",suffix) 
  #Gene set directory
  gene.set.dir <- params$gene_set_dir
  #List of gene sets to evaluate
  gene.set.file <- params$gene_set
  #List of gene sets to evaluate
  gene.set.databases <- as.list(paste0(gene.set.dir,gene.set.file))
  #Output directory
  output.dir <- paste0("~/motrpac-pass/mawg-data/pass1b-06/prot-ph/gsea/",
                       gene.set.file,"/",current.tissue,"/")
  dir.create(output.dir,recursive=T)
  
  #1. PTM-SEA analysis
  #Recommended to run directly on the command line
  
  # only difference from GSEA params: nperm = 1000 instead of 10000
    
  current.wd <- getwd() 
  setwd(output.dir)
  #Setting sample.norm.type to none and correl.type to rank uses the actual signed log p-vals
  gsea.out <- ssGSEA2::run_ssGSEA2(
    input.ds = input.ds, 
    output.prefix = output.prefix, 
    gene.set.databases = gene.set.databases,
    sample.norm.type = "none", 
    weight = 0.75, 
    correl.type = "rank", 
    statistic = "area.under.RES",
    output.score.type = "NES", 
    nperm = 1000, 
    min.overlap = 5, 
    extended.output = TRUE, 
    global.fdr = FALSE,
    log.file = sprintf("%s/run.log", outdir))
  #Original parameters used before 3/30/2021
  # gsea.out <- ssGSEA2(input.ds = input.ds, output.prefix = output.prefix, gene.set.databases = gene.set.databases,
  #                     sample.norm.type = "rank", weight=0.75, correl.type="z.score", statistic="area.under.RES", 
  #                     output.score.type="NES", nperm=1000, min.overlap=5, extended.output=T, global.fdr=F)
  setwd(current.wd)
}


#' Load Uniprot human protein FASTA file 
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
  
  url = "https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/misc/uniprot-reviewed-homo_sapiens_20210203.fasta"
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
