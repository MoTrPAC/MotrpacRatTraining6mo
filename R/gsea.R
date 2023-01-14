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
  gct = new("GCT", mat = mat)
  
  # write to file 
  cmapR::write_gct(gct, ofile = sprintf("%s/%s.gct", outdir, outfile_prefix), appenddim = FALSE)
  return(sprintf("%s/%s.gct", outdir, outfile_prefix))
}


prepare_ptmsea_input = function(tissue){
  
  # #Path to feature-mapping file in local copy of data freeze bucket
  # mapping.dir <- "~/Projects/motrpac-data-freeze-pass/pass1b-06/v1.0/analysis/proteomics-untargeted/prot-ph/normalized-data/"
  # #Path to results in local copy of data freeze bucket
  # dea.dir <- "~/Projects/motrpac-data-freeze-pass/pass1b-06/v1.0/analysis/proteomics-untargeted/prot-ph/dea/"
  # dea.result.file <- paste0("pass1b-06_",current.tissue,"_prot-ph_timewise-dea-fdr_20211006.txt")
  # dea.annot.file <- paste0("motrpac_pass1b-06_",current.tissue,"_prot-ph_normalized-data-feature-annot.txt")
  # #Rat to human phosphosite mapping table
  # rat2human <- read_csv("~/Projects/motrpac-pass/mawg-data/mappings/proteomics/motrpac_pass1b-06_proteomics-ph-rat2human-20211016.csv")
  rat2human = get("RAT_TO_HUMAN_PHOSPHO", envir=as.environment("package:MotrpacRatTraining6moData"))
  # #FASTA database used in searches. Will be used to add flanking sequences if the column VMSiteFlanks is not available.
  # #Rat FASTA
  # # db.str <- '~/motrpac-pass/fasta/RefSeq.20181127_rat_full.fasta'
  # #Human FASTA
  # human.db.str <- '~/Projects/motrpac-pass/mawg-data/pass1b-06/mappings/proteomics/mappings_proteomics_uniprot-reviewed-homo_sapiens_20210203.fasta'
  # 
  # output.dir <- "~/Projects/motrpac-pass/mawg-data/pass1b-06/prot-ph/gct_ptm-centric/"
  # sign.logp.filename <- paste0("pass1b-06_",current.tissue,"_prot-ph",prot_corr_suffix,"_signed-logp_ptm-centric.gct")
  # sign.logp.human.filename <- paste0("pass1b-06_",current.tissue,"_prot-ph",prot_corr_suffix,"_signed-logp_ptm-centric_human-orthologs.gct")
  # zscore.filename <- paste0("pass1b-06_",current.tissue,"_prot-ph",prot_corr_suffix,"_tscore_ptm-centric.gct")
  # zscore.human.filename <- paste0("pass1b-06_",current.tissue,"_prot-ph",prot_corr_suffix,"_tscore_ptm-centric_human-orthologs.gct")
  
  # dea.results <- read_tsv(Sys.glob(file.path(dea.dir,dea.result.file)))
  # dea.annot <- read_tsv(Sys.glob(file.path(mapping.dir,dea.annot.file))) 
  # dea.results <- left_join(dea.results,dea.annot) 
  dea.results = combine_da_results(tissues=tissue, assays="PHOSPHO")
  dea.annot = get("FEATURE_TO_GENE_FILT", envir=as.environment("package:MotrpacRatTraining6moData"))
  dea.results = merge(dea.results, dea.annot, by="feature_ID")
  
  # where is the flanking_sequence column coming from?
  
  sign.logp.table <- dea.results %>%
    mutate(signed_logpval = tscore) %>%
    #Make into wide format
    select(feature_ID,flanking_sequence,sex,comparison_group,signed_logpval)%>%
    pivot_wider(names_from = c("sex","comparison_group"),values_from = signed_logpval) %>%
    #Expand flanking_sequence that are grouped
    separate_rows(flanking_sequence,sep = "\\|") %>%
    #Now we combine rows with the same flanking sequence
    select(-feature_ID) %>% group_by(flanking_sequence) %>% summarise_all(mean,na.rm=T) %>%
    #Make modifications uppercase
    mutate(flanking_sequence = toupper(flanking_sequence)) %>%
    #Substitute dashes with underscores
    mutate(flanking_sequence = gsub("-","_",flanking_sequence)) %>%
    #Add phospho mark to the sequence
    mutate(flanking_sequence = paste0(flanking_sequence,"-p")) %>%
    as.data.frame()
  sign.logp.gct <- new("GCT",
                       mat = sign.logp.table %>% remove_rownames() %>%
                         column_to_rownames(var = 'flanking_sequence') %>% as.matrix(),
                       rdesc = data.frame(id.flanking = as.character(sign.logp.table$flanking_sequence)))
  sign.logp.gct@rdesc$id.flanking <- as.character(sign.logp.gct@rdesc$id.flanking)
  write_gct(sign.logp.gct,ofile = paste0(output.dir,zscore.filename),appenddim = F)
  
  # unclear if there's more to add here?
  # next section says "#3. Create PTM-SEA input with human flanking sequences",
  # but it only uses signed logFCs, and the methods say we used t-scores
}

gsea = function(tissue, assay, gene_set, gene_centric = TRUE, outdir = "."){
  
  if(!requireNamespace("ssGSEA2", quietly = TRUE)){
    stop(
      "Package 'ssGSEA2' must be installed to run 'gsea()'.",
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
  gsea.out <- run_ssGSEA2(input.ds = input.ds, # input.ds is a path to a file
                          output.prefix = output.prefix,
                          gene.set.databases = gene.set.databases, # list of paths to gene set files 
                          sample.norm.type = "none",
                          weight=0.75,
                          correl.type="rank",
                          statistic="area.under.RES",
                          output.score.type="NES",
                          min.overlap=5,
                          extended.output=TRUE,
                          global.fdr=FALSE,
                          nperm=10000,
                          log.file=sprintf("%s/run.log", outdir))

  #Original parameters used before 2/16/2021
  # gsea.out <- ssGSEA2(input.ds = input.ds, output.prefix = output.prefix, gene.set.databases = gene.set.databases,
  #                     sample.norm.type = "rank", weight=0.75, correl.type="z.score", statistic="area.under.RES",
  #                     output.score.type="NES", nperm=1000, min.overlap=5, extended.output=T, global.fdr=F)
  # gsea.out <- ssGSEA2(input.ds = input.ds, output.prefix = output.prefix, gene.set.databases = gene.set.databases,
  #                     sample.norm.type = "rank", weight=0, correl.type="z.score", statistic="area.under.RES",
  #                     output.score.type="NES", nperm=1000, min.overlap=5, extended.output=T, global.fdr=F)
}

ptmsea = function(tissue){
  
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

  current.wd <- getwd() 
  setwd(output.dir)
  #Setting sample.norm.type to none and correl.type to rank uses the actual signed log p-vals
  gsea.out <- run_ssGSEA2(input.ds = input.ds, 
                          output.prefix = output.prefix, 
                          gene.set.databases = gene.set.databases,
                          sample.norm.type = "none", 
                          weight=0.75, 
                          correl.type="rank", 
                          statistic="area.under.RES",
                          output.score.type="NES", 
                          nperm=1000, 
                          min.overlap=5, 
                          extended.output=TRUE, 
                          global.fdr=FALSE)
  #Original parameters used before 3/30/2021
  # gsea.out <- ssGSEA2(input.ds = input.ds, output.prefix = output.prefix, gene.set.databases = gene.set.databases,
  #                     sample.norm.type = "rank", weight=0.75, correl.type="z.score", statistic="area.under.RES", 
  #                     output.score.type="NES", nperm=1000, min.overlap=5, extended.output=T, global.fdr=F)
  setwd(current.wd)
}
