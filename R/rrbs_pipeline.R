#' Merge sites by cluster
#' 
#' Merge an \code{edgeR} [edgeR::DGEList()] object using a clustering of the sites.
#' 
#' @param yall A [edgeR::DGEList()] object, where yall$genes is a metadata data frame with the locus coordinates (see details), 
#'   and these fields at minimum: "Chr", "EntrezID", "Symbol", and "Strand".
#' @param new_clusters A character vector. Contains the clustering solution of the sites in \code{yall}.
#' 
#' @return A new [edgeR::DGEList()] object that represents the clusters
#' 
#' @importFrom edgeR DGEList
#' 
#' @export
#' 
#' @details 
#' The \code{yall} object has a metadata framework \code{yall$genes}. This data frame has either a "Locus" field or 
#' a pair of fields ("LocStart", "LocEnd"). Assuming that there are no clusters that merge sites across different
#' chromosomes and that clusters represent a continuous window in the genome, the function goes over the
#' clustering solutions in new_clusters and merges the sites from the same cluster.
#' The new genomic features contain the sum of counts of their sites, and the merged metadata of the sites
#' (e.g., a comma separated character with all gene symbols associated with the cluster).
#' See example of the full RRBS read count data pre-processing pipeline in [analyze_tile()].
#' 
#' @seealso [analyze_tile()]
merge_sites_by_clusters<-function(yall,new_clusters){
  if (!"LocStart" %in% colnames(yall$genes)){
    cl_start = tapply(yall$genes$Locus,new_clusters,min)
  }
  else{
    cl_start = tapply(yall$genes$LocStart,new_clusters,min)
  }
  if (!"LocEnd" %in% colnames(yall$genes)){
    cl_end = tapply(yall$genes$Locus,new_clusters,max)
  }
  else{
    cl_end = tapply(yall$genes$LocEnd,new_clusters,max)
  }
  # assume that we do not merge sites across different chormosomes
  cl_chr = as.character(tapply(yall$genes$Chr,new_clusters,function(x)x[1]))
  cl_entrez = tapply(yall$genes$EntrezID,new_clusters,unique)
  cl_entrez = sapply(cl_entrez,paste,collapse=",")
  cl_symbol = tapply(yall$genes$Symbol,new_clusters,unique)
  cl_symbol = sapply(cl_symbol,paste,collapse=",")
  cl_strand = tapply(yall$genes$Strand,new_clusters,unique)
  cl_strand = sapply(cl_strand,paste,collapse=",")
  cl_loc = paste(cl_start,cl_end,sep="-")
  cl_sites = tapply(yall$genes$Locus,new_clusters,unique)
  cl_sites = sapply(cl_sites,paste,collapse=",")
  if (!"NumSites" %in% colnames(yall$genes)){
    cl_num_sites = tapply(yall$genes$Strand,new_clusters,length)
  }
  else{
    cl_num_sites = tapply(yall$genes$NumSites,new_clusters,sum)
  }
  # sanity checks: orders must be the same
  if(!(all(names(cl_end)==names(cl_start))) ||
     !(all(names(cl_end)==names(cl_symbol))) ||
     !(all(names(cl_end)==names(cl_entrez))) ||
     !(all(names(cl_end)==names(cl_strand)))
  ){
    stop("Error in tapply by cluster name")
  }
  cluster_df = data.frame(
    Chr = cl_chr, 
    Locus = cl_loc,
    EntrezID = cl_entrez,
    Symbol = cl_symbol,
    Strand = cl_strand,
    Width = cl_end - cl_start,
    NumSites = cl_num_sites,
    Sites = cl_sites,
    LocStart = cl_start, LocEnd = cl_end
  )
  rownames(cluster_df) = names(cl_start)
  # sum counts in regions
  newy = rowsum(yall,new_clusters,reorder=TRUE)
  # sanity check: symbols should be the same:
  if(!all(rownames(newy) == rownames(cluster_df))){stop("Error after rowsum")}
  # update the gene annotation
  newy$genes = cluster_df[rownames(newy$counts),]
  # reorder
  o <- order(newy$genes$Chr, newy$genes$LocStart)
  newy <- newy[o,]
  return(newy)
}


#' Analyze genome tiles 
#' 
#' A wrapper function that analyzes a tile of loci in the genome to obtain correlated clusters.
#' 
#' @param tile_name A character. The name of the current window for analysis.
#' @param tile_l A named list. tile_name must be an item in the list. Each entry contains the M matrix row indices of the loci of the window.
#' @param M A numeric matrix. Contains the "M" methylation values. Rows are loci, columns are samples.
#' @param min_cor A number. The minimal correlation to be considered for connecting a pair of loci.
#' @param inflations A numeric vector. Internal parameters of the MCL algorithm, see Details.
#' @param plotcorr A logical. TRUE: plot the correlation matrix of the window (useful for analysis of specific windows).
#' 
#' @return A character vector. Names are loci (M matrix rows), entries are the names of the loci clusters.
#' 
#' @export
#' 
#' @details 
#' This function implements clustering analysis of a specific window in the genome. 
#' We observed that often adjacent loci in RRBS data may manifest low correlations. This
#' function can be used for identified a set of highly correlated loci within a window, which
#' can then be used to extract new genomic features for downstream analysis. 
#' RRBS data have two columns per sample: "Un" (unmethylated) and "Me"(methylated). This function
#' takes the "M matrix" which contains the integrated values of these columns (see examples). This matrix
#' is used for computing the correlations of the loci, and then we use the Markov Clustering algorithm (MCL)
#' for identifying homogeneous clusters within a window. 
#' MCL(Markov clustering) details from: https://www.micans.org/mcl/intro.html
#' MCL: Expansion coincides with taking the power of a stochastic matrix
#' using the normal matrix product (i.e. matrix squaring). 
#' Inflation corresponds with taking the Hadamard power of a matrix 
#' (taking powers entrywise), followed by a scaling step, such that the resulting matrix is stochastic again, 
#' i.e. the matrix elements (on each column) correspond to probability values.
#' Inflation parameter: strengthen intra-region connections and promote cluster homogeneity.
#' 
#' @examples 
#' \dontrun{
#' # Raw data in RData file is available through Google Cloud.  
#' # The main URL is https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS
#' # The files that are available through this URL are by tissue:
#' # Brown adipose:
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/BAT_raw.RData
#' # Heart: 
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/HEART_raw.RData
#' # Hippocampus: 
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/HIPPOC_raw.RData
#' # Kidney: 
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/KIDNEY_raw.RData
#' # Lung: 
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/LUNG_raw.RData
#' # Liver: 
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/LIVER_raw.RData
#' # Gastrocnemius: 
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/SKMGN_raw.RData
#' # White adipose: 
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/WATSC_raw.RData
#' # download the gastrocnemius data and load the data object into this session
#' f = "https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/SKMGN_raw.RData"
#' system(sprintf("wget %s", f))
#' yall = get(load("SKMGN_raw.RData"))
#' 
#' # TODO: add links to the data using GCP URLs
#' 
#' # remove control samples
#' is_sample = grepl("^9",colnames(yall),perl=TRUE)
#' yall = yall[,is_sample]
#' 
#' # Filtering unassembled chromosomes
#' keep <- rep(TRUE, nrow(yall))
#' Chr <- as.character(yall$genes$Chr)
#' keep[ grep("random",Chr) ] <- FALSE
#' keep[ grep("chrUn",Chr) ] <- FALSE
#' # remove non-chr ones
#' keep[!grepl("chr",Chr) ] <- FALSE
#' # remove M chromosome (otherwise we get error when assigning annotation below)
#' keep[Chr=="chrM"] <- FALSE
#' # keep.lib.sizes=FALSE causes the library sizes to be recomputed
#' yall <- yall[keep,, keep.lib.sizes=FALSE]
#' # rat genome chromosomes:
#' ChrNames <- paste0("chr",c(1:20,"X","Y"))
#' yall$genes$Chr <- factor(yall$genes$Chr, levels=ChrNames)
#' o <- order(yall$genes$Chr, yall$genes$Locus)
#' yall <- yall[o,]
#' print(paste("Counts matrix dim before low counts filter:"))
#' print(dim(yall))
#' gc()
#' 
#' # Remove low count features
#' # The analysis needs to be restricted to CpG sites that have enough coverage for the
#' # methylation level to be measurable in a meaningful way at that site. 
#' # As a conservative rule of thumb, we require a site to have a total count
#' # (both methylated and unmethylated) of at least 8 in every sample.
#' Coverage <- yall$counts[,grepl("-Me",colnames(yall$counts))] + 
#'      yall$counts[, grepl("-Un",colnames(yall$counts))]
#' min_count_coverage = 10
#' min_per_samples = 1
#' # see description of min_count_coverage and min_per_samples above
#' keep <- rowSums(Coverage >= min_count_coverage) >= min_per_samples*ncol(Coverage)
#' # filter the data
#' yall <- yall[keep,, keep.lib.sizes=FALSE]
#' print(paste("Counts matrix dim after low counts filter:"))
#' print(dim(yall))
#' 
#' # get locations, gene ids, etc
#' TSS <- edgeR::nearestTSS(yall$genes$Chr, yall$genes$Locus, species="Rn")
#' yall$genes$EntrezID <- TSS$gene_id
#' yall$genes$Symbol <- TSS$symbol
#' yall$genes$Strand <- TSS$strand
#' yall$genes$Distance <- TSS$distance
#' yall$genes$Width <- TSS$width
#' 
#' # Now we are ready to cluster the genome
#' wsize = 500
#' # step 1: define the genome-level tiles
#' chrs = as.character(yall$genes$Chr)
#' pos = round(yall$genes$Locus/wsize)
#' tiles = paste(chrs,pos,sep="-")
#' names(tiles) = rownames(yall$counts)
#' # step 2: add the M matrix for correlations
#' Me <- yall$counts[, grepl("-Me",colnames(yall$counts))]
#' Un <- yall$counts[, grepl("-Un",colnames(yall$counts))]
#' M <- log2(Me + 2) - log2(Un + 2)
#' # step 3: analyze each tile
#' tile_l = split(1:nrow(M),tiles)
#' tile_l = tile_l[unique(tiles)] # keep the genome order
#' # parallel comp
#' if(require("parallel")){
#'   tile_clusters = data.table::copy(tiles)
#'   for(chr in unique(chrs)){
#'     print(chr)
#'     curr_chr_inds = grepl(paste0(chr,"-"),tiles)
#'     chr_tiles = unique(tiles[curr_chr_inds])
#'     chr_tile_l = tile_l[chr_tiles]
#'     print(system.time({
#'       chr_tile_clusters = parallel::mclapply(chr_tiles,
#'                                              analyze_tile,
#'                                              tile_l = chr_tile_l,
#'                                              M = M, 
#'                                              mc.cores = 5)
#'       chr_tile_clusters = unlist(chr_tile_clusters)
#'       tile_clusters[names(chr_tile_clusters)] = chr_tile_clusters}))
#'     gc()
#'   }
#' }
#' yall = merge_sites_by_clusters(yall,tile_clusters)
#' 
#' # Data normalization
#' # A key difference between BS-seq and other sequencing data is that the pair of libraries 
#' # holding the methylated and unmethylated reads for a particular sample are treated as a unit.
#' # To ensure that the methylated and unmethylated reads for the same sample are treated on the
#' # same scale, we need to set the library sizes to be equal for each pair of libraries. 
#' # We set the library sizes for each sample to be the average of the total read counts for 
#' # the methylated and unmethylated libraries.
#' # Other normalization methods developed for RNA-seq data, 
#' # such as TMM, are not required for BS-seq data.
#' TotalLibSize <- yall$samples$lib.size[grepl("-Me",colnames(yall$counts))] + 
#'   yall$samples$lib.size[grepl("-Un",colnames(yall$counts))]
#' yall$samples$lib.size <- rep(TotalLibSize, each=2)
#' }
analyze_tile <- function(tile_name,
                         tile_l,
                         M,
                         min_cor=0.7,
                         inflations = c(3,2.5,2,1.5),
                         plotcorr=FALSE){
  
  for(pkg in c("corrplot","MCL")){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        sprintf("Package '%s' must be installed to run 'analyze_tile()'."),
        call. = FALSE
      )
    }
  }
  
  tile_inds = tile_l[[tile_name]]
  if(length(tile_inds) < 2){
    v = tile_name
    names(v) = rownames(M)[tile_inds]
    return(v)
  }
  tile_M = M[tile_inds,]
  tile_corrs = stats::cor(t(tile_M))
  if(plotcorr){
    corrplot::corrplot(tile_corrs)
  }
  tile_clusters = NULL
  for(infl in inflations){
    tr = tryCatch({
      tile_clusters = MCL::mcl(tile_corrs > min_cor,
                               addLoops = TRUE,
                               allow1=TRUE,
                               inflation = infl)$Cluster
    }, error = function(x){}
    )
    if(!is.null(tile_clusters)){break}
  }
  if(is.null(tile_clusters)){
    print("Warning: MCL failed, using singletons instead")
    tile_clusters = 1:length(tile_inds)
  }
  
  # Solve mcl's bug of skipped cluster numbering
  tile_clusters = as.numeric(ordered(tile_clusters))
  tile_clusters = paste0(tile_name,"_cluster",tile_clusters)
  names(tile_clusters) = rownames(tile_M)
  return(tile_clusters)
}


#' RRBS differential anlaysis
#' 
#' Timewise and training differential analysis wrapper for RRBS data.
#' 
#' @param y A [edgeR::DGEList()] object. 
#'   yall$genes is a metadata data frame with the locus coordinates (see details), and these fields at minimum:
#'   Chr, EntrezID, Symbol, and Strand
#' @param PHENO A data frame with a row per sample. Contains at least the following columns: "sex", "group".
#'   [MotrpacRatTraining6moData::PHENO] by default. 
#' @param METHYL_META A data frame with a row per sample. 
#'   Contains the RRBS pipeline QA/QC scores. Contains at least the following columns: "pct_Unaligned".
#'   [MotrpacRatTraining6moData::METHYL_META] by default. 
#' @param verbose A logical. TRUE: comments about the pipeline progress are printed.
#' @param samples_to_remove A character vector. Contains the ids of the samples that should be removed
#'   (e.g., identified outliers or failed samples). 
#'   METHYL samples in [MotrpacRatTraining6moData::OUTLIERS] by default. 
#' @param edger_tol A number. An internal parameter of edgeR. Default is 1e-05. Consider increasing if the algorithm takes too long.
#' @param dataset_name A character. The name of the current dataset. Will be added to the output.
#' @param adj_pct_unaligned A logical. TRUE: adjust for percent unaligned reads. Default is FALSE.
#' 
#' @return A list.First item is called timewise and will contain the contrast-specific differential analysis results.
#' Second item is called training and contains the overall training-level significance per locus.
#' 
#' @importFrom metap sumlog
#' @importFrom edgeR estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom limma makeContrasts zscoreT
#' @importFrom stats model.matrix
#' 
#' @export
#'  
#' @examples
#' \dontrun{
#' data(PHENO)
#' data(METHYL_META)
#' 
#' # Raw data in RData file is available through Google Cloud.  
#' # The main URL is https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS
#' # The files that are available through this URL are by tissue:
#' # Brown adipose: 
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/BAT_raw.RData
#' # Heart: 
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/HEART_raw.RData
#' # Hippocampus: 
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/HIPPOC_raw.RData
#' # Kidney: 
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/KIDNEY_raw.RData
#' # Lung: 
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/LUNG_raw.RData
#' # Liver: 
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/LIVER_raw.RData
#' # Gastrocnemius: 
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/SKMGN_raw.RData
#' # White adipose: 
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/WATSC_raw.RData
#' 
#' # Example: download the gastrocnemius data and load the data object into this session
#' f = "https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/raw/RRBS/SKMGN_raw.RData"
#' system(sprintf("wget %s", f))
#' yall = get(load("SKMGN_raw.RData"))
#'
#' # for the simplicity of this example, we subset the data to 5000 loci
#' y = yall[1:5000,]
#' dea_res = rrbs_differential_analysis(y,PHENO,METHYL_META,adj_pct_unaligned=T)
#' head(dea_res$timewise)
#' head(dea_res$training)
#' 
#' # Alternatively, you can use the processed datasets.
#' # These are also available through the Google Cloud directory:
#' # https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/
#' # File formats are `METHYL_XXX_NORM_DATA.rda` and `METHYL_XXX_DA.rda` 
#' # where XXX is the tissue name:
#' # BAT (brown adipose), HEART, HIPPOC (hippocampus), KIDNEY, LUNG, LIVER, 
#' # SKMGN (gastrocnemius skeletal muscle), and WATSC (white adipose)
#'
#' f = "https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_SKMGN_NORM_DATA.rda"
#' system(sprintf("wget %s", f))
#' y_processed = get(load("METHYL_SKMGN_NORM_DATA.rda"))
#' # again, for the simplicity of this example, we subset the data to 5000 loci
#' y = y_processed[1:5000,-c(1:4)]
#' da_res = rrbs_differential_analysis(y, adj_pct_unaligned=TRUE)
#' head(da_res$timewise)
#' head(da_res$training)
#' }
rrbs_differential_analysis = function(y,
                                      PHENO = MotrpacRatTraining6moData::PHENO,
                                      METHYL_META = MotrpacRatTraining6moData::METHYL_META,
                                      verbose=TRUE,
                                      adj_pct_unaligned=FALSE,
                                      samples_to_remove=data.table::data.table(
                                        MotrpacRatTraining6moData::OUTLIERS)[assay=="METHYL",viallabel],
                                      edger_tol=1e-05,
                                      dataset_name=""){
  
  #remove outlier samples before DA
  s_me=paste(samples_to_remove, "Me", sep="-")
  s_un=paste(samples_to_remove, "Un", sep="-")
  y_input_samples = sapply(colnames(y),function(x)strsplit(x,split="-")[[1]][1])
  y <- y[ , ! colnames(y) %in% c(s_me,s_un)]
  curr_samps = sapply(colnames(y),function(x)strsplit(x,split="-")[[1]][1])
  outlier_samples = setdiff(y_input_samples,curr_samps)
  if (length(outlier_samples) == 0){
    remove_samples = ""
  } else{
    remove_samples = paste(outlier_samples,collapse=",")
  }
  
  # Parse the covariates
  PHENO$timepoint = PHENO$group
  PHENO$is_control = PHENO$group=="control"
  rownames(METHYL_META) = as.character(METHYL_META$vial_label)
  covs = data.frame(
    sample = as.factor(curr_samps),
    Me = grepl("Me",colnames(y)),
    sex = PHENO[curr_samps,"sex"],
    timepoint = PHENO[curr_samps,"timepoint"],
    is_control = PHENO[curr_samps,"is_control"],
    pct_unaligned_1 = stats::poly(METHYL_META[curr_samps,"pct_Unaligned"],1)[,1],
    pct_unaligned_2 = stats::poly(METHYL_META[curr_samps,"pct_Unaligned"],2)[,2]
  )
  rownames(covs) = colnames(y)
  tmptp = covs$timepoint
  tmptp[tmptp=="control"] = "8"
  tmptp = gsub("w","",tmptp)
  covs$group_tp=factor(
    paste(tmptp,as.numeric(covs$is_control),sep="_")
  )
  
  #split covariates to males and females
  covs_males = covs[covs$sex == "male", ] 
  covs_females = covs[covs$sex == "female", ]
  samples_males=rownames(covs_males)
  samples_females=rownames(covs_females)
  
  # subset y into males and females
  y_males = y[ ,samples_males]
  y_females = y[ ,samples_females]
  
  #create design matrices for for males and females separately
  # refactor the samples to have the correct levels
  covs_males$sample = as.factor(as.character(covs_males$sample))
  covs_females$sample = as.factor(as.character(covs_females$sample))
  
  # create design matrices for time-wise and F tests
  # (see https://support.bioconductor.org/p/12119/ and the limma guide) 
  #design matrix for males  
  samples_mat_male = stats::model.matrix(~0+sample,data=covs_males)
  tp_mat_male = stats::model.matrix(~0+group_tp,data=covs_males)
  if(adj_pct_unaligned){
    tp_mat_male = stats::model.matrix(~0+group_tp+pct_unaligned_1+pct_unaligned_2,data=covs_males)
  }
  tp_mat_male[!covs_males$Me,] = 0
  des_males = cbind(samples_mat_male,tp_mat_male)
  #make sure the rows of design matrix and columns of counts are identical
  if(sum(!row.names(des_males) == colnames(y_males$counts))>0){
    stop(paste("ERROR in dataset",dataset_name,"check sample order before moving on"))
  }
  
  #design matrix for females
  samples_mat_female = stats::model.matrix(~0+sample,data=covs_females)
  tp_mat_female = stats::model.matrix(~0+group_tp,data=covs_females)
  if(adj_pct_unaligned){
    tp_mat_female = stats::model.matrix(~0+group_tp+pct_unaligned_1+pct_unaligned_2,data=covs_females)
  }
  tp_mat_female[!covs_females$Me,] = 0
  des_females = cbind(samples_mat_female,tp_mat_female)
  if(sum(!row.names(des_females) == colnames(y_females$counts))>0){
    stop(paste("ERROR in dataset",dataset_name,"check sample order before moving on"))
  }
  
  if(adj_pct_unaligned){
    full_model_str = "~0+sample+1me+group_me + stats::poly(pct_unaligned,2)"
    null_model_str = "~0+sample+1me + stats::poly(pct_unaligned,2)"
  }
  else{
    full_model_str = "~0+sample+1me+group_me"
    null_model_str = "~0+sample+1me"
  }
  
  # define the contrasts for males and females the analyses below
  C_ttests_female = limma::makeContrasts(
    group_tp1_0 - group_tp8_1,group_tp2_0 - group_tp8_1,
    group_tp4_0 - group_tp8_1,group_tp8_0 - group_tp8_1,
    levels = des_females
  )
  
  C_ttests_male = limma::makeContrasts(
    group_tp1_0 - group_tp8_1,group_tp2_0 - group_tp8_1,
    group_tp4_0 - group_tp8_1,group_tp8_0 - group_tp8_1,
    levels = des_males
  )
  
  #Estimate dispersions for males and females separately
  print("Estimating dispersion for the male design matrix...")
  y1_males <- edgeR::estimateDisp(y_males, design=des_males,tol = edger_tol)
  if(verbose){
    print("Done")
    print("Running glmQLFit...")
  }
  fit.ttest_males <- edgeR::glmQLFit(y1_males,des_males)
  if(verbose){print("Done")}
  
  print("Estimating dispersion for the female design matrix...")
  y1_females <- edgeR::estimateDisp(y_females, design=des_females,tol = edger_tol)
  if(verbose){
    print("Done")
    print("Running glmQLFit...")
  }
  fit.ttest_females <- edgeR::glmQLFit(y1_females,des_females)
  if(verbose){print("Done")}
  
  # extract contrast info
  if(verbose){print("Extracting timewise results from the contrasts...")}
  #### male contrast
  tissue_tp_res = c()
  for(col in colnames(C_ttests_male)){
    sex_str = "male"
    curr_tp = strsplit(col,split="_")[[1]][2]
    curr_tp = strsplit(curr_tp,split="tp")[[1]][2]
    if(verbose){print(paste("Fitting timewise model for:",col))}
    res <- edgeR::glmQLFTest(fit.ttest_males,contrast=C_ttests_male[,col])
    # extract the results into a table
    edger_res <- edgeR::topTags(res, n=Inf, p.value = 1,
                         adjust.method = "BH",sort.by = "none")$table
    # add z-scores
    edger_res$F[edger_res$F < 0] = 1e-10
    t.stat <- sign(edger_res$logFC) * sqrt(edger_res$F)
    z <- limma::zscoreT(t.stat, df=res$df.total)
    curr_res = data.frame(
      feature_ID = rownames(y_males),
      edger_res[,1:4],
      assay = "epigen-rrbs",
      tissue = dataset_name,
      removed_samples = remove_samples,
      sex = sex_str,
      logFC = edger_res$logFC,
      fscore = edger_res$F,
      zscore = z,
      covariates = NA,
      comparison_group = paste0(curr_tp,"w"),
      p_value = edger_res$PValue,
      adj_p_value = edger_res$FDR
    )
    # add the results
    tissue_tp_res = rbind(tissue_tp_res,curr_res)
  }
  if(verbose){print("Done")}
  ##### female contrast
  for(col in colnames(C_ttests_female)){
    sex_str = "female"
    curr_tp = strsplit(col,split="_")[[1]][2]
    curr_tp = strsplit(curr_tp,split="tp")[[1]][2]
    if(verbose){print(paste("Fitting timewise model for:",col))}
    res <- edgeR::glmQLFTest(fit.ttest_females,contrast=C_ttests_female[,col])
    # extract the results into a table
    edger_res <- edgeR::topTags(res, 
                                n=Inf, 
                                p.value = 1,
                                adjust.method = "BH",
                                sort.by = "none")$table
    # add z-scores
    edger_res$F[edger_res$F < 0] = 1e-10
    t.stat <- sign(edger_res$logFC) * sqrt(edger_res$F)
    z <- limma::zscoreT(t.stat, df=res$df.total)
    curr_res = data.frame(
      feature_ID = rownames(y_females),
      edger_res[,1:4],
      assay = "epigen-rrbs",
      tissue = dataset_name,
      removed_samples = remove_samples,
      sex = sex_str,
      logFC = edger_res$logFC,
      fscore = edger_res$F,
      zscore = z,
      covariates = NA,
      comparison_group = paste0(curr_tp,"w"),
      p_value = edger_res$PValue,
      adj_p_value = edger_res$FDR
    )
    # add the results
    tissue_tp_res = rbind(tissue_tp_res,curr_res)
  }
  if(verbose){print("Done")}

  #males
  if(verbose){print("Fitting the model for the F-tests...")}
  ftest_tp_mat_male = stats::model.matrix(~1+group_tp,data=covs_males)
  if(adj_pct_unaligned){
    ftest_tp_mat_male = stats::model.matrix(~1+pct_unaligned_1+pct_unaligned_2+group_tp,data=covs_males)
  }
  ftest_tp_mat_male[!covs_males$Me,] = 0
  ftest_des_males = cbind(samples_mat_male,ftest_tp_mat_male)
  y2_males <- edgeR::estimateDisp(y_males, design=ftest_des_males,tol = edger_tol)
  fit.ftest_males <- edgeR::glmQLFit(y2_males,ftest_des_males)
  is_group_variable = grepl("group",colnames(ftest_des_males))
  res <- edgeR::glmQLFTest(fit.ftest_males,coef=colnames(ftest_des_males)[is_group_variable])
  ftest_edger_res <- edgeR::topTags(res, 
                                    n=Inf, 
                                    p.value = 1,
                                    adjust.method = "BH",
                                    sort.by ="none")$table
  if(verbose){print("Done")}
  # add the results
  curr_f_test_res_males = data.frame(
    feature_ID = rownames(y_males),
    ftest_edger_res[,1:4],
    assay = "epigen-rrbs",
    tissue = dataset_name,
    #sex = "male",
    p_value_male = ftest_edger_res$PValue,
    #adj_p_value = ftest_edger_res$FDR,
    fscore_male = ftest_edger_res$F,
    full_model = full_model_str,
    reduced_model = null_model_str
  )
  
  #females
  if(verbose){print("Fitting the model for the F-tests for females...")}
  ftest_tp_mat_female = stats::model.matrix(~1+group_tp,data=covs_females)
  if(adj_pct_unaligned){
    ftest_tp_mat_female = stats::model.matrix(~1+pct_unaligned_1+pct_unaligned_2+group_tp,data=covs_females)
  }
  ftest_tp_mat_female[!covs_females$Me,] = 0
  ftest_des_females = cbind(samples_mat_female,ftest_tp_mat_female)
  y2_females <- edgeR::estimateDisp(y_females, design=ftest_des_females,tol = edger_tol)
  fit.ftest_females <- edgeR::glmQLFit(y2_females,ftest_des_females)
  is_group_variable = grepl("group",colnames(ftest_des_females))
  res <- edgeR::glmQLFTest(fit.ftest_females,coef=colnames(ftest_des_females)[is_group_variable])
  ftest_edger_res <- edgeR::topTags(res, 
                                    n=Inf, 
                                    p.value = 1,
                                    adjust.method = "BH",
                                    sort.by ="none")$table
  if(verbose){print("Done")}
  curr_f_test_res_females = data.frame(
    feature_ID = rownames(y_females),
    ftest_edger_res[,1:4],
    assay = "epigen-rrbs",
    tissue = dataset_name,
    #sex = "female",
    p_value_female = ftest_edger_res$PValue,
    #adj_p_value = ftest_edger_res$FDR,
    fscore_female = ftest_edger_res$F,
    full_model = full_model_str,
    reduced_model = null_model_str
  )
  
  ###combine curr_f_test_res from males and females
  if(nrow(curr_f_test_res_females)!=nrow(curr_f_test_res_males)){
    stop(paste("ERROR in dataset",dataset_name,"male and female ftest tables have different nrow"))
  }
  if(any(curr_f_test_res_females$feature_ID != curr_f_test_res_males$feature_ID)){
    stop(paste("ERROR in dataset",dataset_name,"male and female ftest tables have different row names"))
  }
  common_col_names <- intersect(names(curr_f_test_res_males), names(curr_f_test_res_females))
  curr_f_test_res_combined = merge(curr_f_test_res_males,curr_f_test_res_females,
                                   by=common_col_names, all=TRUE)
  
  #Add sumlog p-value
  curr_f_test_res_combined$p_value = 
    apply(curr_f_test_res_combined[,c("p_value_male","p_value_female")],1,function(x)metap::sumlog(x)$p)
  
  # Add removed samples column
  curr_f_test_res_combined$removed_samples = remove_samples
  
  return(list(
    timewise = tissue_tp_res,
    training = curr_f_test_res_combined
  ))
}
