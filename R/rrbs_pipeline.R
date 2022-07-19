#' Merge an edgeR DGEList object using a clustering of the sites
#' 
#' @param yall a DGEList object. 
#'        yall$genes is a metadata data frame with the locus coordinates (see details), and these fields at minimum:
#'        Chr, EntrezID, Symbol, and Strand
#' @param new_clusters a character vector, contains the clustering solution of the sites in yall
#' 
#' @return A new DGEList object that represents the clusters
#' @import edgeR,MotrpacRatTraining6moData
#' 
#' @details The yall object has a metadata framework yall$genes. This data frame has either a Locus field or 
#'          a pair of fields (LocStart,LocEnd). Assuming that there are no clusters that merge sites across different
#'          chromosomes and that clusters represent a continuous window in the genome, the function goes over the
#'          clustering solutions in new_clusters and merges the sites from the same cluster.
#'          The new genomic features contain the sum of counts of their sites, and the merged metadata of the sites
#'          (e.g., a comma separated character with all gene symbols associated with the cluster)
#'
#' @examples TODO
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
  newy = rowsum(yall,new_clusters,reorder=T)
  # sanity check: symbols should be the same:
  if(!all(rownames(newy) == rownames(cluster_df))){stop("Error after rowsum")}
  # update the gene annotation
  newy$genes = cluster_df[rownames(newy$counts),]
  # reorder
  o <- order(newy$genes$Chr, newy$genes$LocStart)
  newy <- newy[o,]
  return(newy)
}


#' A wrapper function that analyzes a tile of loci in the genome to obtain correlated clusters.
#' 
#' @param tile_name
#' @param tile_l
#' @param M
#' @param min_cor
#' @param inflations
#' @param plotcorr
#' 
#' @return 
#' 
#' @import MotrpacRatTraining6moData, MCL
#' 
#' @details 
#'        MCL(Markov clustering) details: from: https://www.micans.org/mcl/intro.html
#'        MCL: Expansion coincides with taking the power of a stochastic matrix
#'        using the normal matrix product (i.e. matrix squaring). 
#'        Inflation corresponds with taking the Hadamard power of a matrix 
#'        (taking powers entrywise), followed by a scaling step, such that the resulting matrix is stochastic again, 
#'        i.e. the matrix elements (on each column) correspond to probability values.
#'        Inflation parameter: strengthen intra-region connections and promote cluster homogeneity.
#' 
#' @examples 
#' # load the RData file from the cloud
#' system("wget xxx")
#' yall = get(load(xxx))
#' # remove control samples
#' is_sample = grepl("^9",colnames(yall),perl=T)
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
#' print(paste("Data from tissue",tissue,"was loaded into session"))
#' print(paste("how many sites were removed (unassembled+chrM)?",sum(!keep)))
#' # keep.lib.sizes=FALSE causes the library sizes to be recomputed
#' yall <- yall[keep,, keep.lib.sizes=FALSE]
#' # fix the levels and order by chromosome
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
#' yall$counts[, grepl("-Un",colnames(yall$counts))]
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
#' TSS <- nearestTSS(yall$genes$Chr, yall$genes$Locus, species="Rn")
#' yall$genes$EntrezID <- TSS$gene_id
#' yall$genes$Symbol <- TSS$symbol
#' yall$genes$Strand <- TSS$strand
#' yall$genes$Distance <- TSS$distance
#' yall$genes$Width <- TSS$width
#' 
#' Now we are ready to cluster the genome
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
#' # single core
#' #tile_clusters = lapply(names(tile_l),analyze_tile,tile_l=tile_l,M=M)
#' # parallel comp
#' tile_clusters = copy(tiles)
#' for(chr in unique(chrs)){
#'   print(chr)
#'   curr_chr_inds = grepl(paste0(chr,"-"),tiles)
#'   chr_tiles = unique(tiles[curr_chr_inds])
#'   chr_tile_l = tile_l[chr_tiles]
#'   print(system.time({
#'       chr_tile_clusters = mclapply(chr_tiles,analyze_tile,
#'           tile_l = chr_tile_l,M = M, mc.cores = 5)
#'       chr_tile_clusters = unlist(chr_tile_clusters)
#'       tile_clusters[names(chr_tile_clusters)] = chr_tile_clusters}))
#'   gc()
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
#' # Other normalization methods developed for RNA-seq data, such as TMM, are not required for BS-seq data.
#' TotalLibSize <- yall$samples$lib.size[grepl("-Me",colnames(yall$counts))] +
#'     +       yall$samples$lib.size[grepl("-Un",colnames(yall$counts))]
#' yall$samples$lib.size <- rep(TotalLibSize, each=2)
analyze_tile<-function(tile_name,tile_l,M,min_cor=0.7,
                       inflations = c(3,2.5,2,1.5),plotcorr=F){
  #print(tile_name)
  tile_inds = tile_l[[tile_name]]
  if(length(tile_inds) < 2){
    v = tile_name
    names(v) = rownames(M)[tile_inds]
    return(v)
  }
  tile_M = M[tile_inds,]
  tile_corrs = cor(t(tile_M))
  if(plotcorr){
    corrplot(tile_corrs)
  }
  tile_clusters = NULL
  for(infl in inflations){
    tr = tryCatch({
      tile_clusters = mcl(tile_corrs > min_cor,addLoops = T,allow1=T,
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
