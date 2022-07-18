#' Merge an edgeR DGEList object using a clustering of the sites
#' 
#' @param yall a DGEList object. 
#'        yall$genes is a metadata data frame with the locus coordinates (see details), and these fields at minimum:
#'        Chr, EntrezID, Symbol, and Strand
#' @param new_clusters a character vector, contains the clustering solution of the sites in yall
#' 
#' @return A new DGEList object that represents the clusters
#' @import edgeR
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

# from: https://www.micans.org/mcl/intro.html
# MCL: Expansion coincides with taking the power of a stochastic matrix
# using the normal matrix product (i.e. matrix squaring). 
# Inflation corresponds with taking the Hadamard power of a matrix 
# (taking powers entrywise), followed by a scaling step, 
# such that the resulting matrix is stochastic again, 
# i.e. the matrix elements (on each column) correspond to probability values.
# Inlation: strengthen intra-region connections

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
#' @details 
#' 
#' @examples TODO
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
