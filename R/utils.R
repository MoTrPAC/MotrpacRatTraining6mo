#' Map vial labels to PIDs
#' 
#' Map sample identifiers (vial labels) to participant IDs (PIDs).
#' Vial labels are unique per sample; PIDs are unique per animal. 
#' Multiple vial labels correspond to different aliquots of the same tissue sample. 
#'
#' @param viallabels vector, list of vial labels  
#'
#' @return named list where names are vial labels and values are PIDs 
#' @export
#' 
#' @import data.table 
#' @import MotrpacRatTraining6moData
#'
#' @examples
#' viallabel_to_pid(c("90416015402", "90416015403", "90416015302"))
#' viallabel_to_pid(c(90416015402, 90416015403, 90416015302))
viallabel_to_pid = function(viallabels){
  pheno = as.data.table(MotrpacRatTraining6moData::PHENO)
  pheno = unique(pheno[,.(viallabel, pid)])
  pheno = pheno[viallabel %in% as.character(viallabels)]
  vl_to_pid = pheno[,pid]
  names(vl_to_pid) = pheno[,viallabel]
  return(vl_to_pid)
}


#' List available data
#' 
#' List available data, including lazily-loaded data. Useful to gather
#' split data frames.
#'
#' @param package optional string to specify a package 
#'
#' @return character vector of names of data objects available to load
#' @export
#'
#' @examples
#' list_available_data()
#' list_available_data("MotrpacBicQC")
list_available_data = function(package=NULL){
  res = data(package=package)
  obj = res$results[,3]
  # remove objects that can't be called directly
  obj = obj[!grepl("\\(", obj)]
  return(obj)
}


#' Check arguments for DEA functions 
#' 
#' Internal function used to check arguments for differential analysis functions
#'
#' @param tissue `r tissue()`
#' @param outfile character, output file 
#' @param overwrite bool, whether to overwrite \code{outfile} if it exists
#' @param outfile_is_rdata bool, whether \code{outfile} is intended to save RData
#'
#' @return NULL
#'
#' @examples
#' check_dea_args("BAT", "~/Desktop/bat.rda", overwrite = FALSE) 
#' 
check_dea_args = function(tissue, outfile, overwrite, outfile_is_rdata = TRUE){
  # check arguments 
  if(length(tissue)>1){
    stop("Please specify a single tissue, e.g., 'BAT' for brown adipose tissue. See 'TISSUE_ABBREV' for options.")
  }
  if(!is.null(outfile)){
    if(!overwrite & file.exists(outfile)){
      stop(sprintf("'%s' already exists and 'overwrite' = FALSE. Specify a new outfile or set 'overwrite' to TRUE to generate new results.",
                   outfile))
    }
    if(outfile_is_rdata){
      if(!any(unlist(lapply(c("\\.rda$", "\\.rdata"), function(pattern){
        grepl(pattern, outfile, ignore.case = TRUE)
      })))){
        stop(sprintf("Outfile '%s' should end with '.rda' or '.RData' to indicate an RData file (not case-sensitive).", outfile))
      }
    }
    # check if path is valid
    dir = dirname(outfile)
    if(dir != "" & !dir.exists(dir)){
      message(sprintf("Creating directory for outfile: %s", dir))
      dir.create(dir, recursive = TRUE) # this will throw an error if it's not a valid path/uncreatable 
    }
  }
  return()
}
