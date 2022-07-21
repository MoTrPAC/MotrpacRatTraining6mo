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


#' Download an RData file from Dropbox
#' 
#' @param target character, public URL specifying an RData file ('.rda' suffix) on Dropbox
#' @param tmpdir character, local path in which to download \code{target}
#' @param redownload boolean, whether or not to download the file if it already exists in \code{tmpdir}
#'
#' @return object contained in \code{target}
#'
#' @examples
#' target = "https://www.dropbox.com/s/rsenv5dgsx0to4i/test.rda?raw=1"
#' obj = load_dropbox_rdata(target, tmpdir)
#' obj = load_dropbox_rdata(target, tmpdir, redownload = T)
load_dropbox_rdata = function(target, tmpdir, redownload = FALSE){
  
  # This is intended to be an internal function only
  accepted_files = c("https://www.dropbox.com/s/rsenv5dgsx0to4i/test.rda?raw=1",
                     "https://www.dropbox.com/s/s6c4otggrkyhh05/ATAC_NORM_DATA.rda?raw=1",
                     "https://www.dropbox.com/s/t92nm792snxzi3c/METHYL_NORM_DATA.rda?raw=1",
                     "https://www.dropbox.com/s/b2wxfvb5q59xtw0/METHYL_RAW_COUNTS.rda?raw=1")
  if(!target %in% accepted_files){
    warning(sprintf("Unrecognized Dropbox link. 'target' is expected to be one of:\n%s",
                    paste0(accepted_files, collapse="\n")))
  }
  
  # Check that it ends in "?raw=1"
  if(!endsWith(target, "?raw=1")){
    stop("'target' should end with '?raw=1'.")
  }
  
  # Check that it ends in ".rda?raw=1"
  if(!endsWith(gsub("\\?raw=1", "", basename(target)), ".rda")){
    stop("'target' should be an RData file with suffix '.rda'.")
  }
  
  dest = sprintf("%s/%s", tmpdir, gsub("\\?raw=1","",basename(target)))
  
  if(!file.exists(dest) | redownload){
    download.file(target,
                  destfile = dest,
                  method = "auto")
  }
  
  if(target=="https://www.dropbox.com/s/s6c4otggrkyhh05/ATAC_NORM_DATA.rda?raw=1"){
    message("Loading and returning an object of size 3.31GB...")
  }else if (target=="https://www.dropbox.com/s/b2wxfvb5q59xtw0/METHYL_RAW_COUNTS.rda?raw=1"){
    message("Loading and returning an object of size 3.70GB...")
  }else if(target=="https://www.dropbox.com/s/t92nm792snxzi3c/METHYL_NORM_DATA.rda?raw=1"){
    message("Loading and returning an object of size 3.70GB...")
  }
  
  load(dest)
  obj_name = gsub("\\.rda","",basename(dest))
  if(!exists(obj_name)){
    stop(sprintf("File '%s' doesn't contain an object called '%s'. 'target' must specify an RData file that contains a single object with the same name as the file.", basename(dest), obj_name))
  }
  
  return(get(obj_name))
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


#' TODO
#' Internal function used to check arguments for differential analysis functions
#'
#' @param tissue tissue abbreviation, one of [MotrpacRatTraining6moData::TISSUE_ABBREV]
#' @param outfile 
#' @param overwrite 
#' @param outfile_is_rdata 
#'
#' @return NULL
#'
#' @examples
#' TODO 
#' 
check_dea_args = function(tissue, outfile, overwrite, outfile_is_rdata = TRUE){
  # check arguments 
  if(length(tissue)>1){
    stop("Please specify a single tissue, e.g., 'BAT' for brown adipose tissue. See 'TISSUE_ABBREV' for options.")
  }
  if(!is.null(rdata_outfile)){
    if(!overwrite & file.exists(rdata_outfile)){
      stop(sprintf("'%s' already exists and 'overwrite' = FALSE. Specify a new outfile or set 'overwrite' to TRUE to generate new results.",
                   outfile))
    }
    if(outfile_is_rdata){
      if(!any(unlist(lapply(c("\\.rda$", "\\.rdata"), function(pattern){
        grepl(pattern, rdata_outfile, ignore.case = TRUE)
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
