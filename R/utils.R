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
