#' Return list of accepted tissue abbreviations 
#' 
#' Return tissue abbreviations for the 20 tissues analyzed in the manuscript. 
#' 
#' @param alphabetical_order bool, whether to alphabetize the tissue abbreviations. 
#'     \code{TRUE} by default. If \code{FALSE}, return them in the biological order they 
#'     appear in manuscript figures. 
#'     
#' @return a vector of strings 
#' 
#' @examples
#' tissue_abbreviations()
#' tissue_abbreviations(alphabetical_order = FALSE)
#' 
#' @export
#' 
tissue_abbreviations = function(alphabetical_order = TRUE){
  tissue_abbr = c('BLOOD',
                  'PLASMA',
                  'HEART',
                  'VENACV',
                  'SPLEEN',
                  'SKM-GN',
                  'SKM-VL',
                  'WAT-SC',
                  'BAT',
                  'LIVER',
                  'LUNG',
                  'KIDNEY',
                  'ADRNL',
                  'CORTEX',
                  'HYPOTH',
                  'HIPPOC',
                  'SMLINT',
                  'COLON',
                  'OVARY',
                  'TESTES')
  if(alphabetical_order){
    return(tissue_abbr[order(tissue_abbr)])
  }else{
    return(tissue_abbr)
  }
}
