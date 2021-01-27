#' Hardy-Weinbery Wrapper
#'
#' Wrapper for main HWE function
#' @param Tab Data frame of genotype files post processing.
#' @param All.ColNames Character vector of Tab object column names.
#' @param Output Data return carryover from main BIGDAWG function
#' @param Verbose Summary display carryover from main BIGDAWG function
#' @note This function is for internal BIGDAWG use only.
HWE.wrapper <- function(Tab,All.ColNames,Output,Verbose) {
  
  HWE <- HWE(Tab,All.ColNames)
  
  if(Output) {
    sink("HWE.txt")
    print(HWE,quote=F)
    sink()
  }
  
  cat("\n> HARDY-WEINBERG ANALYSIS COMPLETED\n")
  if(Verbose) {
    cat("\nControls (Group 0):\n")
    HWE.con <- as.data.frame(HWE[['controls']])
    print(HWE.con,row.names=F,quote=F)
    cat("\nCases (Group 1):\n")
    HWE.cas <- as.data.frame(HWE[['cases']])
    print(HWE.cas,row.names=F,quote=F)
    cat("\n")
  }
  
  return(HWE)
  
}
