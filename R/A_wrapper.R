#' Amino Acid Wrapper
#'
#' Wrapper function for amino acid analysis.
#' @param loci Loci being analyzed.
#' @param loci.ColNames The column names of the loci being analyzed.
#' @param genos Genotype table.
#' @param grp Case/Control or Phenotype groupings.
#' @param EPL Exon protein alignment.
#' @param Cores Number of cores to use for analysis.
#' @param Strict.Bin Logical specify if strict rare cell binning should be used in ChiSq test
#' @param Output Data return carryover from main BIGDAWG function
#' @param Verbose Summary display carryover from main BIGDAWG function
#' @note This function is for internal BIGDAWG use only.
A.wrapper <- function(loci,loci.ColNames,genos,grp,EPL,Cores,Strict.Bin,Output,Verbose) {
  
  cat("\n>>>> STARTING AMINO ACID LEVEL ANALYSIS...\n")
  
  # Define Lists for Per Loci Running Tallies
  AAlog <-  list()
  AminoAcid.binned <- list()
  AminoAcid.freq <- list()
  overall.chisq <- list()
  ORtable <- list()
  Final_binned <- list()
  
  cat("Processing Locus: ")
  
  # Loop Through Loci
  for(x in loci) {
    
    # Get Locus
    Locus <- x
    cat(Locus,".. ")
    
    # Read in Locus Alignment file for Locus specific exons
    ExonAlign <- EPL[[Locus]]; rownames(ExonAlign) <- NULL
    
    # Run Amino Acid Analysis
    A.list <- A(Locus,loci.ColNames,genos,grp,Strict.Bin,ExonAlign,Cores)
    
    # Build Output Lists
    AAlog[[Locus]] <- A.list[['log']]
    AminoAcid.binned[[Locus]] <- A.list[['binned']]
    overall.chisq[[Locus]] <- A.list[['chisq']]
    ORtable[[Locus]] <- A.list[['OR']]
    Final_binned[[Locus]] <- A.list[['table']]
    AminoAcid.freq[[Locus]] <- A.list[['freq']]
    
  }; rm(x) #locus loop
  
  cat("\n\n")
  
  Out <- list()
  Out[['AL']] <- do.call(rbind,AAlog)
  Out[['AB']] <- do.call(rbind,AminoAcid.binned)
  Out[['AF']] <- do.call(rbind,AminoAcid.freq)
  Out[['CS']] <- do.call(rbind,overall.chisq)
  Out[['OR']] <- do.call(rbind,ORtable)
  Out[['FB']] <- do.call(rbind,Final_binned)
  
  if(Output) {
    ## write to file
    write.table(Out[['AL']], file = "AA_log.txt", sep="\t", row.names = F, col.names=T, quote = F)
    write.table(Out[['AF']], file = "AA_freqs.txt", sep="\t", row.names = F, col.names=T, quote = F)
    write.table(Out[['AB']], file = "AA_binned.txt", sep="\t", row.names = F, col.names=T, quote = F)
    write.table(Out[['OR']], file = "AA_OR.txt", sep="\t", row.names = F, col.names=T, quote = F)
    write.table(Out[['CS']], file = "AA_chisq.txt", sep="\t", row.names = F, col.names=T, quote = F)
    write.table(Out[['FB']], file = "AA_table.txt", sep="\t", row.names = F, col.names=T, quote = F)
  }
  
  cat("> AMINO ACID ANALYSIS COMPLETED\n")
  
  if(Verbose){
    cat("Significant Amino Acid Position(s):","\n")
    tmp <- do.call(rbind,overall.chisq); rownames(tmp) <- NULL
    tmp.sig <- tmp[which(tmp[,'sig']=="*"),]; rownames(tmp.sig) <- NULL
    if(nrow(tmp.sig)>0) { print(tmp.sig,row.names=F,quote=F) }
  }
  
  
  return(Out)
  
}