
#' Amino Acid Analysis Function
#'
#' This is the workhorse function for the amino acid analysis.
#' @param Locus Locus being analyzed.
#' @param loci.ColNames The column names of the loci being analyzed.
#' @param genos Genotype table.
#' @param grp Case/Control or Phenotype groupings.
#' @param Strict.Bin Logical specify if strict rare cell binning should be used in ChiSq test.
#' @param ExonAlign Exon protein alignment filtered for locus.
#' @param Cores Number of cores to use for analysis.
#' @note This function is for internal BIGDAWG use only.
A <- function(Locus,loci.ColNames,genos,grp,Strict.Bin,ExonAlign,Cores) {

  # pull out locus specific columns
  getCol <- seq(1,length(loci.ColNames),1)[loci.ColNames %in% Locus]
  HLA_grp <- cbind(grp,genos[,getCol])
  rownames(HLA_grp) <- NULL
  nAllele <- length(na.omit(HLA_grp[,2])) + length(na.omit(HLA_grp[,3]))

  ## extract alleles and counts for Grp1 and Grp0
  ## Use 2F because it represent unique protein molecule.
  Alleles <- unique(c(HLA_grp[,2],HLA_grp[,3]))
  Alleles2F <- sort(unique(as.character(sapply(Alleles, GetField, Res=2))))

  if( length(Alleles)>1 ) {

    # Filter exon alignment matrix for specific alleles
    TabAA <- AlignmentFilter(ExonAlign,Alleles2F,Locus)
    TabAA.list <- lapply(seq_len(ncol(TabAA)), function(x) TabAA[,c(2,x)])
    TabAA.list <- TabAA.list[6:ncol(TabAA)]
    TabAA.names <- colnames(TabAA)[6:ncol(TabAA)]

    # Generate List of Contingency Tables
    ConTabAA.list <- parallel::mclapply(TabAA.list,AAtable.builder,y=HLA_grp,mc.cores=Cores)
    names(ConTabAA.list) <- TabAA.names

    # Apply Contingency Table Checker
    FlagAA.list <- parallel::mclapply(ConTabAA.list,AA.df.check,mc.cores=Cores)

    # Run ChiSq
    csRange <- which(FlagAA.list==FALSE)
    if(Strict.Bin) {
      ChiSqTabAA.list <- parallel::mclapply(ConTabAA.list[csRange],RunChiSq,mc.cores=Cores)
    } else {
      ChiSqTabAA.list <- parallel::mclapply(ConTabAA.list[csRange],RunChiSq_c,mc.cores=Cores)
    }


    # build data frame for 2x2 tables
    Final_binned.list <- lapply(ChiSqTabAA.list,"[[",1)
    ccdat.list <- lapply(Final_binned.list,TableMaker)
    OR.list <- lapply(ccdat.list, cci.pval.list) #OR

  } else {

    # Filter exon alignment matrix for single allele
    TabAA <- AlignmentFilter(ExonAlign,Alleles2F,Locus)
    FlagAA.list <- as.list(rep(TRUE,length(TabAA[6:length(TabAA)])))
    names(FlagAA.list) <- names(TabAA[6:length(TabAA)])

    ChiSqTabAA.list <- NA
    OR.list <- NA
    Final_binned.list <- NA

  }

  ####################################################### Build Output List

  A.tmp <- list()

  ## AAlog_out - Positions with insufficient variation
  csRange <- which(FlagAA.list==T)
  FlagAA.fail <- rownames(do.call(rbind,FlagAA.list[csRange]))
  AAlog.out <- cbind(rep(Locus,length(FlagAA.fail)),
                     FlagAA.fail,
                     rep("Insufficient variation at position.",length(FlagAA.fail)))
  colnames(AAlog.out) <- c("Locus","Position","Comment")
  rownames(AAlog.out) <- NULL
  A.tmp[['log']] <- AAlog.out

  ## AminoAcid.binned_out
  if(length(ChiSqTabAA.list)>1) {
    binned.list <- lapply(ChiSqTabAA.list,"[[",2)
    binned.list <- binned.list[which(lapply(binned.list,is.logical)==F)]
    binned.out <- do.call(rbind,binned.list)
    if(!is.null(nrow(binned.out))) {
      binned.out <- cbind(rep(Locus,nrow(binned.out)),
                          rep(names(binned.list),as.numeric(lapply(binned.list,nrow))),
                          rownames(binned.out),
                          binned.out)
      colnames(binned.out) <- c("Locus","Position","Residue","Group.0","Group.1")
      rownames(binned.out) <- NULL
      A.tmp[['binned']] <- binned.out
    } else {
       binned.out <- cbind(Locus,'Nothing.binned',NA,NA,NA)
       colnames(binned.out) <- c("Locus","Position","Residue","Group.0","Group.1")
       rownames(binned.out) <- NULL
       A.tmp[['binned']] <- binned.out
    }
  } else{
    binned.out <- cbind(Locus,'Nothing.binned',NA,NA,NA)
    colnames(binned.out) <- c("Locus","Position","Residue","Group.0","Group.1")
    rownames(binned.out) <- NULL
    A.tmp[['binned']] <- binned.out
  }

  ## overall.chisq_out
  if(length(ChiSqTabAA.list)>1) {
    ChiSq.list <- lapply(ChiSqTabAA.list,"[[",3)
    ChiSq.out <- do.call(rbind,ChiSq.list)
    if(!is.null(ChiSq.out)){
      ChiSq.out <- cbind(rep(Locus,nrow(ChiSq.out)),
                         rownames(ChiSq.out),
                         ChiSq.out)
      colnames(ChiSq.out) <- c("Locus","Position","X.square","df","p.value","sig")
      rownames(ChiSq.out) <- NULL
      A.tmp[['chisq']] <- ChiSq.out
    } else {
      Names <- c("Locus","Position","X.square","df","p.value","sig")
      A.tmp[['chisq']] <- Create.Null.Table(Locus,Names,nr=1)
    }
  } else {
    Names <- c("Locus","Position","X.square","df","p.value","sig")
    A.tmp[['chisq']] <- Create.Null.Table(Locus,Names,nr=1)
  }

  ## ORtable_out
  if(length(OR.list)>1) {
    OR.out <- do.call(rbind,OR.list)
    if(!is.null(OR.out)) {
      OR.out <- cbind(rep(Locus,nrow(OR.out)),
                      rep(names(OR.list),as.numeric(lapply(OR.list,nrow))),
                      rownames(OR.out),
                      OR.out)
      colnames(OR.out) <- c("Locus","Position","Residue","OR","CI.lower","CI.upper","p.value","sig")
      rownames(OR.out) <- NULL
      A.tmp[['OR']] <- OR.out
    } else {
      Names <- c("Locus","Position","Residue","OR","CI.lower","CI.upper","p.value","sig")
      A.tmp[['OR']] <- Create.Null.Table(Locus,Names,nr=1)
    }
  } else {
    Names <- c("Locus","Position","Residue","OR","CI.lower","CI.upper","p.value","sig")
    A.tmp[['OR']] <- Create.Null.Table(Locus,Names,nr=1)
  }

  ## Final_binned_out (Final Table)
  if( length(Final_binned.list)>1 ) {
    Final_binned.out <- do.call(rbind,Final_binned.list)
    if(!is.null(Final_binned.out)) {
      Final_binned.out <- cbind(rep(Locus,nrow(Final_binned.out)),
                                rep(names(Final_binned.list),as.numeric(lapply(Final_binned.list,nrow))),
                                rownames(Final_binned.out),
                                Final_binned.out)
      colnames(Final_binned.out) <- c("Locus","Position","Residue","Group.0","Group.1")
      rownames(Final_binned.out) <- NULL
      A.tmp[['table']] <- Final_binned.out
    } else {
      Final_binned.out <- NULL
      Names <- c("Locus","Position","Residue","Group.0","Group.1")
      A.tmp[['table']] <- Create.Null.Table(Locus,Names,nr=1)
    }
  } else {
    Final_binned.out <- NULL
    Names <- c("Locus","Position","Residue","Group.0","Group.1")
    A.tmp[['table']] <- Create.Null.Table(Locus,Names,nr=1)
  }

  ## AminoAcid.freq_out
  if( !is.null(Final_binned.out) ) {
    Final_binned.out[,'Group.0'] <- round(as.numeric(Final_binned.out[,'Group.0']) / sum(as.numeric(Final_binned.out),na.rm=T),digits=5)
    Final_binned.out[,'Group.1'] <- round(as.numeric(Final_binned.out[,'Group.1']) / sum(as.numeric(Final_binned.out),na.rm=T),digits=5)
    A.tmp[['freq']] <- Final_binned.out
  } else {
    Names <- c("Locus","Position","Residue","Group.0","Group.1")
    A.tmp[['freq']] <- Create.Null.Table(Locus,Names,nr=1)
  }

  return(A.tmp)

}
