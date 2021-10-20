#' Alignment Filter
#'
#' Filter Protein Exon Alignment File for Specific Alleles.
#' @param Align Protein Alignment Object.
#' @param Alleles to be pulled.
#' @param Locus Locus to be filtered against.
#' @note This function is for internal BIGDAWG use only.
AlignmentFilter <- function(Align, Alleles, Locus) {

  getCols <- c(match(c("Trimmed","Unknowns","NullPositions"),
                     colnames(Align)),
               grep("Pos\\.",colnames(Align)))
  getRows <-  Align[,"Trimmed"] %in% Alleles
  Align.sub <- Align[getRows,getCols]
  Align.sub <- unique(Align.sub)

  if(!is.null(nrow(Align.sub))) {

    Alleles.S <- names(which(table(Align.sub[,'Trimmed'])==1))
    Alleles.M <- names(which(table(Align.sub[,'Trimmed'])>1))

    # Removing Duplicates at 2-Field Level
    if( length(Alleles.M > 0 ) ) {

      Align.tmp <- list()
      for( m in 1:length(Alleles.M) ) {

        Allele <- Alleles.M[m]
        Alignsub.Grp <- Align.sub[which(Align.sub[,"Trimmed"]==Allele),]

        Unknowns.Grp <- which(Alignsub.Grp[,'Unknowns']==min(Alignsub.Grp[,'Unknowns']))
        if(length(Unknowns.Grp)>1) { Unknowns.Grp <- Unknowns.Grp[1] }
        Align.tmp[[Allele]] <- Alignsub.Grp[Unknowns.Grp,]

      }
      Align.tmp <- do.call(rbind,Align.tmp)

      if( length(Alleles.S) > 0 ) {
        AlignMatrix <- rbind(Align.tmp, Align.sub[which(Align.sub[,'Trimmed'] %in% Alleles.S==T),,drop=F])
      } else {
        AlignMatrix <- Align.tmp
      }

    } else {

      AlignMatrix <- Align.sub

    }

    AlignMatrix <- cbind(rep(Locus,nrow(AlignMatrix)),AlignMatrix)
    rownames(AlignMatrix) <- NULL
    colnames(AlignMatrix)[1] <- "Locus"
    colnames(AlignMatrix)[which(colnames(AlignMatrix)=="Trimmed")] <- "Allele.2D"
    AlignMatrix <- AlignMatrix[ order(AlignMatrix[,'Allele.2D']), ]

  } else {

    AlignMatrix <- cbind(rep(Locus,nrow(Align.sub)),Align.sub)
    rownames(AlignMatrix) <- NULL
    colnames(AlignMatrix)[1] <- "Locus"
    colnames(AlignMatrix)[which(colnames(AlignMatrix)=="Trimmed")] <- "Allele.2D"
    AlignMatrix <- AlignMatrix[ order(AlignMatrix[,'Allele.2D']), ]

  }

  return(AlignMatrix)
}

#' Amino Acid Contingency Table Build
#'
#' Build Contingency Tables for Amino Acid Analysis.
#' @param x Filtered alignmnet list element.
#' @param y Phenotype groupings.
#' @note This function is for internal BIGDAWG use only.
AAtable.builder <- function(x,y) {
  #x = list element for filtered alignment
  #y = HLA_grp (case vs control)

  # Create count Grp 0 v Grp 1 (Control v Case)
  x[,2] <- gsub(" ","Null",x[,2])
  x[,2] <- gsub("\\*","Unknown",x[,2])
  x[,2] <- gsub("\\.","InDel",x[,2])

  Residues <- unique(x[,2])
  AminoAcid.df <- mat.or.vec(nr=length(Residues),2)
  colnames(AminoAcid.df) <- c("Group.0", "Group.1")
  rownames(AminoAcid.df) <- Residues

  y[,2:3] <- apply(y[,2:3],MARGIN=c(1,2),GetField,Res=2)

  # Grp 0 (Control)
  Grp0 <- y[which(y[,'grp']==0),]
  Grp0cnt <- table(c(x[match(Grp0[,2],x[,1]),2],
                     x[match(Grp0[,3],x[,1]),2]))
  PutRange <- match(rownames(AminoAcid.df),names(Grp0cnt))
  AminoAcid.df[,'Group.0'] <- Grp0cnt[PutRange]

  # Grp 1 (Case)
  Grp1 <- y[which(y[,'grp']==1),]
  Grp1cnt <- table(c(x[match(Grp1[,2],x[,1]),2],
                     x[match(Grp1[,3],x[,1]),2]))
  PutRange <- match(rownames(AminoAcid.df),names(Grp1cnt))
  AminoAcid.df[,'Group.1'] <- Grp1cnt[PutRange]

  AminoAcid.df[is.na(AminoAcid.df)] <- 0

  # drop unknowns
  rmRow <- which(row.names(AminoAcid.df)=="Unknown")
  if( length(rmRow) > 0 ) { AminoAcid.df <- AminoAcid.df[-rmRow,,drop=F] }

  return(AminoAcid.df)

}

#' Contingency Table Check
#'
#' Checks amino acid contingency table data frame to ensure required variation exists.
#' @param x contingency table.
#' @param Strict.Bin Logical specify if strict rare cell binning should be used in ChiSq test.
#' @note This function is for internal BIGDAWG use only.
AA.df.check <- function(x,Strict.Bin) {
  # Returns true if insufficient variations exists
  # RunChiSq Flag is true is sufficient variant exists

  if( nrow(x)<2 ) {
    return(FALSE)
  } else if (Strict.Bin) {
    return(!(RunChiSq(x)$Flag))
  } else {
    return(!(RunChiSq_c(x)$Flag))
  }

}

#' Contingency Table Amino Acid ChiSq Testing
#'
#' Runs ChiSq test on amino acid contingency table data frames.
#' @param x contingency table.
#' @param Strict.Bin Logical specify if strict rare cell binning should be used in ChiSq test.
#' @note This function is for internal BIGDAWG use only.
AA.df.cs <- function(x,Strict.Bin) {
  # RunChiSq on data frame

  if( nrow(x) < 2 ) {

    tmp.chisq <- data.frame(rbind(rep("NCalc",4)))
    colnames(tmp.chisq) <- c("X.square", "df", "p.value", "sig")
    row.names(tmp.chisq) <- "X-squared"
    chisq.out <- list(Matrix = x,
                      Binned = NA,
                      Test = tmp.chisq,
                      Flag = FALSE)


    return( chisq.out )

  } else if (Strict.Bin) {

    return( RunChiSq(x) )

  } else {

    return( RunChiSq_c(x) )

  }

}

#' Create Empty Table
#'
#' Creates matrix of NA for no result tables.
#' @param Locus Locus being analyzed.
#' @param Names Column names for final matrix.
#' @param nr Number of rows.
#' @note This function is for internal BIGDAWG use only.
Create.Null.Table <- function(Locus,Names,nr) {
  nc <- length(Names)
  Data <- c( Locus,rep(NA,nc-1) )
  tmp <- matrix(Data,nrow=nr,ncol=nc)
  colnames(tmp) <- Names
  rownames(tmp) <- NULL
  return(tmp)
}
