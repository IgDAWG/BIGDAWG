#' File Fetcher
#'
#' Download Protein Alignment and Accessory Files
#' @param Loci HLA Loci to be fetched. Limited Loci available.
#' @note This function is for internal BIGDAWG use only.
GetFiles <- function(Loci) {
  #downloads *_prot.txt alignment files
  #downloads hla_nom_p.txt file

  # Get P-Groups Files
  download.file("ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/wmda/hla_nom_p.txt",destfile="hla_nom_p.txt",method="libcurl")

  # Get Release Version
  download.file("ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/release_version.txt",destfile="release_version.txt",method="libcurl")
  Release <- read.table("release_version.txt",comment.char="",sep="\t")
  Release <- apply(Release,MARGIN=1,FUN= function(x) gsub(": ",":",x))
  RD <- unlist(strsplit(Release[2],split=":"))[2]
  RV <- unlist(strsplit(Release[3],split=":"))[2]
  write.table(c(RD,RV),file="Release.txt",quote=F,col.names=F,row.names=F)
  file.remove("release_version.txt")

  # Get Locus Based Alignments
  for(i in 1:length(Loci)) {
    Locus <- Loci[i]
    FileName <- paste(Locus,"_prot.txt",sep="")
    URL <- paste("ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/",FileName,sep="")
    download.file(URL,destfile = FileName,method="libcurl")
  }

}

#' HLA P group File Formatter
#'
#' Format the hla_nom_p.txt read table object for a specific locus.
#' @param x P group object from read.table command.
#' @param Locus Locus to be filtered on.
#' @note This function is for internal BIGDAWG use only.
PgrpFormat <- function(x,Locus) {

  # Identify for Locus ... change necessary if DRB
  x.sub <- x[which(x[,1]==Locus),]
  rownames(x.sub) <- NULL
  x.sub[,2] <- sapply(x.sub[,2],function(i) paste(paste(Locus,"*",unlist(strsplit(i,"/")),sep=""),collapse="/"))
  colnames(x.sub) <- c("Locus","Allele","P.Group")

  #Expand
  x.list <- list()
  for(i in 1:nrow(x.sub)) {
    if(grepl("/",x.sub[i,'Allele'],fixed=T)) {
      tmp <- unlist(strsplit(x.sub[i,'Allele'],"/"))
      tmp <- cbind(rep(Locus,length(tmp)),tmp,rep(x.sub[i,'P.Group'],length(tmp)))
      colnames(tmp) <- colnames(x.sub)
      x.list[[i]] <- tmp
    } else {
      x.list[[i]] <- x.sub[i,]
    }
  }
  x.list <- do.call(rbind,x.list)
  x.list <- x.list[order(x.list[,'Allele']),]
  rownames(x.list) <- NULL
  colnames(x.list) <- c("Locus","Allele","P.Group")
  return(x.list)
}

#' HLA P group Finder
#'
#' Identify P group for a given allele if exists.
#' @param x Allele of interest.
#' @param y Formatted P groups.
#' @note This function is for internal BIGDAWG use only.
PgrpExtract <- function(x,y) {
  getRow <- grep(x,y[,'Allele'],fixed=T)
  if(length(getRow)>=1) {
    if(length(getRow)>1) { getRow <- getRow[which(sapply(as.character(y[getRow,'Allele']),nchar)==nchar(x))] }
    return(as.character(y[getRow,'P.Group']))
  } else { return("") }
}

#' Protein Exon Alignment Formatter
#'
#' Dynamically creates an alignmnet of Allele exons for Analysis.
#' @param Locus Locus alignment to be formatted.
#' @param RefTab Reference exon protein information for alignment formatting.
#' @note This function is for internal BIGDAWG use only.
ExonPtnAlign.Create <- function(Locus,RefTab) {

  #########################################################################
  # Need to remove if DRB split into single locus files
  if(grepl("DRB",Locus)) { Locus.get <- "DRB" } else { Locus.get <- Locus }
  #########################################################################

  AlignMatrix <- NULL; rm(AlignMatrix)

  #Read in P-Groups
  Pgrps <- read.table("hla_nom_p.txt",fill=T,header=F,sep=";",stringsAsFactors=F,strip.white=T,colClasses="character")
  Pgrps[,1] <- gsub("\\*","",Pgrps[,1])
  Pgrps <- PgrpFormat(Pgrps,Locus)

  #Read in Alignment
  Name <- paste(Locus.get,"_prot.txt",sep="")
  Align <- read.table(Name,fill=T,header=F,sep="\t",stringsAsFactors=F,strip.white=T,colClasses="character")

  #Trim
  Align <- as.matrix(Align[-nrow(Align),]) #Remove Footer

  #Begin Formatting
  Align <- as.matrix(Align[-grep("\\|",Align[,1]),1]) #Remove Pipes
  Align[,1] <- sapply(Align[,1],FUN=sub,pattern=" ",replacement="~")
  Align[,1] <- sapply(Align[,1],FUN=gsub,pattern=" ",replacement="")
  Align <- strsplit(Align[,1],"~")
  Align <- as.matrix(do.call(rbind,Align))

  #Adjust rows to blank where Sequence column == Allele Name
  Align[which(Align[,1]==Align[,2]),2] <- ""

  #Find start of repeating blocks
  Start <- which(Align[,1]=="Prot") + 1
  End <- c(Start[2:length(Start)] - 1,nrow(Align))

  #Rearrange Alignment, build Alignment Block file
  Block.size <- End[1]-Start[1] + 1
  Block <- mat.or.vec(nr=Block.size,nc=length(Start)+1)
  for(i in 1:length(Start) ) {
    if(i==1) {
      Block[,1:2] <- Align[Start[i]:End[i],1:2]
    } else {
      putOrder <- match(Align[Start[i]:End[i],1],Block[,1])
      Block[putOrder,i+1] <- Align[Start[i]:End[i],2]
    }
  }; rm(i)
  Block[,2:ncol(Block)] <- apply(Block[,2:ncol(Block)],MARGIN=c(1,2),FUN=gsub,pattern=0,replacement="")

  #Paste Sequence into Single Column -- Fill in gaps with * to make char lengths even
  Block <- cbind(Block[,1],apply(Block[,2:ncol(Block)],MARGIN=1,paste,collapse=""))
  Block.len <- max(as.numeric(sapply(Block[,2],FUN=nchar)))
  for( i in 1:nrow(Block) ) {
    Block.miss <- Block.len - nchar(Block[i,2])
    if( Block.miss > 0 ) {
      Block[i,2] <- paste(Block[i,2], paste(rep("*",Block.miss-1), collapse=""))
    }
  }; rm(i)


  #Split Allele name into separate Locus and Allele, Send Back to Align object
  AlignAlleles <- do.call(rbind,strsplit(Block[,1],"[*]"))
  AlignAlleles <- cbind(AlignAlleles,apply(AlignAlleles,MARGIN=c(1,2),FUN=GetField,Res=2)[,2])
  rownames(AlignAlleles) <- NULL

  Align <- cbind(AlignAlleles,Block)
  colnames(Align) <- c("Locus","Allele","Trimmed","FullName","Sequence")

  #Define Reference
  RefSeq <- Align[1,]

  #Ensure Locus Specific Rows - Add Reference
  Align <- Align[which(Align[,'Locus']==Locus),]
  Align <- rbind(RefSeq,Align)
  Align[1,1:4]  <- "RefSeq"
  rownames(Align) <- NULL

  #Get Reference Exon Map
  RefExon <- RefTab[which(RefTab[,'Locus']==Locus),'Reference.Peptide']
  RefStart <- as.numeric(RefTab[which(RefTab[,'Locus']==Locus),'Reference.Start'])
  RefAllele <- RefTab[which(RefTab[,'Locus']==Locus),'Reference.Allele']

  #Find Exon Specific Overlap Based on RefExon
  Align.seq <- Align[1,'Sequence'] # map character positions
  Align.map <- 1:nchar(Align.seq)
  Align.map.rm <- unlist(gregexpr("\\.",Align.seq))
  Align.map.sub <- setdiff(Align.map,Align.map.rm) # remap based on removed characters
  Align.seq.sub <- gsub("\\.","",Align.seq)

  Start.sub <- unlist(gregexpr(RefExon,Align.seq.sub))
  End.sub <- Start.sub + nchar(RefExon) - 1
  getMap <- Align.map.sub[Start.sub:End.sub]
  Exon.start <- min(getMap)
  Exon.end <- max(getMap)

  #Split Sub Alignment into composite elements and Extract relavent positions
  Align.split <- strsplit(Align[,'Sequence'],"")
  Align.split <- do.call(rbind,Align.split)
  Align.split <- Align.split[,Exon.start:Exon.end]

  AlignMatrix <- cbind(Align[,1:4],Align.split)
  colnames(AlignMatrix) <- c("Locus","Allele","Trimmed","FullName",paste("Position",seq(RefStart,ncol(AlignMatrix[,5:ncol(AlignMatrix)]) + RefStart - 1),sep="."))
  rownames(AlignMatrix) <- NULL

  #Propagate Consensus Positions
  for(i in 5:ncol(AlignMatrix)) {
    x <- AlignMatrix[,i]
    x[which(x=="-")] <- x[1]
    AlignMatrix[,i] <- x
  }

  #Remove Reference
  AlignMatrix <- AlignMatrix[-1,]

  #Assign P groups
  AlignMatrix <- cbind(AlignMatrix,sapply(AlignMatrix[,'FullName'],PgrpExtract,y=Pgrps))
  colnames(AlignMatrix)[ncol(AlignMatrix)] <- "P group"

  #Add Absent Allele (Absence due to lack of allele and not lack of typing information)
  AlignMatrix <- rbind(c(Locus,"00:00:00:00","00:00",paste(Locus,"*00:00:00:00",sep=""),rep("^",ncol(AlignMatrix)-4)),
                       AlignMatrix)

  #Tally Unknowns as separate column
  AlignMatrix <- cbind(AlignMatrix,
                       apply(AlignMatrix[,5:(ncol(AlignMatrix)-1)],MARGIN=1,FUN=function(x) length(which(unlist(grep("*",x,fixed=T))>0))) )
  colnames(AlignMatrix)[ncol(AlignMatrix)] <- "Unknowns"

  #Tally Null Positions as separate column
  AlignMatrix <- cbind(AlignMatrix,
                       apply(AlignMatrix[,5:(ncol(AlignMatrix)-2)],MARGIN=1,FUN=function(x) length(which(unlist(grep("-",x,fixed=T))>0))) )
  colnames(AlignMatrix)[ncol(AlignMatrix)] <- "NullPositions"

  #Tally InDels
  AlignMatrix <- cbind(AlignMatrix,
                       apply(AlignMatrix[,5:(ncol(AlignMatrix)-3)],MARGIN=1,FUN=function(x) length(which(unlist(grep(".",x,fixed=T))>0))) )
  colnames(AlignMatrix)[ncol(AlignMatrix)] <- "InDels"

  rownames(AlignMatrix) <- NULL
  FileName <- paste("ExonPtnAlign_",Locus,".obj",sep="")

  save(AlignMatrix,file=FileName)

}

#' Alignment Object Creator
#'
#' Synthesize Object for Exon Protein Alignments.
#' @param Loci Loci to be bundled.
#' @param Release IMGT/HLA database release version.
#' @param RefTab Data of reference exons used for protein alignment creation.
#' @note This function is for internal BIGDAWG use only.
AlignObj.Create <- function(Loci,Release,RefTab) {

  AlignMatrix <- NULL; rm(AlignMatrix)

  ExonPtnList <- list()
  for(i in 1:length(Loci)) {
    Locus <- Loci[i]
    FileName <- paste("ExonPtnAlign_",Locus,".obj",sep="")
    load(FileName) #Loads AlignMatrix
    ExonPtnList[[Locus]] <- AlignMatrix
  }

  ExonPtnList[['Release.Version']] <- as.character(Release[2,])
  ExonPtnList[['Release.Date']] <- as.character(Release[1,])
  ExonPtnList[['RefExons']] <- RefTab
  save(ExonPtnList,file="ExonPtnAlign.obj")

}

#' Updated Alignment Object Creator
#'
#' Synthesize Object for Exon Protein Alignments.
#' @param Loci Loci to be bundled.
#' @param Release IMGT/HLA database release version.
#' @param RefTab Data of reference exons used for protein alignment creation.
#' @note This function is for internal BIGDAWG use only.
AlignObj.Update <- function(Loci,Release,RefTab) {

  AlignMatrix <- NULL; rm(AlignMatrix)

  UpdatePtnList <- list()
  for(i in 1:length(Loci)) {
    Locus <- Loci[i]
    FileName <- paste("ExonPtnAlign_",Locus,".obj",sep="")
    load(FileName) #Loads AlignMatrix
    UpdatePtnList[[Locus]] <- AlignMatrix
  }
  UpdatePtnList[['Release']] <- Release
  UpdatePtnList[['RefExons']] <- RefTab
  save(UpdatePtnList,file="UpdatePtnAlign.RData")

}
