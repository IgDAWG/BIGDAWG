#' Case-Control Odds ratio calculation and graphing
#'
#' cci function port epicalc version 2.15.1.0 (Virasakdi Chongsuvivatwong, 2012)
#' @param caseexp Number of cases exposed
#' @param controlex Number of controls exposed
#' @param casenonex Number of cases not exosed
#' @param controlnonex Number of controls not exposed
#' @param cctable A 2-by-2 table. If specified, will supercede the outcome and exposure variables
#' @param graph If TRUE (default), produces an odds ratio plot
#' @param design Specification for graph; can be "case control","case-control", "cohort" or "prospective"
#' @param main main title of the graph
#' @param xlab label on X axis
#' @param ylab label on Y axis
#' @param xaxis two categories of exposure in graph
#' @param yaxis two categories of outcome in graph
#' @param alpha level of significance
#' @param fisher.or whether odds ratio should be computed by the exact method
#' @param exact.ci.or whether confidence limite of the odds ratio should be computed by the exact method
#' @param decimal number of decimal places displayed
#' @note This function is for internal BIGDAWG use only.
cci <- function (caseexp, controlex, casenonex, controlnonex, cctable = NULL, graph = TRUE, design = "cohort", main, xlab, ylab, xaxis, yaxis, alpha = 0.05, fisher.or = FALSE, exact.ci.or = TRUE, decimal = 2) {

  if (is.null(cctable)) {
    frame <- cbind(Outcome <- c(1, 0, 1, 0), Exposure <- c(1, 1, 0, 0), Freq <- c(caseexp, controlex, casenonex, controlnonex))
    Exposure <- factor(Exposure)
    expgrouplab <- c("Non-exposed", "Exposed")
    levels(Exposure) <- expgrouplab
    Outcome <- factor(Outcome)
    outcomelab <- c("Negative", "Positive")
    levels(Outcome) <- outcomelab
    table1 <- xtabs(Freq ~ Outcome + Exposure, data = frame)
  }
  else {
    table1 <- as.table(get("cctable"))
  }
  fisher <- fisher.test(table1)
  caseexp <- table1[2, 2]
  controlex <- table1[1, 2]
  casenonex <- table1[2, 1]
  controlnonex <- table1[1, 1]
  se.ln.or <- sqrt(1/caseexp + 1/controlex + 1/casenonex +
                     1/controlnonex)
  if (!fisher.or) {
    or <- caseexp/controlex/casenonex * controlnonex
    p.value <- chisq.test(table1, correct = FALSE)$p.value
  }
  else {
    or <- fisher$estimate
    p.value <- fisher$p.value
  }
  if (exact.ci.or) {
    ci.or <- as.numeric(fisher$conf.int)
  }
  else {
    ci.or <- or * exp(c(-1, 1) * qnorm(1 - alpha/2) * se.ln.or)
  }
  if (graph == TRUE) {
    caseexp <- table1[2, 2]
    controlex <- table1[1, 2]
    casenonex <- table1[2, 1]
    controlnonex <- table1[1, 1]
    if (!any(c(caseexp, controlex, casenonex, controlnonex) <
             5)) {
      if (design == "prospective" || design == "cohort" ||
          design == "cross-sectional") {

        if (missing(main))
          main <- "Odds ratio from prospective/X-sectional study"
        if (missing(xlab))
          xlab <- ""
        if (missing(ylab))
          ylab <- paste("Odds of being", ifelse(missing(yaxis),
                                                "a case", yaxis[2]))
        if (missing(xaxis))
          xaxis <- c("non-exposed", "exposed")
        axis(1, at = c(0, 1), labels = xaxis)
      }
      else {

        if (missing(main))
          main <- "Odds ratio from case control study"
        if (missing(ylab))
          ylab <- "Outcome category"
        if (missing(xlab))
          xlab <- ""
        if (missing(yaxis))
          yaxis <- c("Control", "Case")
        axis(2, at = c(0, 1), labels = yaxis, las = 2)
        mtext(paste("Odds of ", ifelse(xlab == "", "being exposed",
                                       paste("exposure being", xaxis[2]))), side = 1,
              line = ifelse(xlab == "", 2.5, 1.8))
      }
      title(main = main, xlab = xlab, ylab = ylab)
    }
  }
  if (!fisher.or) {
    results <- list(or.method = "Asymptotic", or = or, se.ln.or = se.ln.or,
                    alpha = alpha, exact.ci.or = exact.ci.or, ci.or = ci.or,
                    table = table1, decimal = decimal)
  }
  else {
    results <- list(or.method = "Fisher's", or = or, alpha = alpha,
                    exact.ci.or = exact.ci.or, ci.or = ci.or, table = table1,
                    decimal = decimal)
  }
  class(results) <- c("cci", "cc")
  return(results)
}

#' Creation of a 2x2 table using the indicated orientation.
#'
#' make2x2 function port epicalc version 2.15.1.0 (Virasakdi Chongsuvivatwong, 2012)
#' @param caseexp Number of cases exposed
#' @param controlex Number of controls exposed
#' @param casenonex Number of cases not exosed
#' @param controlnonex Number of controls not exposed
#' @note This function is for internal BIGDAWG use only.
make2x2 <- function (caseexp, controlex, casenonex, controlnonex)  {

  table1 <- c(controlnonex, casenonex, controlex, caseexp)
  dim(table1) <- c(2, 2)
  rownames(table1) <- c("Non-diseased", "Diseased")
  colnames(table1) <- c("Non-exposed", "Exposed")
  attr(attr(table1, "dimnames"), "names") <- c("Outcome", "Exposure")
  table1

}

#' Table Maker
#'
#' Table construction of per haplotype for odds ratio, confidence intervals, and pvalues
#' @param x Contingency table with binned rare cells.
#' @note This function is for internal BIGDAWG use only.
TableMaker <- function(x) {
  grp1_sum <- sum(x[,'Group.1'])
  grp0_sum <- sum(x[,'Group.0'])
  grp1_exp <- x[,'Group.1']
  grp0_exp <- x[,'Group.0']
  grp1_nexp <- grp1_sum - grp1_exp
  grp0_nexp <- grp0_sum - grp0_exp
  cclist <- cbind(grp1_exp, grp0_exp, grp1_nexp, grp0_nexp)
  tmp <- as.data.frame(t(cclist))
  names(tmp) <- row.names(x)
  return(tmp)
}

#' Case Control Odds Ratio Calculation from Epicalc
#'
#' Calculates odds ratio and pvalues from 2x2 table
#' @param x List of 2x2 matrices for calculation, output of TableMaker.
#' @note This function is for internal BIGDAWG use only.
cci.pval <- function(x) {
  tmp <- list()
  caseEx <- x[1]
  controlEx <- x[2]
  caseNonEx <- x[3]
  controlNonEx <- x[4]
  table1 <- make2x2(caseEx, controlEx, caseNonEx, controlNonEx)
  tmp1 <- cci(cctable=table1, design = "case-control", graph = FALSE)
  tmp[['OR']] <- round(tmp1$or,digits=2)
  tmp[['CI.L']] <- round(tmp1$ci.or[1],digits=2)
  tmp[['CI.U']] <- round(tmp1$ci.or[2],digits=2)
  tmp[['p.value']] <-  format.pval(chisq.test(table1, correct=F)$p.value)
  tmp[['sig']] <- ifelse(chisq.test(table1, correct=F)$p.value <= 0.05,"*","NS")
  return(tmp)
}

#' Case Control Odds Ratio Calculation from Epicalc list variation
#'
#' Variation of the cci.pvalue function
#' @param x List of 2x2 matrices to apply the cci.pvalue function. List output of TableMaker.
#' @note This function is for internal BIGDAWG use only.
cci.pval.list <- function(x) {
  tmp <- lapply(x, cci.pval)
  tmp <- do.call(rbind,tmp)
  colnames(tmp) <- c("OR","CI.lower","CI.upper","p.value","sig")
  return(tmp)
}

#' Contextual Binning Chi-squared Contingency Table Test
#'
#' Calculates chi-squared contingency table tests and bins rare cells at 20% capture rate.
#' @param x Contingency table.
#' @note This function is for internal BIGDAWG use only.
RunChiSq_c <- function(x) {

  ### get expected values for cells
  ExpCnts <- chisq.test(as.matrix(x))$expected

  # Order Counts
  getOrder <- order(ExpCnts[,1],ExpCnts[,2],decreasing=T)
  ExpCnts <- ExpCnts[getOrder,]
  x.sub <- x[getOrder,]

  ### pull out cells that don't need binning, bin remaining
  #unbinned
  OK.rows <- as.numeric(which(apply(ExpCnts,min,MARGIN=1)>=5))
  if(length(OK.rows)>0) {
    if(length(OK.rows)>=2) {
      unbinned <- x.sub[OK.rows,]
    } else {
      unbinned <- do.call(cbind,as.list(x.sub[OK.rows,]))
      rownames(unbinned) <- rownames(x.sub)[OK.rows]
    }
  } else {
    unbinned <- NULL
  }

  #binned based on 20% threshold
  Rare.rows <- as.numeric(which(apply(ExpCnts,min,MARGIN=1)<5))
  minRows <- floor(0.2*nrow(x.sub)) #min number of rows to bin to be under 20% threshold
  if(length(Rare.rows)>=2 && length(Rare.rows)<=minRows) {

    # Contextual Min Rows Equal to Strict Rows to Bin
    binned <- x.sub[Rare.rows,]
    New.df <- rbind(unbinned,colSums(x.sub[Rare.rows,]))
    rownames(New.df)[nrow(New.df)] <- "binned"

  } else if(length(Rare.rows)>=2 && length(Rare.rows)>minRows) {

    # Subsetting Strict Binning to Retain Rows
    Rare.rows.in <- sort(sort(Rare.rows,decreasing=T)[1:minRows]) # to be binned
    Rare.rows.out <- setdiff(Rare.rows,Rare.rows.in) # to be rescued from binning

    if(length(Rare.rows.out)>1) {
      # Rebin in.rows that overlap with the exp counts of out.rows
      out.cnts <- apply(ExpCnts[Rare.rows.out,],MARGIN=1,paste,collapse=":")
    } else {
      out.cnts <- paste(ExpCnts[Rare.rows.out,],collapse=":")
    }


    if(length(Rare.rows.in)>1) {
      in.cnts <- apply(ExpCnts[Rare.rows.in,],MARGIN=1,paste,collapse=":")
      in.cnts.rev <- apply(ExpCnts[Rare.rows.in,c(2,1)],MARGIN=1,paste,collapse=":")
    } else {
      in.cnts <- paste(ExpCnts[Rare.rows.in,],collapse=":")
      in.cnts.rev <- paste(ExpCnts[Rare.rows.in,c(2,1)],collapse=":")
    }

    Rare.rows.final <- c(Rare.rows.in,
                         Rare.rows.out[out.cnts %in% in.cnts],
                         Rare.rows.out[out.cnts %in% in.cnts.rev])

    binned <- x.sub[Rare.rows.final,]
    New.df <- rbind(unbinned,colSums(x.sub[Rare.rows.final,]))
    rownames(New.df)[nrow(New.df)] <- "binned"

  } else {

    binned <- cbind(NA,NA)
    colnames(binned) <- c("Group.0","Group.1")
    New.df <- ExpCnts

  }

  if(nrow(New.df)>1) {

    # flag if final matrix fails Cochran's rule of thumb (more than 20% of exp cells are less than 5)
    # True = OK ; False = Not good for Chi Square
    ExpCnts <- chisq.test(New.df)$expected
    if(sum(ExpCnts<5)==0){
      flag <- FALSE
    } else if( sum(ExpCnts<5)/sum(ExpCnts>=0)<=0.2 && sum(ExpCnts>=1)>length(ExpCnts) ){
      flag <- FALSE
    } else {
      flag <- TRUE
    }

    ## chi square test on binned data
    df.chisq <- chisq.test(New.df)
    Sig <- if(df.chisq$p.value > 0.05) { "NS" } else { "*" }


    ## show results of overall chi-square analysis
    tmp.chisq <- data.frame(cbind(round(df.chisq$statistic,digits=4),
                                  df.chisq$parameter,
                                  format.pval(df.chisq$p.value),
                                  Sig))
    colnames(tmp.chisq) <- c("X.square", "df", "p.value", "sig")

    chisq.out <- list(Matrix = New.df,
                      Binned = binned,
                      Test = tmp.chisq,
                      Flag = flag)

    return(chisq.out)

  } else {

    flag <- TRUE
    tmp.chisq <- data.frame(rbind(rep("NCalc",4)))
    colnames(tmp.chisq) <- c("X.square", "df", "p.value", "sig")
    chisq.out <- list(Matrix = New.df,
                      Binned = binned,
                      Test = tmp.chisq,
                      Flag = flag)

    return(chisq.out)

  }

}

#' Strict Chi-squared Contingency Table Test
#'
#' Calculates chi-squared contingency table tests and bins all rare cells.
#' @param x Contingency table.
#' @note This function is for internal BIGDAWG use only.
RunChiSq <- function(x) {

  ### get expected values for cells
  ExpCnts <- chisq.test(as.matrix(x))$expected

  ## pull out cells that don't need binning, bin remaining
  #unbinned
  OK.rows <- as.numeric(which(apply(ExpCnts,min,MARGIN=1)>=5))
  if(length(OK.rows)>0) {
    if(length(OK.rows)>=2) {
      unbinned <- x[OK.rows,]
    } else {
      unbinned <- do.call(cbind,as.list(x[OK.rows,]))
      rownames(unbinned) <- rownames(x)[OK.rows]
    }
  } else {
    unbinned <- NULL
  }

  #binned
  Rare.rows <- as.numeric(which(apply(ExpCnts,min,MARGIN=1)<5))
  if(length(Rare.rows)>=2) {
    binned <- x[Rare.rows,]
    New.df <- rbind(unbinned,colSums(x[Rare.rows,]))
    rownames(New.df)[nrow(New.df)] <- "binned"
  } else {
    binned <- cbind(NA,NA)
    colnames(binned) <- c("Group.0","Group.1")
    New.df <- x
  }

  if(nrow(New.df)>1) {

    # flag if final matrix fails Cochran's rule of thumb (more than 20% of exp cells are less than 5)
    # True = OK ; False = Not good for Chi Square
    ExpCnts <- chisq.test(New.df)$expected
    if(sum(ExpCnts<5)==0){
      flag <- FALSE
    } else if( sum(ExpCnts<5)/sum(ExpCnts>=0)<=0.2 && sum(ExpCnts>=1)>length(ExpCnts) ){
      flag <- FALSE
    } else {
      flag <- TRUE
    }

    ## chi square test on binned data
    df.chisq <- chisq.test(New.df)
    Sig <- if(df.chisq$p.value > 0.05) { "NS" } else { "*" }


    ## show results of overall chi-square analysis
    tmp.chisq <- data.frame(cbind(round(df.chisq$statistic,digits=4),
                                  df.chisq$parameter,
                                  format.pval(df.chisq$p.value),
                                  Sig))
    colnames(tmp.chisq) <- c("X.square", "df", "p.value", "sig")

    chisq.out <- list(Matrix = New.df,
                      Binned = binned,
                      Test = tmp.chisq,
                      Flag = flag)

    return(chisq.out)

  } else {

    flag <- TRUE
    tmp.chisq <- data.frame(rbind(rep("NCalc",4)))
    colnames(tmp.chisq) <- c("X.square", "df", "p.value", "sig")
    chisq.out <- list(Matrix = New.df,
                      Binned = binned,
                      Test = tmp.chisq,
                      Flag = flag)

  }

}