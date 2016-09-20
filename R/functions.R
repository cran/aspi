##################################################################################
#                                                                                #
#  Analysis of Symmetry of Parasitic Infections (ASPI)                           #
#                                                                                #
#  Version 0.2.0                                                                   #
#                                                                                #
#  Copyright (C) 2016 Matthew Thomas Wayland                                     #
#                                                                                #
#  This file is part of Analysis of Symmetry of Parasitic Infections.            #
#                                                                                #
#  Analysis of Symmetry of Parasitic Infections is free software: you can        #
#  redistribute it and/or modify it under the terms of the GNU General Public    #
#  License as published by the Free Software Foundation, either version 3 of     #
#  the License, or (at your option) any later version.                           #
#                                                                                #
#  Analysis of Symmetry of Parasitic Infections is distributed in the hope that  #
#  it will be useful, but WITHOUT ANY WARRANTY; without even the implied         #
#  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     #
#  GNU General Public License for more details.                                  #
#                                                                                #
#  You should have received a copy of the GNU General Public License along       #
#  with Analysis of Symmetry of Parasitic Infections.                            #
#  If not, see <http://www.gnu.org/licenses/>.                                   #
#                                                                                #
#                                                                                #
##################################################################################

#'@importFrom graphics abline hist plot
#'@importFrom stats binom.test p.adjust pchisq

# check data format and remove uninfected hosts
.prepareData <- function(x){
  x <- as.data.frame(x)

  if(length(x[1,])!=2)
    stop("x must contain exactly 2 columns")

  if(length(x[,1])<1)
    stop("x must contain data for at least one host")

  if(!is.numeric(x[,1]) | !is.numeric(x[,2]))
    stop("Both columns of x must be numeric")

  return(x[apply(x,1,sum)>0,])
}

#'Replicated G-tests of goodness-of-fit
#'@description Perform replicated G-tests of goodness-of-fit to assess symmetry of parasitic infections.
#'@param x a matrix or data frame with two numeric columns;
#'first column is for left-side and 2nd column for right-side.
#'Identifiers for hosts can be provided as row names.
#'@details
#'This function implements Sokal & Rohlf's (1995) G-test for the specific case of an expected 1:1 ratio
#'The function takes as its argument a matrix or data frame with two numeric columns; first column is for
#'left-side and 2nd column for right-side. Identifiers for hosts can be provided as row names. Uninfected
#' hosts (zero count for both left and right sides) are ignored. Cannot be applied to data containing zero
#' counts; use eb.test instead.

#'@return
#'A list containing two data.frames:
#'  \item{summary}{results of total, heterogeneity and pooled G-tests.
#'    Data frame has four columns: test, degrees of freedom, G-statistic and p-value.}
#'  \item{hosts}{results of individual G-tests on distribution of parasites in each host.
#'    Data frame has seven columns: Host (ID), Left (count of parasites on left side),
#'    Right (count of parasites on right side), G (G-statistic), p (p-value), BH (p-value adjusted
#'    using Benjamini and Hochberg's procedure for controlling the false discovery rate) and Holm (p-value
#'    adjusted using Holm's method).}
#'@references
#'R.R. Sokal & F.J. Rohlf (1995) Biometry. 3rd Edition. New York: W.H. Freeman and Company. 887 pp.
#'@examples
#'g.test(diplostomum_eyes_excl_lenses)
#'@export

g.test <- function(x){
  # prepare data
  x <- .prepareData(x)

  # check for zero counts
  if(sum(x==0)!=0)
    stop("Data contains zero counts!")

  calcG <- function(fL, fR, df=1){
    fSum <- fL+fR
    G <- 2*(fL*log(fL/(fSum*0.5)) + fR*log(fR/(fSum*0.5)))
    pval <- pchisq(G, df, lower.tail=F)
    return(c(G, pval))
  }

  # G for each individual host GI
  individualG <- apply(x, 1, function(r){
    return(calcG(r[1],r[2]))
  })

  Holm <- p.adjust(individualG[2,], 'holm')
  BH <- p.adjust(individualG[2,], 'BH')

  individualG <- as.data.frame(cbind(row.names(x), x, t(individualG), BH, Holm))
  names(individualG) <- c("Host", "Left", "Right", "G", "p", "BH", "Holm")

  # Total G (GT)
  GT <- sum(individualG$G)
  dfGT <- length(individualG$G)
  pGT <- pchisq(GT, dfGT, lower.tail=F)

  # Pooled G (GP)
  pooledResult <- calcG(sum(x[,1]), sum(x[,2]))
  GP <- pooledResult[1]
  pGP <- pooledResult[2]

  # Heterogeneity G (GH)
  GH <- GT - GP
  dfGH <- dfGT-1
  pGH <- pchisq(GH, dfGH, lower.tail=F)

  # prepare output
  testSummary <- as.data.frame(cbind(c("Pooled", "Heterogeneity", "Total"), c(1, dfGH, dfGT), c(GP, GH, GT), c(pGP, pGH, pGT)))
  names(testSummary) <- c("Test", "df", "G", "p")
  testSummary$df <- as.numeric(as.character(testSummary$df))
  testSummary$G <- as.numeric(as.character(testSummary$G))
  testSummary$p <- as.numeric(as.character(testSummary$p))

  output <- list(testSummary, individualG)
  names(output) <- c("summary", "hosts")
  return(output)

}

#'Exact binomial tests
#'@description Assess symmetry of parasitic infections by performing exact binomial tests on pooled data and
#'individual hosts.
#'@param x a matrix or data frame with two numeric columns;
#'first column is for left-side and 2nd column for right-side.
#'Identifiers for hosts can be provided as row names.
#'@details
#'This function performs a binomial exact tests with the null hypothesis of a 1:1 ratio. It takes
#' as its argument a matrix or data frame with two numeric columns; first column is for left-side
#' and 2nd column for right-side. Identifiers for hosts can be provided as row names. Uninfected
#' hosts (zero count for both left and right sides) are ignored.

#'@return
#'It returns a list containing two elements:
#'  \item{pooled}{p-value for pooled binomial exact test (null hypothesis: the ratio of the total
#'  number of parasites from each side doesn't differ from 1:1).}
#'  \item{hosts}{data.frame of results of binomial exact tests performed on the distribution of
#'  parasites in each host.}
#'
#'@examples
#'eb.test(diplostomum_lenses)
#'@export

eb.test <- function(x){

  # prepare data
  x <- .prepareData(x)

  binomTestP <- function(fL, fR){
    return(binom.test(fL, fL+fR, p=0.5, alternative='two.sided')$p.value)
  }

  # Binomial tests for individual hosts
  individualP <- apply(x, 1, function(r){
    return(binomTestP(r[1],r[2]))
  })

  Holm <- p.adjust(individualP, 'holm')
  BH <- p.adjust(individualP, 'BH')

  individualP <- as.data.frame(cbind(row.names(x), x, individualP, BH, Holm))
  names(individualP) <- c("Host", "Left", "Right", "p", "BH", "Holm")

  # Binomial tests for pooled hosts
  fLSum <- sum(x[,1])
  fRSum <- sum(x[,2])

  pooledP <- binomTestP(fLSum, fRSum)

  output <- list(pooledP, individualP)
  names(output) <- c("pooled", "hosts")
  return(output)

}

#'Plot histogram
#'@description Creates a histogram showing distribution of fold differences
#'in abundance of parasites between left and right sides of host.
#'@param x a matrix or data frame with two numeric columns;
#'first column is for left-side and 2nd column for right-side.
#'Identifiers for hosts can be provided as row names.
#'@param nBreaks number of cells for the histogram. A suggestion
#'only; breakpoints will be set to pretty values.
#'@param ... optional further arguments and graphical parameters
#'passed to plot.
#'
#'@details
#' plot.Histogram creates a histogram showing distribution of fold differences
#' in abundance of parasites between left and right sides. For each infected host the
#' number of parasites on the right side is divided by the number of parasites
#' on the left side, and the result binary log transformed. The log2 ratio will
#' be negative if there are more parasites on the left than right and positive
#' if there are more parasites on the right than left. A log2 ratio of one
#' corresponds to a one-fold difference, i.e. double the number of parasites.
#' Perfect symmetry is a log2 ratio of zero.
#'
#'@examples
#'plotHistogram(diplostomum_eyes_excl_lenses)
#'plotHistogram(diplostomum_eyes_excl_lenses,nBreaks=20,
#'main="Diplostomum metacercariae in eyes of ruffe")
#'@export
plotHistogram <- function(x, nBreaks=10, ...){
  # prepare data
  x <- .prepareData(x)

  if(sum(x==0)!=0)
    stop("Data contains zero counts!")

  logRatios <- log2(x[,2]/x[,1])
  bins <- pretty(c(min(logRatios), max(logRatios)), n=as.numeric(nBreaks)+1)
  #bins <- pretty(c(min(logRatios()), max(logRatios())), n=as.numeric(input$nBreaks)+1)
  histColours <- c(rep("#4DAF4A", sum(bins<0)), rep("#984EA3", sum(bins>0)))
  hist(logRatios, col=histColours, breaks=bins, xlab=expression(paste("log"[2],"(right/left)")), ...)
}

#'Volcano plot
#'@description Produces scatterplot of statistical significance vs fold difference in parasite abundance
#'between left and right.
#'@param x a matrix or data frame with two numeric columns; first column is for left-side and 2nd column
#'for right-side. Identifiers for hosts can be provided as row names.
#'@param test if set to "G" (default) a G-test is performed; otherwise an exact binomial test is performed.
#'@param pAdj method for correcting p-values for multiple comparisons. If set to "BH" (default), Benjamini
#'& Hochberg's procedure is used to control the false discovery rate (FDR); otherwise Holm's methos is used
#'to control the familywise error rate (FWER).
#'@param sigThresh significance threshold (defaults to 0.05); p-values below this value will be called significant.
#'@param ... optional further arguments and graphical parameters passed to plot.
#'@details plot.Volcano creates a volcano plot, i.e. a scatterplot of statistical significance
#'(-log10(p-value)) vs fold difference (log2 ratio - as calculated for the histogram above)
#'in parasite abundance between left and right. Each point in the scatterplot represents the
#'parasite distribution in an individual host. A dashed horizontal line represents the user-defined
#'p-value threshold for significance. If a parasite distribution deviates significantly from symmetry
#'it is shown as a red square, otherwise as a blue circle.
#'@examples
#'plotVolcano(diplostomum_eyes_excl_lenses)
#'plotVolcano(diplostomum_eyes_excl_lenses, test="G", pAdj="BH", sigThresh=0.1,
#'main="Diplostomum metacercariae in eyes of ruffe")
#'@export
plotVolcano <- function(x, test="G", pAdj="BH", sigThresh=0.05, ...){
  # prepare data
  x <- .prepareData(x)

  if(sum(x==0)!=0)
    stop("Data contains zero counts!")

  if(test=="G"){
    testResults <- g.test(x)
  }else{
    testResults <- eb.test(x)
  }

  if(pAdj=="BH"){
    log10p <- log10(testResults$hosts$BH)
  }else{
    log10p <- log10(testResults$hosts$Holm)
  }

  logRatios <- log2(x[,2]/x[,1])
  plotColours <- ifelse(log10p < log10(as.numeric(sigThresh)), "#E41A1C", "#377EB8")
  plotSymbols <- ifelse(log10p < log10(as.numeric(sigThresh)), 15, 16)
  plot(logRatios, -log10p,  col=plotColours, pch = plotSymbols, xlab=expression(paste("log"[2],"(right/left)")), ylab=expression(paste("-log"[10],"(p-value)")), ...)
  abline(h=-log10(as.numeric(sigThresh)), col="gray10", lty=2)
}


