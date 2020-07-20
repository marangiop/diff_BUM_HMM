### CHANGE THE DIRECTORY TO RUN FROM THE CORRECT FOLDER!!!

### In case you don't have the dependencies:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("Biostrings")
#BiocManager::install("Rcpp", dependencies = TRUE)
#BiocManager::install("BUMHMM")

###

source("../Functions/computePvals.R")
source("../Functions/calculateLDRs.R")
source("../Functions/hmmFwbw_differential_two_betas.R")
source("../Functions/betaParamsEM.R")
source("../Functions/betaParamsMStep.R")
source("../Functions/findPatternPos.R")
source("../Functions/nuclPerm.R")
source("../Functions/scaleDOR.R")
source("../Functions/selectNuclPos.R")
source("../Functions/computeStretches.R")
source("../Functions/stabiliseVariance.R")

suppressPackageStartupMessages({ 
  library(Biostrings)
  library(SummarizedExperiment) })

#make a table: take coverage counts and drop off counts from the working directory
mergedcountswt <- read.table("../Data/Control_datasets/25S_readcounts.txt", comment.char="#",col.names=c("chromosome","position","5S_DMSO_1","5S_DMSO_2","5S_DMSO_3","5S_DMS_1", "5S_DMS_2", "5S_DMS_3"))
mergedstartswt <- read.table("../Data/Control_datasets/25S_dropoffcounts.txt",comment.char="#",col.names=c("chromosome","position","5S_DMSO_1","5S_DMSO_2","5S_DMSO_3","5S_DMS_1", "5S_DMS_2", "5S_DMS_3"))
mergedcountsmut <- read.table("../Data/Control_datasets/25S_readcounts.txt", comment.char="#",col.names=c("chromosome","position","5S_DMSO_1","5S_DMSO_2","5S_DMSO_3","5S_DMS_1", "5S_DMS_2", "5S_DMS_3"))
mergedstartsmut <- read.table("../Data/Control_datasets/25S_dropoffcounts.txt",comment.char="#",col.names=c("chromosome","position","5S_DMSO_1","5S_DMSO_2","5S_DMSO_3","5S_DMS_1", "5S_DMS_2", "5S_DMS_3"))

head(mergedcountswt)

refsequence <- "../Reference_sequences/Xist.seq"

logdropoffswt <- calculateLDRs(mergedcountswt,mergedstartswt,3,refsequence)
logdropoffsmut <- calculateLDRs(mergedcountsmut,mergedstartsmut3,refsequence)

#prints out the null distribution histogram, the argument breaks is use to define 
#number of bins we want to break the data up into
hist(logdropoffswt$LDR_C, breaks = 30, main = 'mature rRNA dataset - 25S LDR null distribution')
hist(logdropoffswt$LDR_CT, breaks = 30, main = 'mature rRNA dataset - 25S treated vs. control sample LDR distribution')
hist(logdropoffsmut$LDR_C, breaks = 30, main = 'Null distribution of LDRs')
hist(logdropoffsmut$LDR_CT, breaks = 30, main = 'Null distribution of LDRs')

## ------------------------------------------------------------------------
###check if the matrices of p-values can be called after the pipeline has been run twice
head(logdropoffswt$LDR_C)
head(logdropoffswt$LDR_CT)
head(logdropoffsmut$LDR_C)
head(logdropoffsmut$LDR_CT)

Nc <- Nt <- 3

strand = "+"

###
empPvals_1 <- computePvals(logdropoffswt$LDR_C,logdropoffswt$LDR_CT, Nc, Nt, strand, logdropoffswt$nuclPosition,
                           logdropoffswt$nuclSelection$analysedC, logdropoffswt$nuclSelection$analysedCT)

empPvals_2 <- computePvals(logdropoffsmut$LDR_C,logdropoffsmut$LDR_CT, Nc, Nt, strand, logdropoffsmut$nuclPosition,
                           logdropoffsmut$nuclSelection$analysedC, logdropoffsmut$nuclSelection$analysedCT)


stretches <-overlapsRanges(logdropoffswt$stretches,logdropoffsmut$stretches)

## Number of nucleotides in the sequence = number of rows in empPvals_1
nNucl <- length(empPvals_1[1, ])


## ------------------------------------------------------------------------
###computes posterior probabilities of all nucleotides in the stretch specified above.
#This is the step where the null distributions and p-values are calculated as well.
#The most important arguments are the LDRs, the positions used to compute the
#null distribution, as well as the positional information of the selected stretches
#and nucleotide pairs where the LDRs were obtained.

###There are two alternatives, unhash as needed:
###Alternative 1: Calculate posteriors on the data at positions specified by stretches only:  
Pv1 <- matrix(ncol = 1,nrow = length(empPvals_1[,1]))
Pv2 <- matrix(ncol = 1,nrow = length(empPvals_2[,1]))
pvaluesstretch <-list(Pv1, Pv2)

for (i in 1:length(stretches)) {
  if (i>1 & i<length(stretches)) {
    ## Extract start and end of a current stretch from the stretches object
    ## This condition is executed only for IRanges objects in stretches that are not the first or last
    stretchStart <- start(stretches)[i]
    stretchEnd <- end(stretches)[i]
    ## This loop checks stretches to see which positions meets the threshold for analysis, and extracts the corresponding pvalues from empPvals.
    ## The positions that are excluded from analysis (not in stretches) are given empty cells
    Pv1 <-cbind(Pv1, matrix(nrow = length(empPvals_1[,1]), ncol = (stretchStart - end(stretches[i-1])-1)))
    Pv2 <-cbind(Pv2, matrix(nrow = length(empPvals_2[,1]), ncol = (stretchStart - end(stretches[i-1])-1)))
    Pv1 <- cbind(Pv1,empPvals_1[,stretchStart:stretchEnd])
    Pv2 <- cbind(Pv2,empPvals_2[,stretchStart:stretchEnd])
    pvaluesstretch <-list(Pv1, Pv2)
    next()
  } else if (i==1) {
    ## Extract start and end of a current stretch from the stretches object
    ## This condition is only executed for the first IRanges object in stretches
    stretchStart <- start(stretches)[i]
    stretchEnd <- end(stretches)[i]
    
    ## This condition assesses if the first object in stretches is the first nucleotide in the RNA molecule. If not, empty columns are added to the
    ## beginning of the Pv1 and 2 matrix.This ensures that the indexing of the output posteriors match the exact positions on the RNA molecule!
    if (stretchStart != 1) {
      Pv1 <- matrix(nrow = length(empPvals_1[,1]), ncol = (stretchStart-1))
      Pv2 <- matrix(nrow = length(empPvals_2[,1]), ncol = (stretchStart-1))
      Pv1 <- cbind(Pv1, empPvals_1[,stretchStart:stretchEnd])
      Pv2 <- cbind(Pv2, empPvals_2[,stretchStart:stretchEnd])
      
    } else if(stretchStart == 1) {
      Pv1 <- empPvals_1[,stretchStart:stretchEnd]
      Pv2 <- empPvals_2[,stretchStart:stretchEnd]
    }
    pvaluesstretch <- list(Pv1,Pv2)
    next()
    
  } else if (i == length(stretches)) {
    ## Extract start and end of a current stretch from the stretches object
    ## This condition is only executed for the last IRanges objects in stretches
    stretchStart <- start(stretches)[i]
    stretchEnd <- end(stretches)[i]
    
    ##In the end we add empty coulumns to the end of the Pv1 and 2 matrix as needed to match the actual length of the RNA molecule.
    Pv1 <-cbind(Pv1, matrix(nrow = length(empPvals_1[,1]), ncol = (stretchStart - end(stretches[i-1])-1)))
    Pv2 <-cbind(Pv2, matrix(nrow = length(empPvals_2[,1]), ncol = (stretchStart - end(stretches[i-1])-1)))
    Pv1 <- cbind(Pv1,empPvals_1[,stretchStart:stretchEnd], matrix(nrow = length(empPvals_1[,1]), ncol = (length(empPvals_1[1,])-stretchEnd)))
    Pv2 <- cbind(Pv2,empPvals_2[,stretchStart:stretchEnd], matrix(nrow = length(empPvals_2[,1]), ncol = (length(empPvals_2[1,])-stretchEnd)))
    pvaluesstretch <- list(Pv1,Pv2)
  }
  
}

##TEST if pvaluesstretch contains p-values
#pvaluesstretch [[1]][,100:200]

posteriors_diff <- hmmFwbw_differential_two_betas(pvaluesstretch)
colnames(posteriors_diff) <- c("UU","UM","MU","MM")
head(posteriors_diff)

#save this for trying to implement stretches
#for (i in 1:length(stretches)) {
## Extract start and end of a current stretch
#stretchStart <- start(stretches)[i]
#stretchEnd <- end(stretches)[i]
#stretchpvalues <- list(empPvals_1[, stretchStart:stretchEnd],empPvals_2[, stretchStart:stretchEnd])
#if (strand == '+') {
#selection <- hmmFwbw_differential_two_betas(stretchpvalues)
#posteriors_diff[stretchStart:stretchEnd, ] <- t(selection)
#}
#}


### Alternative 2 - Calculating posteriors on all the data without considering stretches:

#posteriors_diff <- hmmFwbw_differential_two_betas(list(empPvals_1,empPvals_2))
#colnames(posteriors_diff) <- c(" ","UU","UM","MU","MM")

#head(posteriors_diff)

## ------------------------------------------------------------------------
#head(posteriors)


## All posterior probabilities need to be shifted by 1 nt because the RT
## is believed to stop 1 nucleotide before the modified nucleodide.
## So below a new matrix is made containing the values shifted by one position.
## ------------------------------------------------------------------------
shifted_posteriors <- matrix(, nrow=dim(posteriors_diff)[1], ncol=4)
shifted_posteriors[1:(length(shifted_posteriors[,1]) - 1), ] <- posteriors_diff[2:(length(shifted_posteriors[,1])), ]
colnames(shifted_posteriors) <- c("UU","UM","MU","MM")

head(shifted_posteriors)
head(posteriors_diff)
posteriors_diff

differentiallymod <- shifted_posteriors[,2] + shifted_posteriors[,3]

## ------------------------------------------------------------------------
plot(differentiallymod, xlab = 'Nucleotide position',
     ylab = 'Probability of differential modification',
     main = 'Mature yeast rRNA dataset: Probabilites of differential modification - 5.8S',
     ylim = c(0,1))

## ----eval=FALSE----------------------------------------------------------
## ## Call the function with the additonal tolerance parameter
## posteriors <- computeProbs(LDR_C, LDR_CT, Nc, Nt, '+', nuclPosition,
##                            nuclSelection$analysedC, nuclSelection$analysedCT,
##                            stretches, 0.001)

## ------------------------------------------------------------------------
shifted_posteriors <- replace(shifted_posteriors,is.na(shifted_posteriors),-999)
write.table(shifted_posteriors,sep="\t",quote=FALSE,file="25S_diffBUM_HMM_samedataset.txt",col.names = c("UU","UM","MU","MM"), row.names = TRUE)
