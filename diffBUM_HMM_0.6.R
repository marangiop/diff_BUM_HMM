### CHANGE THE DIRECTORY TO RUN FROM THE CORRECT FOLDER!!!

working_directory <-"/Users/maran/Desktop/diff_BUM_HMM_Project/Github/diff_BUM_HMM/"
setwd(working_directory)  


#biocLite("BUMHMM")
#biocLite("Biostrings")
#install.packages("Rcpp", dependencies = TRUE)

source("https://bioconductor.org/biocLite.R")
source("Functions/computePvals.R")
source("Functions/calculateLDRs.R")
source("Functions/hmmFwbw_differential_two_betas.R")
source("Functions/betaParamsEM.R")
source("Functions/betaParamsMStep.R")
source("Functions/findPatternPos.R")
source("Functions/nuclPerm.R")
source("Functions/scaleDOR.R")
source("Functions/selectNuclPos.R")
source("Functions/computeStretches.R")
source("Functions/stabiliseVariance.R")

suppressPackageStartupMessages({ 
  library(Biostrings)
  library(SummarizedExperiment) })

noreplicates <- 2

ref_seq_directory <- paste(working_directory, "Reference sequences/" ,sep="")
setwd(ref_seq_directory)  


refsequence <- "Xist.seq"

setwd(working_directory)  	

outputfilename <-paste0('Xist_in vivo_vs_ex vivo','_diff_BUM_HMM_analysed_0.5noise','.txt')

mergedcountswt <- read.table("Data/Xist_invivo_reads.txt", comment.char="#",col.names=c("chromosome","position","DMSO_1","DMSO_2","1M7_1","1M7_2"))
mergedstartswt <- read.table("Data/Xist_invivo_substitutions.txt",comment.char="#",col.names=c("chromosome","position","DMSO_1","DMSO_2","1M7_1","1M7_2"))
mergedcountsmut <- read.table("Data/Xist_exvivo_reads.txt", comment.char="#",col.names=c("chromosome","position","DMSO_1","DMSO_2","1M7_1","1M7_2"))
mergedstartsmut <- read.table("Data/Xist_exvivo_substitutions.txt",comment.char="#",col.names=c("chromosome","position","DMSO_1","DMSO_2","1M7_1","1M7_2"))

logdropoffswt <- calculateLDRs(mergedcountswt,mergedstartswt, noreplicates, refsequence, working_directory)
logdropoffsmut <- calculateLDRs(mergedcountsmut,mergedstartsmut, noreplicates, refsequence, working_directory)
#prints out the null distribution histogram, the argument breaks is use to define 
#number of bins we want to break the data up into
hist(logdropoffswt$LDR_C, breaks = 30, main = 'Null distribution of LDRs')

## ------------------------------------------------------------------------
###check if the matrices of p-values can be called after the pipeline has been run twice
head(logdropoffswt$LDR_C)
head(logdropoffswt$LDR_CT)
head(logdropoffsmut$LDR_C)
head(logdropoffsmut$LDR_CT)

Nc <- Nt <- noreplicates

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
pvaluesstretch<-list(Pv1, Pv2)
for (i in 1:length(stretches)) {
  if (i>1 & i<=length(stretches)) {
    ## Extract start and end of a current stretch
    stretchStart <- start(stretches)[i]
    stretchEnd <- end(stretches)[i]
    Pv1 <-cbind(Pv1, matrix(nrow = length(empPvals_1[,1]), ncol = (stretchStart - end(stretches[i-1])-1)))
    Pv2 <-cbind(Pv2, matrix(nrow = length(empPvals_2[,1]), ncol = (stretchStart - end(stretches[i-1])-1)))
    Pv1 <- cbind(Pv1,empPvals_1[,stretchStart:stretchEnd])
    Pv2 <- cbind(Pv2,empPvals_2[,stretchStart:stretchEnd])
    pvaluesstretch <-list(Pv1, Pv2)
    next()
  } else {
    ## Extract start and end of a current stretch
    stretchStart <- start(stretches)[i]
    stretchEnd <- end(stretches)[i]
    Pv1 <- cbind(Pv1,empPvals_1[,stretchStart:stretchEnd])
    Pv2 <- cbind(Pv2,empPvals_2[,stretchStart:stretchEnd])
    pvaluesstretch <- list(Pv1,Pv2)
    next()
  }
  return(pvaluesstretch)
}

##TEST if pvaluesstretch contains p-values
#pvaluesstretch [[1]][,100:200]
noise=0.5
posteriors_diff <- hmmFwbw_differential_two_betas(pvaluesstretch, noise)
#colnames(posteriors_diff) <- c("  ","UU","UM","MU","MM") - The way Toby had this line
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
#colnames(shifted_posteriors) <- c(" ","UU","UM","MU","MM")
colnames(shifted_posteriors) <- c("UU","UM","MU","MM")


head(shifted_posteriors)
head(posteriors_diff)
posteriors_diff

differentiallymod <- shifted_posteriors[,2] + shifted_posteriors[,3]

## ------------------------------------------------------------------------
plot(differentiallymod, xlab = 'Nucleotide position',
     ylab = 'Probability of modification',
     main = 'diffBUMHMM output: Probabilites of differential modification between Erb1 and delta 5',
     ylim = c(0,1))

## ----eval=FALSE----------------------------------------------------------
## ## Call the function with the additonal tolerance parameter
## posteriors <- computeProbs(LDR_C, LDR_CT, Nc, Nt, '+', nuclPosition,
##                            nuclSelection$analysedC, nuclSelection$analysedCT,
##                            stretches, 0.001)

## ------------------------------------------------------------------------
shifted_posteriors <- replace(shifted_posteriors,is.na(shifted_posteriors),-999)
write.table(shifted_posteriors,sep="\t",quote=FALSE,file=outputfilename,col.names = c("UU","UM","MU","MM"), row.names = TRUE)
