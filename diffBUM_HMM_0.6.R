### CHANGE THE DIRECTORY TO RUN FROM THE CORRECT FOLDER!!!

working_directory <-"/Users/maran/Desktop/diff_BUM_HMM_Project/Github/diff_BUM_HMM/"
setwd(working_directory)  


#biocLite("BUMHMM")
#biocLite("Biostrings")
#install.packages("Rcpp", dependencies = TRUE)

source("https://bioconductor.org/biocLite.R")
source("Functions/computePvals.R")
source("Functions/calculateLDRs.R")
#source("Functions/hmmFwbw_differential_two_betas.R")
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

#noise=-0.5
refsequence <- "Xist.seq"

setwd(working_directory)  	

outputfilename <-paste0('Xist_in vivo_vs_ex vivo','_diff_BUM_HMM_analysed_gaussian_noise','.txt')

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
pValues<-list(Pv1, Pv2)
for (i in 1:length(stretches)) {
  if (i>1 & i<=length(stretches)) {
    ## Extract start and end of a current stretch
    stretchStart <- start(stretches)[i]
    stretchEnd <- end(stretches)[i]
    Pv1 <-cbind(Pv1, matrix(nrow = length(empPvals_1[,1]), ncol = (stretchStart - end(stretches[i-1])-1)))
    Pv2 <-cbind(Pv2, matrix(nrow = length(empPvals_2[,1]), ncol = (stretchStart - end(stretches[i-1])-1)))
    Pv1 <- cbind(Pv1,empPvals_1[,stretchStart:stretchEnd])
    Pv2 <- cbind(Pv2,empPvals_2[,stretchStart:stretchEnd])
    pValues <-list(Pv1, Pv2)
    next()
  } else {
    ## Extract start and end of a current stretch
    stretchStart <- start(stretches)[i]
    stretchEnd <- end(stretches)[i]
    Pv1 <- cbind(Pv1,empPvals_1[,stretchStart:stretchEnd])
    Pv2 <- cbind(Pv2,empPvals_2[,stretchStart:stretchEnd])
    pValues <- list(Pv1,Pv2)
    next()
  }
  return(pValues)
}

##TEST if pvaluesstretch contains p-values
#pvaluesstretch [[1]][,100:200]

#posteriors_diff <- hmmFwbw_differential_two_betas(pvaluesstretch)


#hmmFwbw_differential_two_betas <- function(pValues){
    
    ## Settings for HMM
    ## Set the transition matrix: try amending this for performance
a = 0.9025 + rnorm(1, mean=0, sd =0.01)
b = 0.0475 + rnorm(1, mean=0, sd =0.01)
c = 0.0475 + rnorm(1, mean=0, sd =0.01) 
d = 1 - (a + b + c)

e = 0.19 + rnorm(1, mean=0, sd =0.01)
f = 0.76 + rnorm(1, mean=0, sd =0.01)
g = 0.01 + rnorm(1, mean=0, sd =0.01)
h = 1 - (e + f + g)

i = 0.19 + rnorm(1, mean=0, sd =0.01)
j = 0.01 + rnorm(1, mean=0, sd =0.01)
k = 0.76 + rnorm(1, mean=0, sd =0.01)
l = 1 - (i + j + k)

m = 0.04 + rnorm(1, mean=0, sd =0.01)
n = 0.16 + rnorm(1, mean=0, sd =0.01)
o = 0.16 + rnorm(1, mean=0, sd =0.01)
p = 1 - (m + n + o)

trans <- matrix(c(a, e, i, m,
                  b, f, j, n,
                  c, g, k, o,
                  d, h, l, p), nrow = 4, ncol = 4, byrow = TRUE)

## Set the values for Beta shape parameters in the emission mixture
## model
alpha_P1 = 1
beta_P1 = 10

alpha_P2 = 1
beta_P2 = 10

## Set initial probability to 0.25 in each state
initialProb = c(0.25, 0.25, 0.25, 0.25)
## Number of experiments (treatment-control comparisons)
nexp <- length(pValues[[1]][, 1])

# calculates number of nucleotides
nBins <- length(pValues[[1]][1, ]) 

# calculates number of states
nStates <- length(trans[, 1]) 

## Log-likelihood of observations given state
obsLike <- matrix(1, ncol = nBins, nrow = nStates) ##generates a matrix where rows are states and columns are nucleotides

fwdMessage <- matrix(0, ncol = nBins, nrow = nStates)
bwdMessage <- matrix(0, ncol = nBins, nrow = nStates)

## Non independence assumption
## Calculation of likelihoods (sum over replicates of likelihoods of each experiment)
#iterates over a sequence of values from 1 to the number of experimental comparisons
for (index in 1:nexp ) { 
    #iterates over a sequence of values from 1 to the number of nucleotides
    for (index2 in 1:nBins) { 
        # Likelihod for UU
        if (is.na(pValues[[1]][index, index2]) || is.na(pValues[[2]][index, index2])) {
            obsLike[, index2] <- obsLike[, index2]
        } else { 
            obsLike[1, index2] <- obsLike[1, index2] * stats::dbeta(pValues[[1]][index, index2], shape1 = 1, shape2 = 1) * 
                stats::dbeta(pValues[[2]][index, index2], shape1 = 1, shape2 = 1)
            # Likelihod for MU
            obsLike[2, index2] <- obsLike[2, index2] * stats::dbeta(pValues[[1]][index, index2], shape1 = 1, shape2 = 1) *
                stats::dbeta(pValues[[2]][index, index2], shape1 = alpha_P2, shape2 = beta_P2)
            # Likelihod for UM
            obsLike[3, index2] <- obsLike[3, index2] * stats::dbeta(pValues[[1]][index, index2], shape1 = alpha_P1, shape2 = beta_P1) *
                stats::dbeta(pValues[[2]][index, index2], shape1 = 1, shape2 = 1)
            # Likelihod for MM
            obsLike[4, index2] <- obsLike[4, index2] * stats::dbeta(pValues[[1]][index, index2], shape1 = alpha_P1, shape2 = beta_P1) *
                stats::dbeta(pValues[[2]][index, index2], shape1 = alpha_P2, shape2 = beta_P2)
        }
    }
}

## Calculation of the forward messages
fwdMessage[, 1] <- initialProb * obsLike[, 1] 
fwdMessage[, 1] <- fwdMessage[, 1] / sum(fwdMessage[, 1]) 

for (index in 2:nBins) {
    fwdMessage[, index] <- (trans %*% fwdMessage[, index - 1]) *
        obsLike[, index]     
    fwdMessage[, index] <- fwdMessage[, index] / sum(fwdMessage[, index]) 
}

## Calculation of the backward message
bwdMessage[, nBins] <- 1

for (index in (nBins - 1):1) {
    bwdMessage[, index] <- (trans %*% bwdMessage[, index + 1]) * obsLike[, index]
    bwdMessage[, index] <- bwdMessage[, index] / sum(bwdMessage[, index])
}

## Calculation of posteriors
posteriors_diff <- fwdMessage * bwdMessage
posteriors_diff <- posteriors_diff/ (matrix(1, nrow = length(posteriors_diff[, 1]))
                          %*% colSums(posteriors_diff))
posteriors_diff <- t(posteriors_diff)







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
