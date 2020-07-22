#### ------- PACKAGES INSTALLATION AND IMPORT OF HELPER FUNCTIONS ------ ######

# This script assumes: R version 4.0.0 (2020-04-24); RStudio Version 1.2.5001

install.packages("rstudioapi")
library(rstudioapi)

install.packages("BiocManager")
install.packages("formattable")

BiocManager::install(c("Biostrings", "SummarizedExperiment"), version = "3.11")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
###

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

#### ------- LOADING DATA FOR ALL 4 rRNAs ------ ######

noreplicates <- 3
setwd("Reference_sequences")
refsequence25S <- "25S_refseq.txt"
refsequence18S <- "18S_refseq.txt"
refsequence5.8S <- "5.8S_refseq.txt"
refsequence5S <- "5S_refseq.txt"
setwd('../Data/Control_rRNA_dataset')


## input the coverage and drop off counts for each of the rRNAs from the working directory:
# wt:wild type; mut: mutant, for this control dataset wt and mut are the same, to test for false positives

# 25S
mergedcountswt25S <- read.table("25S_readcounts.txt", comment.char="#",col.names=c("chromosome","position","25S_DMSO_1","25S_DMSO_2","25S_DMSO_3","25S_DMS_1", "25S_DMS_2", "25S_DMS_3"))
mergedstartswt25S <- read.table("25S_dropoffcounts.txt",comment.char="#",col.names=c("chromosome","position","25S_DMSO_1","25S_DMSO_2","25S_DMSO_3","25S_DMS_1", "25S_DMS_2", "25S_DMS_3"))
mergedcountsmut25S <- read.table("25S_readcounts.txt", comment.char="#",col.names=c("chromosome","position","25S_DMSO_1","25S_DMSO_2","25S_DMSO_3","25S_DMS_1", "25S_DMS_2", "25S_DMS_3"))
mergedstartsmut25S <- read.table("25S_dropoffcounts.txt",comment.char="#",col.names=c("chromosome","position","25S_DMSO_1","25S_DMSO_2","25S_DMSO_3","25S_DMS_1", "25S_DMS_2", "25S_DMS_3"))

# 18S
mergedcountswt18S <- read.table("18S_readcounts.txt", comment.char="#",col.names=c("chromosome","position","18S_DMSO_1","18S_DMSO_2","18S_DMSO_3","18S_DMS_1", "18S_DMS_2", "18S_DMS_3"))
mergedstartswt18S <- read.table("18S_dropoffcounts.txt",comment.char="#",col.names=c("chromosome","position","18S_DMSO_1","18S_DMSO_2","18S_DMSO_3","18S_DMS_1", "18S_DMS_2", "18S_DMS_3"))
mergedcountsmut18S <- read.table("18S_readcounts.txt", comment.char="#",col.names=c("chromosome","position","18S_DMSO_1","18S_DMSO_2","18S_DMSO_3","18S_DMS_1", "18S_DMS_2", "18S_DMS_3"))
mergedstartsmut18S <- read.table("18S_dropoffcounts.txt",comment.char="#",col.names=c("chromosome","position","18S_DMSO_1","18S_DMSO_2","18S_DMSO_3","18S_DMS_1", "18S_DMS_2", "18S_DMS_3"))

# 5.8S
mergedcountswt5.8S <- read.table("5.8S_readcounts.txt", comment.char="#",col.names=c("chromosome","position","5.8S_DMSO_1","5.8S_DMSO_2","5.8S_DMSO_3","5.8S_DMS_1", "5.8S_DMS_2", "5.8S_DMS_3"))
mergedstartswt5.8S <- read.table("5.8S_dropoffcounts.txt",comment.char="#",col.names=c("chromosome","position","5.8S_DMSO_1","5.8S_DMSO_2","5.8S_DMSO_3","5.8S_DMS_1", "5.8S_DMS_2", "5.8S_DMS_3"))
mergedcountsmut5.8S <- read.table("5.8S_readcounts.txt", comment.char="#",col.names=c("chromosome","position","5.8S_DMSO_1","5.8S_DMSO_2","5.8S_DMSO_3","5.8S_DMS_1", "5.8S_DMS_2", "5.8S_DMS_3"))
mergedstartsmut5.8S <- read.table("5.8S_dropoffcounts.txt",comment.char="#",col.names=c("chromosome","position","5.8S_DMSO_1","5.8S_DMSO_2","5.8S_DMSO_3","5.8S_DMS_1", "5.8S_DMS_2", "5.8S_DMS_3"))

# 5S
mergedcountswt5S <- read.table("5S_readcounts.txt", comment.char="#",col.names=c("chromosome","position","5S_DMSO_1","5S_DMSO_2","5S_DMSO_3","5S_DMS_1", "5S_DMS_2", "5S_DMS_3"))
mergedstartswt5S <- read.table("5S_dropoffcounts.txt",comment.char="#",col.names=c("chromosome","position","5S_DMSO_1","5S_DMSO_2","5S_DMSO_3","5S_DMS_1", "5S_DMS_2", "5S_DMS_3"))
mergedcountsmut5S <- read.table("5S_readcounts.txt", comment.char="#",col.names=c("chromosome","position","5S_DMSO_1","5S_DMSO_2","5S_DMSO_3","5S_DMS_1", "5S_DMS_2", "5S_DMS_3"))
mergedstartsmut5S <- read.table("5S_dropoffcounts.txt",comment.char="#",col.names=c("chromosome","position","5S_DMSO_1","5S_DMSO_2","5S_DMSO_3","5S_DMS_1", "5S_DMS_2", "5S_DMS_3"))

setwd("./../../")
#### ------- DATA PRE-PROCESSING (CALCULATING LOG RATIOS OF MUTATION RATES AND EMPIRICAL P-VALUES ) ------- ######

# Calculating log of drop off rate ratios

# 25S
logdropoffswt25S <- calculateLDRs(mergedcountswt25S,mergedstartswt25S,noreplicates,refsequence25S)
logdropoffsmut25S <- calculateLDRs(mergedcountsmut25S,mergedstartsmut25S,noreplicates,refsequence25S)

# 18S
logdropoffswt18S <- calculateLDRs(mergedcountswt18S,mergedstartswt18S,noreplicates,refsequence18S)
logdropoffsmut18S <- calculateLDRs(mergedcountsmut18S,mergedstartsmut18S,noreplicates,refsequence18S)

# 5.8S
logdropoffswt5.8S <- calculateLDRs(mergedcountswt5.8S,mergedstartswt5.8S,noreplicates,refsequence5.8S)
logdropoffsmut5.8S <- calculateLDRs(mergedcountsmut5.8S,mergedstartsmut5.8S,noreplicates,refsequence5.8S)

# 5S
logdropoffswt5S <- calculateLDRs(mergedcountswt5S,mergedstartswt5S,noreplicates,refsequence5S)
logdropoffsmut5S <- calculateLDRs(mergedcountsmut5S,mergedstartsmut5S,noreplicates,refsequence5S)


## ------- QUALITY CONTROL: INSPECTION OF LOG DROP OFF RATE RATIOS DISTRIBUTION (OPTIONAL)------ ##

# TO RUN,UNCOMMENT LINES 98-136 WITH: CTRL + SHIFT + C On Windows or command + SHIFT + C on Mac OS 
# WT AND MUT LDRs WILL BE THE SAME SINCE THEY ARE IDENTICAL SETS OF SAMPLES, HERE ONLY THE CODE FOR ASSESSING WT IS SHOWN.

# setwd("Analysis/LMR_and_LDR_plots/mature rRNA")
# 
# # 25S
# pdf('LDR null distribution-mature rRNA-25S.pdf',width=6,height=4,paper='special')
# hist(logdropoffswt25S$LDR_C, breaks = 30, main = '25S - LDR null distribution')
# dev.off()
# 
# pdf('CT LDR distribution-mature rRNA-25S.pdf',width=6,height=4,paper='special')
# hist(logdropoffswt25S$LDR_CT, breaks = 30, main = '25S - treated vs.control sample LDR distribution')
# dev.off()
# 
# # 18S
# pdf('LDR null distribution-mature rRNA-18S.pdf',width=6,height=4,paper='special')
# hist(logdropoffswt18S$LDR_C, breaks = 30, main = '18S - LDR null distribution')
# dev.off()
# 
# pdf('CT LDR distribution-mature rRNA-18S.pdf',width=6,height=4,paper='special')
# hist(logdropoffswt18S$LDR_CT, breaks = 30, main = '18S - treated vs.control sample LDR distribution')
# dev.off()
# 
# # 5.8S
# pdf('LDR null distribution-mature rRNA-5.8S.pdf',width=6,height=4,paper='special')
# hist(logdropoffswt5.8S$LDR_C, breaks = 30, main = '5.8S - LDR null distribution')
# dev.off()
# 
# pdf('CT LDR distribution-mature rRNA-5.8S.pdf',width=6,height=4,paper='special')
# hist(logdropoffswt5.8S$LDR_CT, breaks = 30, main = '5.8S - treated vs.control sample LDR distribution')
# dev.off()
# 
# # 5S
# pdf('LDR null distribution-mature rRNA-5S.pdf',width=6,height=4,paper='special')
# hist(logdropoffswt5S$LDR_C, breaks = 30, main = '5S - LDR null distribution')
# dev.off()
# 
# pdf('CT LDR distribution-mature rRNA-5S.pdf',width=6,height=4,paper='special')
# hist(logdropoffswt5S$LDR_CT, breaks = 30, main = '5S - treated vs.control sample LDR distribution')
# dev.off()

# setwd('./../../..')


## ------------------------------------------------------------------------##
head(logdropoffswt25S$LDR_C)
head(logdropoffswt18S$LDR_CT)
head(logdropoffsmut5.8S$LDR_C)
head(logdropoffsmut5S$LDR_CT)

Nc <- Nt <- noreplicates

strand = "+"

##Calculation of empirical p values for the two sets of samples, which normally represent differential probing conditions

# 25S
empPvals_25S_A <- computePvals(logdropoffswt25S$LDR_C,logdropoffswt25S$LDR_CT, Nc, Nt, strand, logdropoffswt25S$nuclPosition,
                           logdropoffswt25S$nuclSelection$analysedC, logdropoffswt25S$nuclSelection$analysedCT)

empPvals_25S_B <- computePvals(logdropoffsmut25S$LDR_C,logdropoffsmut25S$LDR_CT, Nc, Nt, strand, logdropoffsmut25S$nuclPosition,
                           logdropoffsmut25S$nuclSelection$analysedC, logdropoffsmut25S$nuclSelection$analysedCT)

#stretches contain the selection of nulceotide positions that have valid log drop off rate ratios
#overlapping the stretches of nucleotides selected in each group of samples to get a representative set
stretches_25S <-overlapsRanges(logdropoffswt25S$stretches,logdropoffsmut25S$stretches)


# 18S
empPvals_18S_A <- computePvals(logdropoffswt18S$LDR_C,logdropoffswt18S$LDR_CT, Nc, Nt, strand, logdropoffswt18S$nuclPosition,
                               logdropoffswt18S$nuclSelection$analysedC, logdropoffswt18S$nuclSelection$analysedCT)

empPvals_18S_B <- computePvals(logdropoffsmut18S$LDR_C,logdropoffsmut18S$LDR_CT, Nc, Nt, strand, logdropoffsmut18S$nuclPosition,
                               logdropoffsmut18S$nuclSelection$analysedC, logdropoffsmut18S$nuclSelection$analysedCT)

stretches_18S <-overlapsRanges(logdropoffswt18S$stretches,logdropoffsmut18S$stretches)



# 5.8S
empPvals_5.8S_A <- computePvals(logdropoffswt5.8S$LDR_C,logdropoffswt5.8S$LDR_CT, Nc, Nt, strand, logdropoffswt5.8S$nuclPosition,
                               logdropoffswt5.8S$nuclSelection$analysedC, logdropoffswt5.8S$nuclSelection$analysedCT)

empPvals_5.8S_B <- computePvals(logdropoffsmut5.8S$LDR_C,logdropoffsmut5.8S$LDR_CT, Nc, Nt, strand, logdropoffsmut5.8S$nuclPosition,
                               logdropoffsmut5.8S$nuclSelection$analysedC, logdropoffsmut5.8S$nuclSelection$analysedCT)

stretches_5.8S <-overlapsRanges(logdropoffswt5.8S$stretches,logdropoffsmut5.8S$stretches)



# 5S
empPvals_5S_A <- computePvals(logdropoffswt5S$LDR_C,logdropoffswt5S$LDR_CT, Nc, Nt, strand, logdropoffswt5S$nuclPosition,
                               logdropoffswt5S$nuclSelection$analysedC, logdropoffswt5S$nuclSelection$analysedCT)

empPvals_5S_B <- computePvals(logdropoffsmut5S$LDR_C,logdropoffsmut5S$LDR_CT, Nc, Nt, strand, logdropoffsmut5S$nuclPosition,
                               logdropoffsmut5S$nuclSelection$analysedC, logdropoffsmut5S$nuclSelection$analysedCT)

stretches_5S <-overlapsRanges(logdropoffswt5S$stretches,logdropoffsmut5S$stretches)



####### ------- HMM ANALYSIS (CALCULATION OF POSTERIOR PROBABILITIES FOR THE FOUR HIDDEN STATES) ------- ########
###computes posterior probabilities of all nucleotides in the stretches specified above.
#The most important arguments are the LDRs, the positions used to compute the
#null distribution, as well as the positions of the selected stretches
#and nucleotide pairs where the LDRs were obtained.


##Before computing posterior probabilities, we set up matrices to hold the required p-values
selectPvalues <- function(empPvals_1, empPvals_2, stretches){
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
  return(pvaluesstretch)
}

pvaluesstretch_25S <- selectPvalues(empPvals_25S_A, empPvals_25S_B, stretches_25S)
pvaluesstretch_18S <- selectPvalues(empPvals_18S_A, empPvals_18S_B, stretches_18S)
pvaluesstretch_5.8S <- selectPvalues(empPvals_5.8S_A, empPvals_5.8S_B, stretches_5.8S)
pvaluesstretch_5S <- selectPvalues(empPvals_5S_A, empPvals_5S_B, stretches_5S)


## ------- QUALITY CONTROL: INSPECTION OF P-VALUES (OPTIONAL) ------- ##
#TO RUN,UNCOMMENT LINES 270-294 WITH: CTRL + SHIFT + C On Windows or command + SHIFT + C on Mac OS 

# setwd("Analysis/pvalues_plots/mature rRNA")
# 
# write.table(t(pvaluesstretch_25S[[1]]), file="pvalues_25S.txt", row.names=TRUE, col.names=TRUE)
# write.table(t(pvaluesstretch_18S[[1]]), file="pvalues_18S.txt", row.names=TRUE, col.names=TRUE)
# write.table(t(pvaluesstretch_5.8S[[1]]), file="pvalues_5.8S.txt", row.names=TRUE, col.names=TRUE)
# write.table(t(pvaluesstretch_5S[[1]]), file="pvalues_5S.txt", row.names=TRUE, col.names=TRUE)
# 
# 
# pdf('pvalues-maturerRNA-25S.pdf',width=6,height=4,paper='special')
# hist(pvaluesstretch_25S[[1]], breaks = 30, main = 'mature rRNA dataset ??? 25S p???value distribution')
# dev.off()
# 
# pdf('pvalues-maturerRNA-18S.pdf',width=6,height=4,paper='special')
# hist(pvaluesstretch_18S[[1]], breaks = 30, main = 'mature rRNA dataset ??? 18S p???value distribution')
# dev.off()
# 
# pdf('pvalues-maturerRNA-5.8S.pdf',width=6,height=4,paper='special')
# hist(pvaluesstretch_5.8S[[1]], breaks = 30, main = 'mature rRNA dataset ??? 5.8S p???value distribution')
# dev.off()
# 
# pdf('pvalues-maturerRNA-5S.pdf',width=6,height=4,paper='special')
# hist(pvaluesstretch_5S[[1]], breaks = 30, main = 'mature rRNA dataset ??? 5S p???value distribution')
# dev.off()
# 
# setwd('./../../..')

## ------------------------------------------------------------------ ##

#The selected empirical p values are input into the HMM, to give posterior probabilities
posteriors_diff_25S <- hmmFwbw_differential_two_betas(pvaluesstretch_25S)
posteriors_diff_18S <- hmmFwbw_differential_two_betas(pvaluesstretch_18S)
posteriors_diff_5.8S <- hmmFwbw_differential_two_betas(pvaluesstretch_5.8S)
posteriors_diff_5S <- hmmFwbw_differential_two_betas(pvaluesstretch_5S)

colnames(posteriors_diff_25S) <- c("UU","UM","MU","MM")
colnames(posteriors_diff_18S) <- c("UU","UM","MU","MM")
colnames(posteriors_diff_5.8S) <- c("UU","UM","MU","MM")
colnames(posteriors_diff_5S) <- c("UU","UM","MU","MM")
head(posteriors_diff_25S)



## ------- SHIFTING POSTERIORS (OPTIONAL) ------- ##
## Applicable only to RT-stop probing chemistries.
## All posterior probabilities need to be shifted by 1 nt because the RT
## is believed to stop 1 nucleotide before the modified nucleodide.
## So below a new matrix is made containing the values shifted by one position.
## SKIP THIS PART IF THE INPUT DATA IS GENERATED WITH RT-mutate probing techniques(eg. SHAPE-MaP)!!!

shifted_posteriors_25S <- matrix(, nrow=dim(posteriors_diff_25S)[1], ncol=4)
shifted_posteriors_25S[1:(length(shifted_posteriors_25S[,1]) - 1), ] <- posteriors_diff_25S[2:(length(shifted_posteriors_25S[,1])), ]
colnames(shifted_posteriors_25S) <- c("UU","UM","MU","MM")

shifted_posteriors_18S <- matrix(, nrow=dim(posteriors_diff_18S)[1], ncol=4)
shifted_posteriors_18S[1:(length(shifted_posteriors_18S[,1]) - 1), ] <- posteriors_diff_18S[2:(length(shifted_posteriors_18S[,1])), ]
colnames(shifted_posteriors_18S) <- c("UU","UM","MU","MM")

shifted_posteriors_5.8S <- matrix(, nrow=dim(posteriors_diff_5.8S)[1], ncol=4)
shifted_posteriors_5.8S[1:(length(shifted_posteriors_5.8S[,1]) - 1), ] <- posteriors_diff_5.8S[2:(length(shifted_posteriors_5.8S[,1])), ]
colnames(shifted_posteriors_5.8S) <- c("UU","UM","MU","MM")

shifted_posteriors_5S <- matrix(, nrow=dim(posteriors_diff_5S)[1], ncol=4)
shifted_posteriors_5S[1:(length(shifted_posteriors_5S[,1]) - 1), ] <- posteriors_diff_5S[2:(length(shifted_posteriors_5S[,1])), ]
colnames(shifted_posteriors_5S) <- c("UU","UM","MU","MM")

head(posteriors_diff_25S)
head(shifted_posteriors_25S)


##### ------------------ EXPORTING POSTERIORS AND PLOTS -------------------------------#####

#Plotting and outputting the shifted posterior probabilities of differential modification
differentiallymod_25S <- shifted_posteriors_25S[,2] + shifted_posteriors_25S[,3]
differentiallymod_18S <- shifted_posteriors_18S[,2] + shifted_posteriors_18S[,3]
differentiallymod_5.8S <- shifted_posteriors_5.8S[,2] + shifted_posteriors_5.8S[,3]
differentiallymod_5S <- shifted_posteriors_5S[,2] + shifted_posteriors_5S[,3]

setwd("Analysis/diffBUM-HMM")

# 25S
pdf('rRNA_25S_sum_of_diff_states_diff_BUM_HMM.pdf', width = 10)
plot(differentiallymod_25S, xlab = 'Nucleotide position',
     ylab = 'Probability of differential modification',
     main = 'Mature yeast rRNA dataset: Probabilites of differential modification - 25S',
     ylim = c(0,1))
dev.off()

# 18S
pdf('rRNA_18S_sum_of_diff_states_diff_BUM_HMM.pdf', width = 10)
plot(differentiallymod_18S, xlab = 'Nucleotide position',
     ylab = 'Probability of differential modification',
     main = 'Mature yeast rRNA dataset: Probabilites of differential modification - 18S',
     ylim = c(0,1))
dev.off()

# 5.8S
pdf('rRNA_5.8S_sum_of_diff_states_diff_BUM_HMM.pdf', width = 10)
plot(differentiallymod_5.8S, xlab = 'Nucleotide position',
     ylab = 'Probability of differential modification',
     main = 'Mature yeast rRNA dataset: Probabilites of differential modification - 5.8S',
     ylim = c(0,1))
dev.off()

# 5S
pdf('rRNA_5S_sum_of_diff_states_diff_BUM_HMM.pdf', width = 10)
plot(differentiallymod_5S, xlab = 'Nucleotide position',
     ylab = 'Probability of differential modification',
     main = 'Mature yeast rRNA dataset: Probabilites of differential modification - 5S',
     ylim = c(0,1))
dev.off()


#Outputting all posterior probabilities
shifted_posteriors_25S <- replace(shifted_posteriors_25S,is.na(shifted_posteriors_25S),-999)
shifted_posteriors_18S <- replace(shifted_posteriors_18S,is.na(shifted_posteriors_18S),-999)
shifted_posteriors_5.8S <- replace(shifted_posteriors_5.8S,is.na(shifted_posteriors_5.8S),-999)
shifted_posteriors_5S <- replace(shifted_posteriors_5S,is.na(shifted_posteriors_5S),-999)

write.table(shifted_posteriors_25S,sep="\t",quote=FALSE,file="mature_rRNA_25S_control_identical_conditions_diff_BUM_HMM.txt",col.names = c("UU","UM","MU","MM"), row.names = TRUE)
write.table(shifted_posteriors_18S,sep="\t",quote=FALSE,file="mature_rRNA_18S_control_identical_conditions_diff_BUM_HMM.txt",col.names = c("UU","UM","MU","MM"), row.names = TRUE)
write.table(shifted_posteriors_5.8S,sep="\t",quote=FALSE,file="mature_rRNA_5.8S_control_identical_conditions_diff_BUM_HMM.txt",col.names = c("UU","UM","MU","MM"), row.names = TRUE)
write.table(shifted_posteriors_5S,sep="\t",quote=FALSE,file="mature_rRNA_5S_control_identical_conditions_diff_BUM_HMM.txt",col.names = c("UU","UM","MU","MM"), row.names = TRUE)
