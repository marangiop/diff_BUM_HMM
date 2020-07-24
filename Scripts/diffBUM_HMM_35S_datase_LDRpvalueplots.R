#### ------- PACKAGES INSTALLATION AND IMPORT OF HELPER FUNCTIONS ------ ######

# This script assumes: R version 4.0.0 (2020-04-24); RStudio Version 1.2.5001

install.packages("rstudioapi")
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

install.packages("BiocManager")

BiocManager::install(c("Biostrings", "SummarizedExperiment"), version = "3.11")

setwd('..')

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


#### ------- LOADING DATA ------ ######

noreplicates <- 2


setwd("Reference_sequences")
refsequence <- "35S_pre-rRNA_refseq.seq"
setwd('../Data/35S_data/')
outputfilename <-paste0('35S','_diffBUM_HMM_WT_vs_Erb1_diff_BUM_HMM_analysed','.txt')

# input the coverage and drop off counts from the working directory:
# wt:wild type; mut: mutant

mergedcountswt <- read.table("35S_control_delta5_merged_reads.sgr", comment.char="#",col.names=c("chromosome","position","35S_control_1_trimmed_plus_strand_startpositions","35S_control_2_trimmed_plus_strand_startpositions","35S_1M7_1_trimmed_plus_strand_startpositions","35S_1M7_2_trimmed_plus_strand_startpositions"))
mergedstartswt <- read.table("35S_control_delta5_merged_dropoffcounts.sgr",comment.char="#",col.names=c("chromosome","position","35S_control_1_trimmed_plus_strand_startpositions","35S_control_2_trimmed_plus_strand_startpositions","35S_1M7_1_trimmed_plus_strand_startpositions","35S_1M7_2_trimmed_plus_strand_startpositions"))
mergedcountsmut <- read.table("35S_control_Erb1_merged_reads.sgr", comment.char="#",col.names=c("chromosome","position","35S_control_1_trimmed_plus_strand_startpositions","35S_control_2_trimmed_plus_strand_startpositions","35S_1M7_1_trimmed_plus_strand_startpositions","35S_1M7_2_trimmed_plus_strand_startpositions"))
mergedstartsmut <- read.table("35S_control_Erb1_merged_dropoffcounts.sgr",comment.char="#",col.names=c("chromosome","position","35S_control_1_trimmed_plus_strand_startpositions","35S_control_2_trimmed_plus_strand_startpositions","35S_1M7_1_trimmed_plus_strand_startpositions","35S_1M7_2_trimmed_plus_strand_startpositions"))

setwd("./../../")
#### ------- DATA PRE-PROCESSING (CALCULATION LOG RATIOS OF MUTATION RATES AND EMPIRICAL P-VALUES ) ------- ######

## Calculating dropoff rate ratios
logdropoffswt <- calculateLDRs(mergedcountswt,mergedstartswt,noreplicates, refsequence)  
logdropoffsmut <- calculateLDRs(mergedcountsmut,mergedstartsmut,noreplicates, refsequence)


### ------- QUALITY CONTROL: INSPECTION OF LOG DROP OFF RATE RATIOS------ ######


setwd("Analysis/LMR_and_LDR_plots/35S")

pdf('LDR-35S-Control1-Control2-comparison_delta5.pdf',width=6,height=4,paper='special')
hist(logdropoffswt$LDR_C, breaks = 30, main = 'Null distribution of LDRs - delta 5')
dev.off()


pdf('LDR-35S-Control1-Control2-comparison_erb1.pdf',width=6,height=4,paper='special')
hist(logdropoffsmut$LDR_C, breaks = 30, main = 'Null distribution of LDRs - erb1')
dev.off()

ldr_ct_delta5 <-- logdropoffswt$LDR_CT
ldr_ct_erb1 <-- logdropoffsmut$LDR_CT

pdf('LDR-35S-Treatment1-Control1-comparison_delta5.pdf',width=6,height=4,paper='special')
hist(ldr_ct_delta5[ , 1:1], breaks = 30, main = 'LDR T1 - C1 distribution delta 5')
dev.off()

pdf('LDR-35S-Treatment1-Control2-comparison_delta5.pdf',width=6,height=4,paper='special')
hist(ldr_ct_delta5[ , 2:2], breaks = 30, main = 'LDR T1 - C2 distribution delta 5')
dev.off()

pdf('LDR-35S-Treatment2-Control1-comparison_delta5.pdf',width=6,height=4,paper='special')
hist(ldr_ct_delta5[ , 3:3], breaks = 30, main = 'LDR T2 - C1 distribution delta 5')
dev.off()

pdf('LDR-35S-Treatment2-Control2-comparison_delta5.pdf',width=6,height=4,paper='special')
hist(ldr_ct_delta5[ , 4:4], breaks = 30, main = 'LDR T2 - C2 distribution delta 5')
dev.off()


pdf('LDR-35S-Treatment1-Control1-comparison_erb1.pdf',width=6,height=4,paper='special')
hist(ldr_ct_erb1[ , 1:1], breaks = 30, main = 'LDR T1 - C1 distribution erb1')
dev.off()

pdf('LDR-35S-Treatment1-Control2-comparison_erb1.pdf',width=6,height=4,paper='special')
hist(ldr_ct_erb1[ , 2:2], breaks = 30, main = 'LDR T1 - C2 distribution erb1')
dev.off()

pdf('LDR-35S-Treatment2-Control1-comparison_erb1.pdf',width=6,height=4,paper='special')
hist(ldr_ct_erb1[ , 3:3], breaks = 30, main = 'LDR T2 - C1 distribution erb1')
dev.off()

pdf('LDR-35S-Treatment2-Control2-comparison_erb1.pdf',width=6,height=4,paper='special')
hist(ldr_ct_erb1[ , 4:4], breaks = 30, main = 'LDR T2 - C2 distribution erb1')
dev.off()

setwd('./../../..')


## ------------------------------------------------------------------------##

Nc <- Nt <- noreplicates

strand = "+"

##Calculation of empirical p values for the two sets of samples: the WT and mutant.

empPvals_1 <- computePvals(logdropoffswt$LDR_C,logdropoffswt$LDR_CT, Nc, Nt, strand, logdropoffswt$nuclPosition,
                           logdropoffswt$nuclSelection$analysedC, logdropoffswt$nuclSelection$analysedCT)

empPvals_2 <- computePvals(logdropoffsmut$LDR_C,logdropoffsmut$LDR_CT, Nc, Nt, strand, logdropoffsmut$nuclPosition,
                           logdropoffsmut$nuclSelection$analysedC, logdropoffsmut$nuclSelection$analysedCT)


#stretches contain the selection of nulceotide positions that have valid log drop off rate ratios
#overlapping the stretches of nucleotides selected in each group of samples to get a representative set
stretches <-overlapsRanges(logdropoffswt$stretches,logdropoffsmut$stretches)


####### ------- HMM ANALYSIS (CALCULATION OF POSTERIOR PROBABILITIES FOR THE FOUR HIDDEN STATES) ------- ########
###computes posterior probabilities of all nucleotides in the stretches specified above.
#The most important arguments are the LDRs, the positions used to compute the
#null distribution, as well as the positions of the selected stretches
#and nucleotide pairs where the LDRs were obtained.

##Before computing posterior probabilities, we set up matrices to hold the required p-values
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


#### ------- QUALITY CONTROL: INSPECTION OF P-VALUES ------ ######

setwd("Analysis/pvalues_plots/35S")

write.table(t(Pv1), file="pvalues_delta5.txt", row.names=TRUE, col.names=TRUE)
write.table(t(Pv2), file="pvalues_deltaerb1.txt", row.names=TRUE, col.names=TRUE)


pdf('pvalues-35S-Treatment1-Control1-comparison_delta5.pdf',width=6,height=4,paper='special')
hist(Pv1[1:1,], breaks = 30, main = 'pvalues T1 - C1 distribution delta 5')
dev.off()

pdf('pvalues-35S-Treatment1-Control2-comparison_delta5.pdf',width=6,height=4,paper='special')
hist(Pv1[2:2,], breaks = 30, main = 'pvalues T1 - C2 distribution delta 5')
dev.off()

pdf('pvalues-35S-Treatment2-Control1-comparison_delta5.pdf',width=6,height=4,paper='special')
hist(Pv1[3:3,], breaks = 30, main = 'pvalues T2 - C1 distribution delta 5')
dev.off()

pdf('pvalues-35S-Treatment2-Control2-comparison_delta5.pdf',width=6,height=4,paper='special')
hist(Pv1[4:4,], breaks = 30, main = 'pvalues T2 - C2 distribution delta 5')
dev.off()



pdf('pvalues-35S-Treatment1-Control1-comparison_erb1.pdf',width=6,height=4,paper='special')
hist(Pv2[1:1,], breaks = 30, main = 'pvalues T1 - C1 distribution erb1')
dev.off()

pdf('pvalues-35S-Treatment1-Control2-comparison_erb1.pdf',width=6,height=4,paper='special')
hist(Pv2[2:2,], breaks = 30, main = 'pvalues T1 - C2 distribution erb1')
dev.off()

pdf('pvalues-35S-Treatment2-Control1-comparison_erb1.pdf',width=6,height=4,paper='special')
hist(Pv2[3:3,], breaks = 30, main = 'pvalues T2 - C1 distribution erb1')
dev.off()

pdf('pvalues-35S-Treatment2-Control2-comparison_erb1.pdf',width=6,height=4,paper='special')
hist(Pv2[4:4,], breaks = 30, main = 'pvalues T2 - C2 distribution erb1')
dev.off()

setwd('./../../..')



