#### ------- PACKAGES INSTALLATION AND IMPORT OF HELPER FUNCTIONS ------ ######

# This sripts assumes: R version 3.6.3 (2020-02-29); RStudio Version 1.1.442

install.packages("rstudioapi")
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

install.packages("BiocManager")

BiocManager::install(c("Biostrings", "SummarizedExperiment"), version = "3.10")

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
setwd('..')

outputfilename <-paste0('35S','_diffBUM_HMM_WT_vs_Erb1_diff_BUM_HMM_analysed','.txt')

mergedcountswt <- read.table("Data/35S_data/35S_control_delta5_merged_dropoffcounts.sgr", comment.char="#",col.names=c("chromosome","position","35S_control_1_trimmed_plus_strand_startpositions","35S_control_2_trimmed_plus_strand_startpositions","35S_1M7_1_trimmed_plus_strand_startpositions","35S_1M7_2_trimmed_plus_strand_startpositions"))
mergedstartswt <- read.table("Data/35S_data/35S_control_delta5_merged_reads.sgr",comment.char="#",col.names=c("chromosome","position","35S_control_1_trimmed_plus_strand_startpositions","35S_control_2_trimmed_plus_strand_startpositions","35S_1M7_1_trimmed_plus_strand_startpositions","35S_1M7_2_trimmed_plus_strand_startpositions"))
mergedcountsmut <- read.table("Data/35S_data/35S_control_Erb1_merged_dropoffcounts.sgr", comment.char="#",col.names=c("chromosome","position","35S_control_1_trimmed_plus_strand_startpositions","35S_control_2_trimmed_plus_strand_startpositions","35S_1M7_1_trimmed_plus_strand_startpositions","35S_1M7_2_trimmed_plus_strand_startpositions"))
mergedstartsmut <- read.table("Data/35S_data/35S_control_Erb1_merged_reads.sgr",comment.char="#",col.names=c("chromosome","position","35S_control_1_trimmed_plus_strand_startpositions","35S_control_2_trimmed_plus_strand_startpositions","35S_1M7_1_trimmed_plus_strand_startpositions","35S_1M7_2_trimmed_plus_strand_startpositions"))

head(mergedcountswt)

## Calculating dropoff rates
logdropoffswt <- calculateLDRs(mergedcountswt,mergedstartswt,noreplicates, refsequence)  
logdropoffsmut <- calculateLDRs(mergedcountsmut,mergedstartsmut,noreplicates, refsequence)


#### ------- QUALITY CONTROL: INSPECTION OF LOG MUTATION RATES (OPTIONAL)------ ######
#TO RUN,UNCOMMENT LINES 58-107 WITH: CTRL + ALT + C On Windows or command + SHIFT + C in Mac OS 

#setwd("Analysis/LMR_and_LDR_plots/35S")

# pdf('LDR-35S-Control1-Control2-comparison_delta5.pdf',width=6,height=4,paper='special')
# hist(logdropoffswt$LDR_C, breaks = 30, main = 'Null distribution of LDRs - delta 5')
# dev.off()
# 
# 
# pdf('LDR-35S-Control1-Control2-comparison_erb1.pdf',width=6,height=4,paper='special')
# hist(logdropoffsmut$LDR_C, breaks = 30, main = 'Null distribution of LDRs - erb1')
# dev.off()
# 
# ldr_ct_delta5 <-- logdropoffswt$LDR_CT
# ldr_ct_erb1 <-- logdropoffsmut$LDR_CT
# 
# pdf('LDR-35S-Treatment1-Control1-comparison_delta5.pdf',width=6,height=4,paper='special')
# hist(ldr_ct_delta5[ , 1:1], breaks = 30, main = 'LDR T1 - C1 distribution delta 5')
# dev.off()
# 
# pdf('LDR-35S-Treatment1-Control2-comparison_delta5.pdf',width=6,height=4,paper='special')
# hist(ldr_ct_delta5[ , 2:2], breaks = 30, main = 'LDR T1 - C2 distribution delta 5')
# dev.off()
# 
# pdf('LDR-35S-Treatment2-Control1-comparison_delta5.pdf',width=6,height=4,paper='special')
# hist(ldr_ct_delta5[ , 3:3], breaks = 30, main = 'LDR T2 - C1 distribution delta 5')
# dev.off()
# 
# pdf('LDR-35S-Treatment2-Control2-comparison_delta5.pdf',width=6,height=4,paper='special')
# hist(ldr_ct_delta5[ , 4:4], breaks = 30, main = 'LDR T2 - C2 distribution delta 5')
# dev.off()
# 
# 
# pdf('LDR-35S-Treatment1-Control1-comparison_erb1.pdf',width=6,height=4,paper='special')
# hist(ldr_ct_erb1[ , 1:1], breaks = 30, main = 'LDR T1 - C1 distribution erb1')
# dev.off()
# 
# pdf('LDR-35S-Treatment1-Control2-comparison_erb1.pdf',width=6,height=4,paper='special')
# hist(ldr_ct_erb1[ , 2:2], breaks = 30, main = 'LDR T1 - C2 distribution erb1')
# dev.off()
# 
# pdf('LDR-35S-Treatment2-Control1-comparison_erb1.pdf',width=6,height=4,paper='special')
# hist(ldr_ct_erb1[ , 3:3], breaks = 30, main = 'LDR T2 - C1 distribution erb1')
# dev.off()
# 
# pdf('LDR-35S-Treatment2-Control2-comparison_erb1.pdf',width=6,height=4,paper='special')
# hist(ldr_ct_erb1[ , 4:4], breaks = 30, main = 'LDR T2 - C2 distribution erb1')
# dev.off()
#
#setwd('..')
#setwd('..')
#setwd('..')


Nc <- Nt <- noreplicates

strand = "+"


#logdropoffswt <- logdropoffs_incell
#logdropoffsmut <- logdropoffs_exvivo

###
empPvals_1 <- computePvals(logdropoffswt$LDR_C,logdropoffswt$LDR_CT, Nc, Nt, strand, logdropoffswt$nuclPosition,
                           logdropoffswt$nuclSelection$analysedC, logdropoffswt$nuclSelection$analysedCT)

empPvals_2 <- computePvals(logdropoffsmut$LDR_C,logdropoffsmut$LDR_CT, Nc, Nt, strand, logdropoffsmut$nuclPosition,
                           logdropoffsmut$nuclSelection$analysedC, logdropoffsmut$nuclSelection$analysedCT)


stretches <-overlapsRanges(logdropoffswt$stretches,logdropoffsmut$stretches)

## Number of nucleotides in the sequence = number of rows in empPvals_1
nNucl <- length(empPvals_1[1, ])


## ------------------------------------------------------------------------
###computes posterior probabilities of all nucleotides in the stretches specified above.
#This is the step where the null distributions and p-values are calculated as well.
#The most important arguments are the LDRs, the positions used to compute the
#null distribution, as well as the positional information of the selected stretches
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


#### ------- QUALITY CONTROL: INSPECTION OF P-VALUES (OPTIONAL)------ ######
#TO RUN,UNCOMMENT LINES 206-246 WITH: CTRL + ALT + C On Windows or command + SHIFT + C in Mac OS 

#setwd("Analysis/pvalues_plots/35S")

#write.table(Pv1, file="pvalues_delta5.txt", row.names=TRUE, col.names=TRUE)
#write.table(Pv2, file="pvalues_deltaerb1.txt", row.names=TRUE, col.names=TRUE)

# 
# pdf('pvalues-35S-Treatment1-Control1-comparison_delta5.pdf',width=6,height=4,paper='special')
# hist(Pv1[1:1,], breaks = 30, main = 'pvalues T1 - C1 distribution delta 5')
# dev.off()
# 
# pdf('pvalues-35S-Treatment1-Control2-comparison_delta5.pdf',width=6,height=4,paper='special')
# hist(Pv1[2:2,], breaks = 30, main = 'pvalues T1 - C2 distribution delta 5')
# dev.off()
# 
# pdf('pvalues-35S-Treatment2-Control1-comparison_delta5.pdf',width=6,height=4,paper='special')
# hist(Pv1[3:3,], breaks = 30, main = 'pvalues T2 - C1 distribution delta 5')
# dev.off()
# 
# pdf('pvalues-35S-Treatment2-Control2-comparison_delta5.pdf',width=6,height=4,paper='special')
# hist(Pv1[4:4,], breaks = 30, main = 'pvalues T2 - C2 distribution delta 5')
# dev.off()
# 
# 
# 
# pdf('pvalues-35S-Treatment1-Control1-comparison_erb1.pdf',width=6,height=4,paper='special')
# hist(Pv2[1:1,], breaks = 30, main = 'pvalues T1 - C1 distribution erb1')
# dev.off()
# 
# pdf('pvalues-35S-Treatment1-Control2-comparison_erb1.pdf',width=6,height=4,paper='special')
# hist(Pv2[2:2,], breaks = 30, main = 'pvalues T1 - C2 distribution erb1')
# dev.off()
# 
# pdf('pvalues-35S-Treatment2-Control1-comparison_erb1.pdf',width=6,height=4,paper='special')
# hist(Pv2[3:3,], breaks = 30, main = 'pvalues T2 - C1 distribution erb1')
# dev.off()
# 
# pdf('pvalues-xist-Treatment2-Control2-comparison_erb1.pdf',width=6,height=4,paper='special')
# hist(Pv2[4:4,], breaks = 30, main = 'pvalues T2 - C2 distribution erb1')
# dev.off()
#
#
#setwd('..')
#setwd('..')
#setwd('..')

posteriors_diff <- hmmFwbw_differential_two_betas(pvaluesstretch)
colnames(posteriors_diff) <- c("UU","UM","MU","MM")
head(posteriors_diff)

## All posterior probabilities need to be shifted by 1 nt because the RT
## is believed to stop 1 nucleotide before the modified nucleodide.
## So below a new matrix is made containing the values shifted by one position.
## ------------------------------------------------------------------------
shifted_posteriors <- matrix(, nrow=dim(posteriors_diff)[1], ncol=4)
shifted_posteriors[1:(length(shifted_posteriors[,1]) - 1), ] <- posteriors_diff[2:(length(shifted_posteriors[,1])), ]
colnames(shifted_posteriors) <- c("UU","UM","MU","MM")


head(shifted_posteriors)
head(posteriors_diff)


differentiallymod <- shifted_posteriors[,2] + shifted_posteriors[,3]

## ------------------------------------------------------------------------

setwd("Analysis/diffBUM-HMM")

png("35S_sum_of_diff_states_diff_BUM_HMM_.png")
plot(differentiallymod, xlab = 'Nucleotide position',
     ylab = 'Probability of modification (UM+MU)',
     main = 'diffBUMHMM output: ProbabilITY of differential modification between delta5 and erb1',
     ylim = c(0,1))
dev.off()


shifted_posteriors <- replace(shifted_posteriors,is.na(shifted_posteriors),-999)
write.table(shifted_posteriors,sep="\t",quote=FALSE,file=outputfilename,col.names = c("UU","UM","MU","MM"), row.names = TRUE)

