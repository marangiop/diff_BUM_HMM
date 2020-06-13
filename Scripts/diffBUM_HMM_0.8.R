#USER HAS TO ALWAYS MANUALLY SET THE WORKING DIRECTORY TO THE CLONED DIFFBUM-HMM FOLDER
#ON RSTUDIO BEFORE RUNNING THE SCRIPT 

#setwd("C://Users/User/Desktop/diff_BUM_HMM/")



##STEP FOR CORRECT EXECUTION
#FIRST STEP: OPEN R STUDIO, GO TO THE SCRIPTS FOLDER AND SET THAT AS WORKING DIRECTORY, THEN RUN THE REST OF THE PIPELINE 
setwd(dirname(getwd()))


#library(Rmpfr)

#IF YOU HAVEN'T DONE IT YET, READ THE FOLLOWING FILE AND FOLLOW ITS INSTRUCTIONS:
# Bioconductor_March2020_InstallationBug_Fixed.txt

#Line 11 is commented out because after done what it says on line 5 
#(i.e. updating R and Bioconductor to versions 3.6.3 and 3.1.0, respectively) 
#there is no need to run line 8 every time the script is run. 
#source("https://bioconductor.org/biocLite.R")
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

#getwd()

noreplicates <- 2

#cat(working_directory)

#ref_seq_directory <- paste(working_directory, "Reference sequences/Xist.seq",molecule  ,sep="/")
refsequence <- "Xist.seq"
#cat(refsequence)

outputfilename <-paste0('Xist in vivo vs. ex vivo scaled by dc_toby','_diff_BUM_HMM_analysis','.txt')

table1_incell <- read.delim("Data/XIST_1M7_in-cell_rep1.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","in_cell_DMSO1_read_count","in_cell_DMSO1_mutation_rate","in_cell_1M71_read_count","in_cell_1M71_mutation_rate"))
table2_incell <- read.delim("Data/XIST_1M7_in-cell_rep2.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","in_cell_DMSO2_read_count","in_cell_DMSO2_mutation_rate","in_cell_1M72_read_count","in_cell_1M72_mutation_rate"))

table1_exvivo <- read.delim("Data/XIST_1M7_ex-vivo_rep1.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","ex_vivo_DMSO1_read_count","ex_vivo_DMSO1_mutation_rate","ex_vivo_1M71_read_count","ex_vivo_1M71_mutation_rate"))
table2_exvivo <- read.delim("Data/XIST_1M7_ex-vivo_rep2.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","ex_vivo_DMSO2_read_count","ex_vivo_DMSO2_mutation_rate","ex_vivo_1M72_read_count","ex_vivo_1M72_mutation_rate"))

head(table1_incell["in_cell_DMSO1_read_count"])

dc_incell <- read.delim("Old/Xist_1M7_in-cell_wDC.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","in_cell_DMSO1_read_count","in_cell_DMSO1_mutation_rate","in_cell_1M71_read_count","in_cell_1M71_mutation_rate","DC_read_count" ,"DC_mutation_rate"))
dc_exvivo <- read.delim("Old/XIST_1M7_ex-vivo_wDC.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","ex_vivo_DMSO1_read_count","ex_vivo_DMSO1_mutation_rate","ex_vivo_1M71_read_count","ex_vivo_1M71_mutation_rate","DC_read_count" ,"DC_mutation_rate"))



#READ COUNTS  i.e COVERAGE
incell_counts <- data.frame("in_cell_DMSO1_read_count" = table1_incell["in_cell_DMSO1_read_count"],"in_cell_DMSO2_read_count" = table2_incell["in_cell_DMSO2_read_count"],"X1M7_1_read_count" = table1_incell["in_cell_1M71_read_count"],"X1M7_2_read_count"= table2_incell["in_cell_1M72_read_count"])
exvivo_counts <- data.frame("ex_vivo_DMSO1_read_count" = table1_exvivo["ex_vivo_DMSO1_read_count"],"ex_vivo_DMSO2_read_count" = table2_exvivo["ex_vivo_DMSO2_read_count"],"X1M7_1_read_count" = table1_exvivo["ex_vivo_1M71_read_count"],"X1M7_2_read_count"= table2_exvivo["ex_vivo_1M72_read_count"])

head(incell_counts)
head(exvivo_counts)

#MUTATION RATES
incell_rates <- data.frame("in_cell_DMSO1_mutation_rate" = table1_incell["in_cell_DMSO1_mutation_rate"],"in_cell_DMSO2_mutation_rate" = table2_incell["in_cell_DMSO2_mutation_rate"],"X1M7_1_mutation_rate" = table1_incell["in_cell_1M71_mutation_rate"],"X1M7_2_mutation_rate"= table2_incell["in_cell_1M72_mutation_rate"])
exvivo_rates <- data.frame("ex_vivo_DMSO1_mutation_rate" = table1_exvivo["ex_vivo_DMSO1_mutation_rate"],"ex_vivo_DMSO2_mutation_rate" = table2_exvivo["ex_vivo_DMSO2_mutation_rate"],"X1M7_1_mutation_rate" = table1_exvivo["ex_vivo_1M71_mutation_rate"],"X1M7_2_mutation_rate"= table2_exvivo["ex_vivo_1M72_mutation_rate"])

# DENATURED COTNROLS 
dc_incell_column <- data.frame("DC_mutation_rate"=dc_incell["DC_mutation_rate"] )
dc_exvivo_column <- data.frame("DC_mutation_rate"=dc_exvivo["DC_mutation_rate"])


#$in_cell_DMSO1_mutation_rate <- mpfr(incell_rates$in_cell_DMSO1_mutation_rate, 200)
#incell_rates$in_cell_DMSO2_mutation_rate <- mpfr(incell_rates$in_cell_DMSO2_mutation_rate, 200)
#incell_rates$in_cell_1M71_mutation_rate <- mpfr(incell_rates$in_cell_1M71_mutation_rate, 200)
#incell_rates$in_cell_1M72_mutation_rate <- mpfr(incell_rates$in_cell_1M72_mutation_rate, 200)

library(formattable)
incell_rates = formattable(incell_rates,digits = 8, format = "f" )
exvivo_rates = formattable(exvivo_rates,digits = 8, format = "f" )

dc_incell_column = formattable(dc_incell_column, digits = 8, format = "f" )
dc_exvivo_column = formattable(dc_exvivo_column,digits = 8, format = "f" )

#incell_rates_vector <- as.vector(incell_rates[['in_cell_DMSO1_mutation_rate']])
#y = formatC(incell_rates_vector, digits = 8, format = "f")
#y = formattable(incell_rates_vector,digits = 8, format = "f" )

#dc_incell_column_vector <- as.vector(dc_incell_column[['DC_mutation_rate']])
#w = formattable(dc_incell_column_vector,digits = 8, format = "f" )
#
#a = y/w

#form_a = formattable(a,digits = 8, format = "f" )
#z <- mpfr(incell_rates_vector, precBits = 10)

#dc_incell_column$DC_mutation_rate <- mpfr(dc_incell_column$DC_mutation_rate,200)
#dc_exvivo_column$DC_mutation_rate <- mpfr(dc_exvivo_column$DC_mutation_rate,200)

#Dividing DMSO and 1MT columns for both in cell and ex vivo by respective DC replicate
#Note:only 1 single replicate for DC was available, for each condition. 
scaled_incell_rates <- incell_rates[,] / dc_incell_column
scaled_exvivo_rates <- exvivo_rates[,] / dc_exvivo_column

#Removing NaN and Inf values
is.na(scaled_incell_rates)<-sapply(scaled_incell_rates, is.infinite)
scaled_incell_rates[is.na(scaled_incell_rates)]<-0

is.na(scaled_exvivo_rates)<-sapply(scaled_exvivo_rates, is.infinite)
scaled_exvivo_rates[is.na(scaled_exvivo_rates)]<-0

#scaled_incell_rates[,] <- sapply(scaled_incell_rates[,], as.numeric)
#scaled_exvivo_rates[,] <- sapply(scaled_exvivo_rates[,], as.numeric)


#WORKING
attributes(scaled_incell_rates) <- NULL
attributes(scaled_exvivo_rates) <- NULL

#df <-data.frame(t(sapply(scaled_incell_rates,c)))
#df_2 <-data.frame(t(sapply(scaled_exvivo_rates,c)))

n <- length(scaled_incell_rates[[1]])
scaled_incell_rates_df  <- structure(scaled_incell_rates,  row.names = c(NA, -n), class = "data.frame")
colnames(scaled_incell_rates_df) <- c("in_cell_DMSO1_mutation_rate", "in_cell_DMSO2_mutation_rate", "in_cell_1M71_mutation_rate", "in_cell_1M72_mutation_rate" )
scaled_exvivo_rates_df  <- structure(scaled_exvivo_rates,  row.names = c(NA, -n), class = "data.frame")
colnames(scaled_exvivo_rates_df) <- c("in_cell_DMSO1_mutation_rate", "in_cell_DMSO2_mutation_rate", "in_cell_1M71_mutation_rate", "in_cell_1M72_mutation_rate" )

#Setting regions to 0, based on sanity check

scaled_incell_rates_df[2500:4500,]=0 
scaled_exvivo_rates_df[2500:4500,]=0 

scaled_incell_rates_df[1:78,]=0 
scaled_exvivo_rates_df[1:78,]=0 

scaled_incell_rates_df[2451:2599,]=0 
scaled_exvivo_rates_df[2451:2599,]=0 

scaled_incell_rates_df[17801:17918,]=0 
scaled_exvivo_rates_df[17801:17918,]=0 


#scaled_exvivo_rates_df  <- structure(scaled_exvivo_rates, col.names = c("ex_vivo_DMSO1_mutation_rate", "ex_vivo_DMSO2_mutation_rate", "ex_vivo_1M71_mutation_rate", "ex_vivo_1M72_mutation_rate" ), row.names = c(NA, -n), class = "data.frame")



#scaled_incell_rates_df <- data.frame(matrix(unlist(scaled_incell_rates), byrow=FALSE, ncol=length(scaled_incell_rates[[1]])),stringsAsFactors=FALSE)
#scaled_incell_rates_df <- data.frame(matrix(unlist(scaled_exvivo_rates)),stringsAsFactors=FALSE)


#Option 1: Leave mutation counts untouched, using non-scaled rated
#mutation_counts_in_cell <-  incell_counts * incell_rates
#mutation_counts_ex_vivo <-  exvivo_counts * exvivo_rates

#Option 2: calculate counts using scaled rates 
mutation_counts_in_cell <-  incell_counts * scaled_incell_rates_df
mutation_counts_ex_vivo <-  exvivo_counts * scaled_exvivo_rates_df

mutation_counts_in_cell <- cbind(position=0,mutation_counts_in_cell)
mutation_counts_ex_vivo <- cbind(position=0,mutation_counts_ex_vivo)
mutation_counts_in_cell <- cbind(gene="Xist",mutation_counts_in_cell)
mutation_counts_ex_vivo <- cbind(gene="Xist",mutation_counts_ex_vivo)

incell_counts <- cbind(position=0,incell_counts)
exvivo_counts <- cbind(position=0,exvivo_counts)
incell_counts <- cbind(gene="Xist",incell_counts)
exvivo_counts <- cbind(gene="Xist",exvivo_counts)

head(mutation_counts_ex_vivo)
head(exvivo_counts)



## Calculating dropoff rates
logdropoffs_incell <- calculateLDRs(incell_counts,mutation_counts_in_cell, noreplicates, refsequence)
logdropoffs_exvivo <- calculateLDRs(exvivo_counts,mutation_counts_ex_vivo, noreplicates, refsequence)


pdf('LMR-xist-Control1-Control2-comparison_invivo.pdf',width=6,height=4,paper='special')
hist(logdropoffs_incell$LDR_C, breaks = 30, main = 'Null distribution of LDRs - in vivo')
dev.off()


pdf('LMR-xist-Control1-Control2-comparison_exvivo.pdf',width=6,height=4,paper='special')
hist(logdropoffs_exvivo$LDR_C, breaks = 30, main = 'Null distribution of LDRs - ex vivo')
dev.off()

ldr_ct_invivo <-- logdropoffs_incell$LDR_CT
ldr_ct_exvivo <-- logdropoffs_exvivo$LDR_CT

pdf('LMR-xist-Treatment1-Control1-comparison_invivo.pdf',width=6,height=4,paper='special')
hist(ldr_ct_invivo[ , 1:1], breaks = 30, main = 'LDR T1 - C1 distribution in vivo')
dev.off()

pdf('LMR-xist-Treatment1-Control2-comparison_invivo.pdf',width=6,height=4,paper='special')
hist(ldr_ct_invivo[ , 2:2], breaks = 30, main = 'LDR T1 - C2 distribution in vivo')
dev.off()

pdf('LMR-xist-Treatment2-Control1-comparison_invivo.pdf',width=6,height=4,paper='special')
hist(ldr_ct_invivo[ , 3:3], breaks = 30, main = 'LDR T2 - C1 distribution in vivo')
dev.off()

pdf('LMR-xist-Treatment2-Control2-comparison_invivo.pdf',width=6,height=4,paper='special')
hist(ldr_ct_invivo[ , 4:4], breaks = 30, main = 'LDR T2 - C2 distribution in vivo')
dev.off()


pdf('LMR-xist-Treatment1-Control1-comparison_exvivo.pdf',width=6,height=4,paper='special')
hist(ldr_ct_exvivo[ , 1:1], breaks = 30, main = 'LDR T1 - C1 distribution ex vivo')
dev.off()

pdf('LMR-xist-Treatment1-Control2-comparison_exvivo.pdf',width=6,height=4,paper='special')
hist(ldr_ct_exvivo[ , 2:2], breaks = 30, main = 'LDR T1 - C2 distribution ex vivo')
dev.off()

pdf('LMR-xist-Treatment2-Control1-comparison_exvivo.pdf',width=6,height=4,paper='special')
hist(ldr_ct_exvivo[ , 3:3], breaks = 30, main = 'LDR T2 - C1 distribution in vivo')
dev.off()

pdf('LMR-xist-Treatment2-Control2-comparison_exvivo.pdf',width=6,height=4,paper='special')
hist(ldr_ct_exvivo[ , 4:4], breaks = 30, main = 'LDR T2 - C2 distribution in vivo')
dev.off()


## ------------------------------------------------------------------------
###check if the matrices of p-values can be called after the pipeline has been run twice
#head(logdropoffswt$LDR_C)
#head(logdropoffswt$LDR_CT)
#head(logdropoffsmut$LDR_C)
#head(logdropoffsmut$LDR_CT)

head(logdropoffs_incell$LDR_C)
head(logdropoffs_incell$LDR_CT)
head(logdropoffs_exvivo$LDR_C)
head(logdropoffs_exvivo$LDR_CT)

Nc <- Nt <- noreplicates

strand = "+"


#logdropoffswt <- logdropoffs_incell
#logdropoffsmut <- logdropoffs_exvivo

###
empPvals_1 <- computePvals(logdropoffs_incell$LDR_C,logdropoffs_incell$LDR_CT, Nc, Nt, strand, logdropoffs_incell$nuclPosition,
                           logdropoffs_incell$nuclSelection$analysedC, logdropoffs_incell$nuclSelection$analysedCT)

empPvals_2 <- computePvals(logdropoffs_exvivo$LDR_C,logdropoffs_exvivo$LDR_CT, Nc, Nt, strand, logdropoffs_exvivo$nuclPosition,
                           logdropoffs_exvivo$nuclSelection$analysedC, logdropoffs_exvivo$nuclSelection$analysedCT)


stretches <-overlapsRanges(logdropoffs_incell$stretches,logdropoffs_exvivo$stretches)

s## Number of nucleotides in the sequence = number of rows in empPvals_1
nNucl <- length(empPvals_1[1, ])


## ------------------------------------------------------------------------
###computes posterior probabilities of all nucleotides in the stretches specified above.
#This is the step where the null distributions and p-values are calculated as well.
#The most important arguments are the LDRs, the positions used to compute the
#null distribution, as well as the positional information of the selected stretches
#and nucleotide pairs where the LDRs were obtained.

##Before computing posterior probabilities, we set up matrices to hold the required p-values

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

write.table(Pv1, file="pvalues_invivo-latest-version-after-toby-bugfix.txt", row.names=TRUE, col.names=TRUE)
write.table(Pv2, file="pvalues_exvivo-latest-version-after-toby-bugfix.txt", row.names=TRUE, col.names=TRUE)


pdf('pvalues-xist-Treatment1-Control1-comparison_invivo.pdf',width=6,height=4,paper='special')
hist(Pv1[1:1,], breaks = 30, main = 'pvalues T1 - C1 distribution in vivo')
dev.off()

pdf('pvalues-xist-Treatment1-Control2-comparison_invivo.pdf',width=6,height=4,paper='special')
hist(Pv1[2:2,], breaks = 30, main = 'pvalues T1 - C2 distribution in vivo')
dev.off()

pdf('pvalues-xist-Treatment2-Control1-comparison_invivo.pdf',width=6,height=4,paper='special')
hist(Pv1[3:3,], breaks = 30, main = 'pvalues T2 - C1 distribution in vivo')
dev.off()

pdf('pvalues-xist-Treatment2-Control2-comparison_invivo.pdf',width=6,height=4,paper='special')
hist(Pv1[4:4,], breaks = 30, main = 'pvalues T2 - C2 distribution in vivo')
dev.off()



pdf('pvalues-xist-Treatment1-Control1-comparison_exvivo.pdf',width=6,height=4,paper='special')
hist(Pv2[1:1,], breaks = 30, main = 'pvalues T1 - C1 distribution ex vivo')
dev.off()

pdf('pvalues-xist-Treatment1-Control2-comparison_exvivo.pdf',width=6,height=4,paper='special')
hist(Pv2[2:2,], breaks = 30, main = 'pvalues T1 - C2 distribution ex vivo')
dev.off()

pdf('pvalues-xist-Treatment2-Control1-comparison_exvivo.pdf',width=6,height=4,paper='special')
hist(Pv2[3:3,], breaks = 30, main = 'pvalues T2 - C1 distribution ex vivo')
dev.off()

pdf('pvalues-xist-Treatment2-Control2-comparison_exvivo.pdf',width=6,height=4,paper='special')
hist(Pv2[4:4,], breaks = 30, main = 'pvalues T2 - C2 distribution ex vivo')
dev.off()



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


differentiallymod <- shifted_posteriors[,2] + shifted_posteriors[,3]

## ------------------------------------------------------------------------
png("Xist_sum_of_diff_states_diff_BUM_HMM_scaled_w_dc.png")
plot(differentiallymod, xlab = 'Nucleotide position',
     ylab = 'Probability of modification (UM+MU)',
     main = 'diffBUMHMM output: ProbabilITY of differential modification between in vivo and ex vivo',
     ylim = c(0,1))
dev.off()
## ----eval=FALSE----------------------------------------------------------
## ## Call the function with the additonal tolerance parameter
## posteriors <- computeProbs(LDR_C, LDR_CT, Nc, Nt, '+', nuclPosition,
##                            nuclSelection$analysedC, nuclSelection$analysedCT,
##                            stretches, 0.001)

## ------------------------------------------------------------------------
shifted_posteriors <- replace(shifted_posteriors,is.na(shifted_posteriors),-999)
write.table(shifted_posteriors,sep="\t",quote=FALSE,file=outputfilename,col.names = c("UU","UM","MU","MM"), row.names = TRUE)

