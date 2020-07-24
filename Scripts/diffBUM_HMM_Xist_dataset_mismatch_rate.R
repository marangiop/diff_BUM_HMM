#### ------- PACKAGES INSTALLATION AND IMPORT OF HELPER FUNCTIONS ------ ######

# This script assumes: R version 4.0.0 (2020-04-24); RStudio Version 1.2.5001

install.packages("rstudioapi")
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

install.packages("BiocManager")

BiocManager::install(c("Biostrings", "SummarizedExperiment"), version = "3.11")

setwd('..')

library(formattable)
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
refsequence <- "Xist.seq"
setwd('..')


outputfilename <-paste0('Xist','_diff_BUM_HMM_analysis_withscaling_gaussiannoiseadded','.txt')

#ADDED CODE FOR PROCESSING NEW DATA 


table1_incell <- read.delim("Data/Xist_dataset/XIST_1M7_in-cell_rep1.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","in_cell_DMSO1_read_count","in_cell_DMSO1_mutation_rate","in_cell_1M71_read_count","in_cell_1M71_mutation_rate"))
table2_incell <- read.delim("Data/Xist_dataset/XIST_1M7_in-cell_rep2.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","in_cell_DMSO2_read_count","in_cell_DMSO2_mutation_rate","in_cell_1M72_read_count","in_cell_1M72_mutation_rate"))

table1_exvivo <- read.delim("Data/Xist_dataset/XIST_1M7_ex-vivo_rep1.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","ex_vivo_DMSO1_read_count","ex_vivo_DMSO1_mutation_rate","ex_vivo_1M71_read_count","ex_vivo_1M71_mutation_rate"))
table2_exvivo <- read.delim("Data/Xist_dataset/XIST_1M7_ex-vivo_rep2.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","ex_vivo_DMSO2_read_count","ex_vivo_DMSO2_mutation_rate","ex_vivo_1M72_read_count","ex_vivo_1M72_mutation_rate"))

head(table1_incell["in_cell_DMSO1_read_count"])

dc_incell <- read.delim("Data/Xist_dataset/Xist_1M7_in-cell_wDC.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","in_cell_DMSO1_read_count","in_cell_DMSO1_mutation_rate","in_cell_1M71_read_count","in_cell_1M71_mutation_rate","DC_read_count" ,"DC_mutation_rate"))
dc_exvivo <- read.delim("Data/Xist_dataset/XIST_1M7_ex-vivo_wDC.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","ex_vivo_DMSO1_read_count","ex_vivo_DMSO1_mutation_rate","ex_vivo_1M71_read_count","ex_vivo_1M71_mutation_rate","DC_read_count" ,"DC_mutation_rate"))



# READ COUNTS  i.e COVERAGE
incell_counts <- data.frame("in_cell_DMSO1_read_count" = table1_incell["in_cell_DMSO1_read_count"],"in_cell_DMSO2_read_count" = table2_incell["in_cell_DMSO2_read_count"],"X1M7_1_read_count" = table1_incell["in_cell_1M71_read_count"],"X1M7_2_read_count"= table2_incell["in_cell_1M72_read_count"])
exvivo_counts <- data.frame("ex_vivo_DMSO1_read_count" = table1_exvivo["ex_vivo_DMSO1_read_count"],"ex_vivo_DMSO2_read_count" = table2_exvivo["ex_vivo_DMSO2_read_count"],"X1M7_1_read_count" = table1_exvivo["ex_vivo_1M71_read_count"],"X1M7_2_read_count"= table2_exvivo["ex_vivo_1M72_read_count"])

head(incell_counts)
head(exvivo_counts)

# MUTATION RATES
incell_rates <- data.frame("in_cell_DMSO1_mutation_rate" = table1_incell["in_cell_DMSO1_mutation_rate"],"in_cell_DMSO2_mutation_rate" = table2_incell["in_cell_DMSO2_mutation_rate"],"X1M7_1_mutation_rate" = table1_incell["in_cell_1M71_mutation_rate"],"X1M7_2_mutation_rate"= table2_incell["in_cell_1M72_mutation_rate"])
exvivo_rates <- data.frame("ex_vivo_DMSO1_mutation_rate" = table1_exvivo["ex_vivo_DMSO1_mutation_rate"],"ex_vivo_DMSO2_mutation_rate" = table2_exvivo["ex_vivo_DMSO2_mutation_rate"],"X1M7_1_mutation_rate" = table1_exvivo["ex_vivo_1M71_mutation_rate"],"X1M7_2_mutation_rate"= table2_exvivo["ex_vivo_1M72_mutation_rate"])

# DENATURED CONTROLS 
dc_incell_column <- data.frame("DC_mutation_rate"=dc_incell["DC_mutation_rate"] )
dc_exvivo_column <- data.frame("DC_mutation_rate"=dc_exvivo["DC_mutation_rate"])


# NORMALIZATION BY DENATURED CONTROLS 
incell_rates = formattable(incell_rates,digits = 8, format = "f" )
exvivo_rates = formattable(exvivo_rates,digits = 8, format = "f" )

dc_incell_column = formattable(dc_incell_column, digits = 8, format = "f" )
dc_exvivo_column = formattable(dc_exvivo_column,digits = 8, format = "f" )


#Dividing DMSO and 1MT columns for both in cell and ex vivo by respective DC replicate
#Note:only 1 single replicate for DC was available, for each condition. 
scaled_incell_rates <- incell_rates[,] / dc_incell_column
scaled_exvivo_rates <- exvivo_rates[,] / dc_exvivo_column

#Removing NaN and Inf values
is.na(scaled_incell_rates)<-sapply(scaled_incell_rates, is.infinite)
scaled_incell_rates[is.na(scaled_incell_rates)]<-0

is.na(scaled_exvivo_rates)<-sapply(scaled_exvivo_rates, is.infinite)
scaled_exvivo_rates[is.na(scaled_exvivo_rates)]<-0

attributes(scaled_incell_rates) <- NULL
attributes(scaled_exvivo_rates) <- NULL

n <- length(scaled_incell_rates[[1]])
scaled_incell_rates_df  <- structure(scaled_incell_rates,  row.names = c(NA, -n), class = "data.frame")
colnames(scaled_incell_rates_df) <- c("in_cell_DMSO1_mutation_rate", "in_cell_DMSO2_mutation_rate", "in_cell_1M71_mutation_rate", "in_cell_1M72_mutation_rate" )
scaled_exvivo_rates_df  <- structure(scaled_exvivo_rates,  row.names = c(NA, -n), class = "data.frame")
colnames(scaled_exvivo_rates_df) <- c("in_cell_DMSO1_mutation_rate", "in_cell_DMSO2_mutation_rate", "in_cell_1M71_mutation_rate", "in_cell_1M72_mutation_rate" )

#Setting regions to 0 

scaled_incell_rates_df[2500:4500,]=0 
scaled_exvivo_rates_df[2500:4500,]=0 

scaled_incell_rates_df[1:78,]=0 
scaled_exvivo_rates_df[1:78,]=0 

scaled_incell_rates_df[2451:2599,]=0 
scaled_exvivo_rates_df[2451:2599,]=0 

scaled_incell_rates_df[17801:17918,]=0 
scaled_exvivo_rates_df[17801:17918,]=0 

# Back-calculate mutation counts using scaled rates 
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

#### ------- DATA PRE-PROCESSING (CALCULATING LOG RATIOS OF MUTATION RATES AND EMPIRICAL P-VALUES ) ------- ######

# Calculating log mutation rate ratios (denoted as drop offs from here on)
logdropoffs_incell <- calculateLDRs(incell_counts,mutation_counts_in_cell, noreplicates, refsequence)
logdropoffs_exvivo <- calculateLDRs(exvivo_counts,mutation_counts_ex_vivo, noreplicates, refsequence)


## ------------------------------------------------------------------------##

head(logdropoffs_incell$LDR_C)
head(logdropoffs_incell$LDR_CT)
head(logdropoffs_exvivo$LDR_C)
head(logdropoffs_exvivo$LDR_CT)

Nc <- Nt <- noreplicates

strand = "+"

#Calculation of empirical p values for the two probing conditions
empPvals_1 <- computePvals(logdropoffs_incell$LDR_C,logdropoffs_incell$LDR_CT, Nc, Nt, strand, logdropoffs_incell$nuclPosition,
                           logdropoffs_incell$nuclSelection$analysedC, logdropoffs_incell$nuclSelection$analysedCT)

empPvals_2 <- computePvals(logdropoffs_exvivo$LDR_C,logdropoffs_exvivo$LDR_CT, Nc, Nt, strand, logdropoffs_exvivo$nuclPosition,
                           logdropoffs_exvivo$nuclSelection$analysedC, logdropoffs_exvivo$nuclSelection$analysedCT)

#stretches contain the selection of nulceotide positions that have valid log mutation rate ratios
#overlapping the stretches of nucleotides selected in each group of samples to get a representative set
stretches <-overlapsRanges(logdropoffs_incell$stretches,logdropoffs_exvivo$stretches)



####### ------- HMM ANALYSIS (CALCULATION OF POSTERIOR PROBABILITIES FOR THE FOUR HIDDEN STATES) ------- ########
###computes posterior probabilities of all nucleotides in the stretches specified above.
#The most important arguments are the LMRs, the positions used to compute the
#null distribution, as well as the positions of the selected stretches
#and nucleotide pairs where the LMRs were obtained.


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

#b= -0.9

if (a > 0) {
    if (b > 0) {
        if (c > 0) {
            if (d > 0) {
                if (e > 0) {
                    if (f > 0) {
                        if (g > 0) {
                            if (h > 0) {
                                if (i > 0) {
                                    if (j > 0) {
                                        if (k > 0) {
                                            if (l > 0) {
                                                if (m > 0) {
                                                    if (n > 0) {
                                                        if (o > 0) {
                                                            if (p > 0) {
                                                                
                                                                trans <- matrix(c(a, e, i, m,
                                                                                  b, f, j, n,
                                                                                  c, g, k, o,
                                                                                  d, h, l, p), nrow = 4, ncol = 4, byrow = TRUE)
                                                                
                                                                
                                                                posteriors_diff <- hmmFwbw_differential_two_betas(pvaluesstretch)
                                                                colnames(posteriors_diff) <- c("UU","UM","MU","MM")
                                                                #head(posteriors_diff)
                                                                
                                                                setwd("Analysis/diffBUM-HMM")
                                                                posteriors_diff <- replace(posteriors_diff,is.na(posteriors_diff),-999)
                                                                write.table(posteriors_diff,sep="\t",quote=FALSE,file=outputfilename,col.names = c("UU","UM","MU","MM"), row.names = TRUE)

                                                                
                                                                
                                                                
                                                                
                                                                
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        
                                    }
                                    
                                }
                                
                            }
                            
                        }
                        
                        
                    }
                    
                    
                }
                
            }
            
            
        }
        
        
    }

    
}


    