#### ------- PACKAGES INSTALLATION AND IMPORT OF HELPER FUNCTIONS ------ ######

# This sripts assumes: R version 3.6.3 (2020-02-29); RStudio Version 1.1.442

install.packages("devtools")
library(devtools)
install_github("AviranLab/dStruct")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("formattable")
install.packages("rstudioapi")

library(dStruct)
library(ggplot2)
library(reshape2)
library(formattable)

library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

setwd('..')
getwd()

#### ------- LOADING DATA ------ ######

table1_incell <- read.delim("Data/Xist_dataset/XIST_1M7_in-cell_rep1.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","in_cell_DMSO1_read_count","in_cell_DMSO1_mutation_rate","in_cell_1M71_read_count","in_cell_1M71_mutation_rate"))
table2_incell <- read.delim("Data/Xist_dataset/XIST_1M7_in-cell_rep2.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","in_cell_DMSO2_read_count","in_cell_DMSO2_mutation_rate","in_cell_1M72_read_count","in_cell_1M72_mutation_rate"))

table3_exvivo <- read.delim("Data/Xist_dataset/XIST_1M7_ex-vivo_rep1.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","ex_vivo_DMSO1_read_count","ex_vivo_DMSO1_mutation_rate","ex_vivo_1M71_read_count","ex_vivo_1M71_mutation_rate"))
table4_exvivo <- read.delim("Data/Xist_dataset/XIST_1M7_ex-vivo_rep2.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","ex_vivo_DMSO2_read_count","ex_vivo_DMSO2_mutation_rate","ex_vivo_1M72_read_count","ex_vivo_1M72_mutation_rate"))




table1_incell_DMSO_1_rate <- table1_incell[,c(4)]
table2_incell_DMSO_1_rate <- table2_incell[,c(4)]
table1_incell_X1M7_1_rate <- table1_incell[,c(6)]
table2_incell_X1M7_1_rate <- table2_incell[,c(6)]

incell_rate <- data.frame("DMSO_1" = table1_incell_DMSO_1_rate,"DMSO_2" = table2_incell_DMSO_1_rate,"X1M7_1" = table1_incell_X1M7_1_rate,"X1M7_2"= table2_incell_X1M7_1_rate)

table3_exvivo_DMSO_1_rate <- table3_exvivo[,c(4)]
table4_exvivo_DMSO_1_rate <- table4_exvivo[,c(4)]
table3_exvivo_X1M7_1_rate <- table3_exvivo[,c(6)]
table4_exvivo_X1M7_1_rate <- table4_exvivo[,c(6)]

exvivo_rate <- data.frame("DMSO_1" = table3_exvivo_DMSO_1_rate,"DMSO_2" = table4_exvivo_DMSO_1_rate,"X1M7_1" = table3_exvivo_X1M7_1_rate,"X1M7_2"= table4_exvivo_X1M7_1_rate)


incell_rate[2500:4500,]=0 
exvivo_rate[2500:4500,]=0 

incell_rate[1:78,]=0 
exvivo_rate[1:78,]=0 

incell_rate[2451:2599,]=0 
exvivo_rate[2451:2599,]=0 

incell_rate[17801:17918,]=0 
exvivo_rate[17801:17918,]=0 



a_rates <- incell_rate
b_rates <- exvivo_rate

a_raw_reac <- (a_rates[, 3:4] - a_rates[, 1:2])/(1 - a_rates[, 1:2])
b_raw_reac <- (b_rates[, 3:4] - b_rates[, 1:2])/(1 - b_rates[, 1:2])

reac <- cbind(a_raw_reac, b_raw_reac)
reac[reac<0] <- 0
reac <- as.data.frame(reac)
colnames(reac) <- c("A1", "A2", "B1", "B2")
reac <- apply(reac, 2, two.eight.normalize)


result <- dStruct(reac, reps_A = 2, reps_B = 2, min_length = 11) #Change the search length here


res <- subset(result, FDR < 0.05) #Change the FDR level here.

setwd("Analysis/dStruct/Xist")

write.table(res,sep="\t",quote=FALSE,file='output_dStruct_no_scaling_Xist_res_table_11nt.txt', row.names = FALSE)



