library(dStruct)
library(ggplot2)
library(reshape2)

working_directory <-"/Users/maran/Desktop/diff_BUM_HMM_Project/Github/diff_BUM_HMM/"
setwd(working_directory)

#---------------
#Reading data
cond <- c("delta5", "Erb1")
#a_cnts <- read.table("Data/35S_control_delta5_merged_dropoffcounts.sgr",
#                     sep = "\t", header= T)
#a_cov <- read.table("Data/35S_control_delta5_merged_reads.sgr",
#                    sep = "\t", header= T)

#b_cnts <- read.table("Data/35S_control_Erb1_merged_dropoffcounts.sgr",
#                     sep = "\t", header= T)
#b_cov <- read.table("Data/35S_control_Erb1_merged_reads.sgr",
#                    sep = "\t", header= T)


a_cnts <- read.delim("Data/35S_control_delta5_merged_dropoffcounts.sgr", stringsAsFactors=FALSE, col.names= c("chromosome","position","DMSO_1","DMSO_2","1M7_1","1M7_2"))

a_cov <- read.delim("Data/35S_control_delta5_merged_reads.sgr", stringsAsFactors=FALSE, col.names= c("chromosome","position","DMSO_1","DMSO_2","1M7_1","1M7_2"))


b_cnts <- read.delim("Data/35S_control_Erb1_merged_dropoffcounts.sgr", stringsAsFactors=FALSE, col.names= c("chromosome","position","DMSO_1","DMSO_2","1M7_1","1M7_2"))

b_cov <- read.delim("Data/35S_control_Erb1_merged_reads.sgr", stringsAsFactors=FALSE, col.names= c("chromosome","position","DMSO_1","DMSO_2","1M7_1","1M7_2"))




#----------------
#Calculating reactivities
a_rates <- a_cnts[, 3:6]/a_cov[, 3:6]
b_rates <- b_cnts[, 3:6]/b_cov[, 3:6]
a_raw_reac <- (a_rates[, 3:4] - a_rates[, 1:2])/(1 - a_rates[, 1:2])
b_raw_reac <- (b_rates[, 3:4] - b_rates[, 1:2])/(1 - b_rates[, 1:2])

reac <- cbind(a_raw_reac, b_raw_reac)
reac[reac<0] <- 0
reac <- as.data.frame(reac)
colnames(reac) <- c("A1", "A2", "B1", "B2")
reac <- apply(reac, 2, two.eight.normalize)

write.table(reac,sep="\t",quote=FALSE,file='normalised_SHAPE_35S.txt', row.names = TRUE)


ref_seq_directory <- paste(working_directory, "Reference sequences/" ,sep="")
setwd(ref_seq_directory)

seq <- gsub("[\r\n\"]", "", readChar('35S pre-rRNA_refseq.seq', file.info('35S pre-rRNA_refseq.seq')$size))

setwd(working_directory)
seq_split <- strsplit(seq, "")[[1]]

m1 <- matrix(seq_split , ncol=1, byrow=TRUE)
nucleotide_column<- as.data.frame(m1, stringsAsFactors=FALSE)


position_column <- a_cnts[,c(2)]

a1 <- reac[,c(1)]
a2 <- reac[,c(2)]
b1 <- reac[,c(3)]
b2 <- reac[,c(4)]

####
df <- data.frame("Position" = position_column ,"A1" = a1, "nucleotide" = nucleotide_column)
df[,c('standard_error')] <- 0
df <- df[,c(1,2,4,3)]

df[is.na(df)] <- -999.0
df[,2:2][df[, 2:2] == 0] <- -999.0
####

df2 <- data.frame("Position" = position_column ,"A2" = a2, "nucleotide" = nucleotide_column)
df2[,c('standard_error')] <- 0
df2 <- df2[,c(1,2,4,3)]

df2[is.na(df2)] <- -999.0
df2[,2:2][df2[, 2:2] == 0] <- -999.0
####

####
df3 <- data.frame("Position" = position_column ,"B1" = b1, "nucleotide" = nucleotide_column)
df3[,c('standard_error')] <- 0
df3 <- df3[,c(1,2,4,3)]

df3[is.na(df3)] <- -999.0
df3[,2:2][df3[, 2:2] == 0] <- -999.0
####

df4 <- data.frame("Position" = position_column ,"B2" = b2, "nucleotide" = nucleotide_column)
df4[,c('standard_error')] <- 0
df4 <- df4[,c(1,2,4,3)]

df4[is.na(df4)] <- -999.0
df4[,2:2][df4[, 2:2] == 0] <- -999.0
####


write.table(df,sep="\t",quote=FALSE,file='35s_delta5_rep1.map', row.names = FALSE, col.names=FALSE)
write.table(df2,sep="\t",quote=FALSE,file='35s_delta5_rep2.map', row.names = FALSE, col.names=FALSE)

write.table(df3,sep="\t",quote=FALSE,file='35s_deltaerb1_rep1.map', row.names = FALSE, col.names=FALSE)
write.table(df4,sep="\t",quote=FALSE,file='35s_deltaerb1_rep2.map', row.names = FALSE, col.names=FALSE)
