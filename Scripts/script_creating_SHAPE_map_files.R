
working_directory <-"/Users/maran/Desktop/diff_BUM_HMM_Project/Github/diff_BUM_HMM/"
setwd(working_directory)


delta5_shape <- read.table("Data/SHAPE_35S/35S_SHAPE_reactivities.txt",
                    sep = "\t", col.names= c("nucleotide", "rep1", "rep2"))
deltaerb1_shape <- read.table("Data/SHAPE_35S/35S_Erb1_dep_SHAPE_reactivities.txt",
                           sep = "\t", col.names= c("nucleotide","rep1", "rep2"))

delta5_shape[nrow(delta5_shape)+1,] <- c(6868, -999.0,-999.0)
deltaerb1_shape[nrow(deltaerb1_shape)+1,] <- c(6868, -999.0,-999.0)



ref_seq_directory <- paste(working_directory, "Reference sequences/" ,sep="")
setwd(ref_seq_directory)

seq <- gsub("[\r\n\"]", "", readChar('35S pre-rRNA_refseq.seq', file.info('35S pre-rRNA_refseq.seq')$size))

setwd(working_directory)
seq_split <- strsplit(seq, "")[[1]]

m1 <- matrix(seq_split , ncol=1, byrow=TRUE)
nucleotide_column <- as.data.frame(m1, stringsAsFactors=FALSE)


position_column <- delta5_shape[,c(1)]

a1 <- delta5_shape[,c(2)]
a2 <- delta5_shape[,c(3)]
b1 <- deltaerb1_shape[,c(2)]
b2 <- deltaerb1_shape[,c(3)]

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
