library(dStruct)
library(ggplot2)
library(reshape2)

working_directory <-"/Users/maran/Desktop/diff_BUM_HMM_Project/Github/diff_BUM_HMM/"
setwd(working_directory)

#---------------
#Reading data
#cond <- c("delta5", "Erb1")

table1_incell <- read.delim("Data/XIST_1M7_in-cell_rep1.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","DMSO_1","DMSO_2","1M7_1","1M7_2"))
table2_incell <- read.delim("Data/XIST_1M7_in-cell_rep2.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","DMSO_1","DMSO_2","1M7_1","1M7_2"))

table3_exvivo <- read.delim("Data/XIST_1M7_ex-vivo_rep1.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","DMSO_1","DMSO_2","1M7_1","1M7_2"))
table4_exvivo <- read.delim("Data/XIST_1M7_ex-vivo_rep2.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","DMSO_1","DMSO_2","1M7_1","1M7_2"))

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






#a_cnts <- read.table("Data/35S_control_delta5_merged_dropoffcounts.sgr",
	#	                            sep = "\t", header= T)
#a_cov <- read.table("Data/35S_control_delta5_merged_reads.sgr",
#		                        sep = "\t", header= T)

#b_cnts <- read.table("Data/35S_control_Erb1_merged_dropoffcounts.sgr",
		                #          sep = "\t", header= T)
#b_cov <- read.table("Data/35S_control_Erb1_merged_reads.sgr",
#		                        sep = "\t", header= T)

a_rates <- incell_rate
b_rates <- exvivo_rate
#----------------
#Calculating reactivities
#a_rates <- a_cnts[, 3:6]/a_cov[, 3:6]
#b_rates <- b_cnts[, 3:6]/b_cov[, 3:6]

a_raw_reac <- (a_rates[, 3:4] - a_rates[, 1:2])/(1 - a_rates[, 1:2])
b_raw_reac <- (b_rates[, 3:4] - b_rates[, 1:2])/(1 - b_rates[, 1:2])

reac <- cbind(a_raw_reac, b_raw_reac)
reac[reac<0] <- 0
reac <- as.data.frame(reac)
colnames(reac) <- c("A1", "A2", "B1", "B2")
reac <- apply(reac, 2, two.eight.normalize)

write.table(reac,sep="\t",quote=FALSE,file='SHAPE_values_Xist.txt',col.names = c("DMSO_1","DMSO_2","1M7_1","1M7_2"), row.names = TRUE)
#This is deltaSHAPE values - export his out 

#----------------------
#Differential analysis

#REPEAT WITH 0.20 FDR
result <- dStruct(reac, reps_A = 2, reps_B = 2, 
		                    min_length = 1) #Change the search length here.
res <- subset(result, 
	                    FDR < 0.45) #Change the FDR level here.


#--------------
#Plotting results
df <- melt(data.frame(reac, n = 1:nrow(reac)), id.vars = "n")
for (i in 1:nrow(res)) {
	  ggsave(paste0(res$Start[i], "_", res$Stop[i], ".pdf"), 
		          print(ggplot(df, aes(x= n, y = value)) + geom_bar(stat = "identity") +
				          xlab("Nucleotide") + ylab("Normalized reactivity") +
					      facet_grid(variable~.) + 
					            coord_cartesian(ylim = c(0, 3), xlim = c(res$Start[i], res$Stop[i]))),
		     width = 7, height = 7, units = "in")
}

write.table(reac,sep="\t",quote=FALSE,file='output_dStruct_Xist_new_data_reac_table.txt', row.names = TRUE)
write.table(res,sep="\t",quote=FALSE,file='output_dStruct_Xist_new_data_res_table.txt', row.names = FALSE)
