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
# 
# delta5_counts <- read.delim("Data/35S_data/35S_control_delta5_merged_dropoffcounts.sgr", stringsAsFactors=FALSE, col.names= c("chromosome","position","in_cell_DMSO1_read_count","in_cell_DMSO1_mutation_rate","in_cell_1M71_read_count","in_cell_1M71_mutation_rate"))
# delta5_coverage <- read.delim("Data/35S_data/35S_control_delta5_merged_reads.sgr", stringsAsFactors=FALSE, col.names= c("chromosome","position","in_cell_DMSO2_read_count","in_cell_DMSO2_mutation_rate","in_cell_1M72_read_count","in_cell_1M72_mutation_rate"))
# 
# deltaerb1_counts<- read.delim("Data/35S_data/35S_control_Erb1_merged_dropoffcounts.sgr", stringsAsFactors=FALSE, col.names= c("chromosome","position","ex_vivo_DMSO1_read_count","ex_vivo_DMSO1_mutation_rate","ex_vivo_1M71_read_count","ex_vivo_1M71_mutation_rate"))
# deltaerb1_coverage<- read.delim("Data/35S_data/35S_control_Erb1_merged_reads.sgr", stringsAsFactors=FALSE, col.names= c("chromosome","position","ex_vivo_DMSO2_read_count","ex_vivo_DMSO2_mutation_rate","ex_vivo_1M72_read_count","ex_vivo_1M72_mutation_rate"))
# 
a_cnts <- read.table("Data/35S_data/35S_control_delta5_merged_dropoffcounts.sgr",
                     sep = "\t", header= T)
a_cov <- read.table("Data/35S_data/35S_control_delta5_merged_reads.sgr",
                    sep = "\t", header= T)

b_cnts <- read.table("Data/35S_data/35S_control_Erb1_merged_dropoffcounts.sgr",
                     sep = "\t", header= T)
b_cov <- read.table("Data/35S_data/35S_control_Erb1_merged_reads.sgr",
                    sep = "\t", header= T)


a_rates <- a_cnts[, 3:6]/a_cov[, 3:6]
b_rates <- b_cnts[, 3:6]/b_cov[, 3:6]

a_raw_reac <- (a_rates[, 3:4] - a_rates[, 1:2])/(1 - a_rates[, 1:2])
b_raw_reac <- (b_rates[, 3:4] - b_rates[, 1:2])/(1 - b_rates[, 1:2])



a_raw_reac_rates_2_8_norm <- apply(a_raw_reac , 2, two.eight.normalize)
b_raw_reac_rates_2_8_norm <- apply(b_raw_reac, 2, two.eight.normalize)


reac <- cbind(a_raw_reac_rates_2_8_norm , b_raw_reac_rates_2_8_norm)

colnames(reac) <- c("A1", "A2", "B1", "B2")
reac[reac<0] <- 0
reac <- as.data.frame(reac)

setwd("Analysis/dStruct/35S")

#------  WRITE TO FILE AND INSPECT SHAPE REACTIVITIES (OPTIONAL) --------
#write.table(res,sep="\t",quote=FALSE,file='SHAPE_reactivities_dStruct_35S_1nt.txt', row.names = FALSE)


result <- dStruct(reac, reps_A = 2, reps_B = 2, min_length = 1) #Change the search length here

res <- subset(result, FDR < 0.20) #Change the FDR level here.

write.table(res,sep="\t",quote=FALSE,file='output_dStruct_35S_res_table_11nt.txt', row.names = FALSE)


#----- PLOTTING RESULTS OF DSTRUCT (OPTIONAL) -----
#
#df <- melt(data.frame(reac, n = 1:nrow(reac)), id.vars = "n")
#for (i in 1:nrow(res)) {
#	  ggsave(paste0(res$Start[i], "_", res$Stop[i], ".pdf"), #
#		          print(ggplot(df, aes(x= n, y = value)) + geom_bar(stat = "identity") +
#				          xlab("Nucleotide") + ylab("Normalized reactivity") +
#					      facet_grid(variable~.) + 
#					            coord_cartesian(ylim = c(0, 3), xlim = c(res$Start[i], res$Stop[i]))),
#		     width = 7, height = 7, units = "in")
#}

