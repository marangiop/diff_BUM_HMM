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

table1_exvivo <- read.delim("Data/Xist_dataset/XIST_1M7_ex-vivo_rep1.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","ex_vivo_DMSO1_read_count","ex_vivo_DMSO1_mutation_rate","ex_vivo_1M71_read_count","ex_vivo_1M71_mutation_rate"))
table2_exvivo <- read.delim("Data/Xist_dataset/XIST_1M7_ex-vivo_rep2.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","ex_vivo_DMSO2_read_count","ex_vivo_DMSO2_mutation_rate","ex_vivo_1M72_read_count","ex_vivo_1M72_mutation_rate"))

head(table1_incell["in_cell_DMSO1_read_count"])

dc_incell <- read.delim("Data/Xist_dataset/Xist_1M7_in-cell_wDC.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","in_cell_DMSO1_read_count","in_cell_DMSO1_mutation_rate","in_cell_1M71_read_count","in_cell_1M71_mutation_rate","DC_read_count" ,"DC_mutation_rate"))
dc_exvivo <- read.delim("Data/Xist_dataset/XIST_1M7_ex-vivo_wDC.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","ex_vivo_DMSO1_read_count","ex_vivo_DMSO1_mutation_rate","ex_vivo_1M71_read_count","ex_vivo_1M71_mutation_rate","DC_read_count" ,"DC_mutation_rate"))

#READ COUNTS  i.e COVERAGE
incell_counts <- data.frame("in_cell_DMSO1_read_count" = table1_incell["in_cell_DMSO1_read_count"],"in_cell_DMSO2_read_count" = table2_incell["in_cell_DMSO2_read_count"],"X1M7_1_read_count" = table1_incell["in_cell_1M71_read_count"],"X1M7_2_read_count"= table2_incell["in_cell_1M72_read_count"])
exvivo_counts <- data.frame("ex_vivo_DMSO1_read_count" = table1_exvivo["ex_vivo_DMSO1_read_count"],"ex_vivo_DMSO2_read_count" = table2_exvivo["ex_vivo_DMSO2_read_count"],"X1M7_1_read_count" = table1_exvivo["ex_vivo_1M71_read_count"],"X1M7_2_read_count"= table2_exvivo["ex_vivo_1M72_read_count"])

head(incell_counts)
head(exvivo_counts)

#MUTATION RATES
incell_rates <- data.frame("in_cell_DMSO1_mutation_rate" = table1_incell["in_cell_DMSO1_mutation_rate"],"in_cell_DMSO2_mutation_rate" = table2_incell["in_cell_DMSO2_mutation_rate"],"X1M7_1_mutation_rate" = table1_incell["in_cell_1M71_mutation_rate"],"X1M7_2_mutation_rate"= table2_incell["in_cell_1M72_mutation_rate"])
exvivo_rates <- data.frame("ex_vivo_DMSO1_mutation_rate" = table1_exvivo["ex_vivo_DMSO1_mutation_rate"],"ex_vivo_DMSO2_mutation_rate" = table2_exvivo["ex_vivo_DMSO2_mutation_rate"],"X1M7_1_mutation_rate" = table1_exvivo["ex_vivo_1M71_mutation_rate"],"X1M7_2_mutation_rate"= table2_exvivo["ex_vivo_1M72_mutation_rate"])

# DENATURED CONTROLS 
dc_incell_column <- data.frame("DC_mutation_rate"=dc_incell["DC_mutation_rate"] )
dc_exvivo_column <- data.frame("DC_mutation_rate"=dc_exvivo["DC_mutation_rate"])


incell_DMSO= data.frame("DMSO_1"=table1_incell["in_cell_DMSO1_mutation_rate"], "DMSO_2"=table2_incell["in_cell_DMSO2_mutation_rate"])
incell_1M7 = data.frame("1M7_1"=table1_incell["in_cell_1M71_mutation_rate"], "DMSO_2"=table2_incell["in_cell_1M72_mutation_rate"])


exvivo_DMSO= data.frame("DMSO_1"=table1_exvivo["ex_vivo_DMSO1_mutation_rate"], "DMSO_2"=table2_exvivo["ex_vivo_DMSO2_mutation_rate"])
exvivo_1M7 = data.frame("1M7_1"=table1_exvivo["ex_vivo_1M71_mutation_rate"], "DMSO_2"=table2_exvivo["ex_vivo_1M72_mutation_rate"])


incell_substracted=incell_1M7 - incell_DMSO
exvivo_substracted=exvivo_1M7 - exvivo_DMSO


incell_rates = formattable(incell_substracted,digits = 8, format = "f" )
exvivo_rates = formattable(exvivo_substracted,digits = 8, format = "f" )

dc_incell_column = formattable(dc_incell_column, digits = 8, format = "f" )
dc_exvivo_column = formattable(dc_exvivo_column,digits = 8, format = "f" )

scaled_incell_rates <- incell_rates[,] / dc_incell_column
scaled_exvivo_rates <- exvivo_rates[,] / dc_exvivo_column


is.na(scaled_incell_rates)<-sapply(scaled_incell_rates, is.infinite)
scaled_incell_rates[is.na(scaled_incell_rates)]<-0

is.na(scaled_exvivo_rates)<-sapply(scaled_exvivo_rates, is.infinite)
scaled_exvivo_rates[is.na(scaled_exvivo_rates)]<-0

attributes(scaled_incell_rates) <- NULL
attributes(scaled_exvivo_rates) <- NULL

n <- length(scaled_incell_rates[[1]])
scaled_incell_rates_df  <- structure(scaled_incell_rates,  row.names = c(NA, -n), class = "data.frame")
colnames(scaled_incell_rates_df) <- c("in_cell_mutation_rate_rep1", "in_cell_mutation_rate_rep2")

scaled_exvivo_rates_df  <- structure(scaled_exvivo_rates,  row.names = c(NA, -n), class = "data.frame")
colnames(scaled_exvivo_rates_df) <- c("ex_vivo_mutation_rate_rep1", "ex_vivo_mutation_rate_rep2")


scaled_incell_rates_df[2500:4500,]=0 
scaled_exvivo_rates_df[2500:4500,]=0 

scaled_incell_rates_df[1:78,]=0 
scaled_exvivo_rates_df[1:78,]=0 

scaled_incell_rates_df[2451:2599,]=0 
scaled_exvivo_rates_df[2451:2599,]=0 

scaled_incell_rates_df[17801:17918,]=0 
scaled_exvivo_rates_df[17801:17918,]=0 


scaled_incell_rates_2_8_norm <- apply(scaled_incell_rates_df, 2, two.eight.normalize)
scaled_exvivo_rates_2_8_norm <- apply(scaled_exvivo_rates_df, 2, two.eight.normalize)


reac <- cbind(scaled_incell_rates_2_8_norm, scaled_exvivo_rates_2_8_norm)

colnames(reac) <- c("A1", "A2", "B1", "B2")
reac[reac<0] <- 0
reac <- as.data.frame(reac)


result <- dStruct(reac, reps_A = 2, reps_B = 2, min_length = 11) #Change the search length here


res <- subset(result, FDR < 0.05) #Change the FDR level here.

setwd("Analysis/dStruct/Xist")

write.table(res,sep="\t",quote=FALSE,file='output_dStruct_with_scaling_Xist_res_table_11nt.txt', row.names = FALSE)


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

