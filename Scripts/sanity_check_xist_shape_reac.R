working_directory <- getwd()

table1_incell <- read.delim("Data/XIST_1M7_in-cell_rep1.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","in_cell_DMSO1_read_count","in_cell_DMSO1_mutation_rate","in_cell_1M71_read_count","in_cell_1M71_mutation_rate"))
table2_incell <- read.delim("Data/XIST_1M7_in-cell_rep2.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","in_cell_DMSO2_read_count","in_cell_DMSO2_mutation_rate","in_cell_1M72_read_count","in_cell_1M72_mutation_rate"))

table1_exvivo <- read.delim("Data/XIST_1M7_ex-vivo_rep1.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","ex_vivo_DMSO1_read_count","ex_vivo_DMSO1_mutation_rate","ex_vivo_1M71_read_count","ex_vivo_1M71_mutation_rate"))
table2_exvivo <- read.delim("Data/XIST_1M7_ex-vivo_rep2.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","ex_vivo_DMSO2_read_count","ex_vivo_DMSO2_mutation_rate","ex_vivo_1M72_read_count","ex_vivo_1M72_mutation_rate"))

head(table1_incell["in_cell_DMSO1_read_count"])

dc_incell <- read.delim("Xist_1M7_in-cell_wDC.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","in_cell_DMSO1_read_count","in_cell_DMSO1_mutation_rate","in_cell_1M71_read_count","in_cell_1M71_mutation_rate","DC_read_count" ,"DC_mutation_rate"))
dc_exvivo <- read.delim("XIST_1M7_ex-vivo_wDC.txt", stringsAsFactors=FALSE, col.names= c("chromosome","position","ex_vivo_DMSO1_read_count","ex_vivo_DMSO1_mutation_rate","ex_vivo_1M71_read_count","ex_vivo_1M71_mutation_rate","DC_read_count" ,"DC_mutation_rate"))

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


incell_DMSO= data.frame("DMSO_1"=table1_incell["in_cell_DMSO1_mutation_rate"], "DMSO_2"=table2_incell["in_cell_DMSO2_mutation_rate"])
incell_1M7 = data.frame("1M7_1"=table1_incell["in_cell_1M71_mutation_rate"], "DMSO_2"=table2_incell["in_cell_1M72_mutation_rate"])


exvivo_DMSO= data.frame("DMSO_1"=table1_exvivo["ex_vivo_DMSO1_mutation_rate"], "DMSO_2"=table2_exvivo["ex_vivo_DMSO2_mutation_rate"])
exvivo_1M7 = data.frame("1M7_1"=table1_exvivo["ex_vivo_1M71_mutation_rate"], "DMSO_2"=table2_exvivo["ex_vivo_1M72_mutation_rate"])


incell_substracted=incell_1M7 - incell_DMSO
exvivo_substracted=exvivo_1M7 - exvivo_DMSO


library(formattable)
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

#install.packages("devtools")

#library(devtools)
#install_github("AviranLab/dStruct")

scaled_incell_rates_df[2500:4500,]=0 
scaled_exvivo_rates_df[2500:4500,]=0 

scaled_incell_rates_df[1:78,]=0 
scaled_exvivo_rates_df[1:78,]=0 

scaled_incell_rates_df[2451:2599,]=0 
scaled_exvivo_rates_df[2451:2599,]=0 

scaled_incell_rates_df[17801:17918,]=0 
scaled_exvivo_rates_df[17801:17918,]=0 

library(dStruct)

scaled_incell_rates_2_8_norm <- apply(scaled_incell_rates_df, 2, two.eight.normalize)
scaled_exvivo_rates_2_8_norm <- apply(scaled_exvivo_rates_df, 2, two.eight.normalize)


write.table(scaled_incell_rates_2_8_norm,sep="\t",quote=FALSE,file="incell_recalc_SHAPE_reac.txt",col.names = c("in_cell_recalc_SHAPE_reac_rep1","in_cell_recalc_SHAPE_reac_rep2"), row.names = TRUE)
write.table(scaled_exvivo_rates_2_8_norm,sep="\t",quote=FALSE,file="exvivo_recalc_SHAPE_reac.txt",col.names = c("ex_vivo_recalc_SHAPE_reac_rep1","ex_vivo_recalc_SHAPE_reac_rep2"), row.names = TRUE)



