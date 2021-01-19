#### ------- PACKAGES INSTALLATION (EXECUTE ONCE)------ ######

# This sripts assumes: R version 3.6.3 (2020-02-29); RStudio Version 1.1.442

install.packages("devtools")
devtools::install_github("AviranLab/dStruct", force = TRUE)
install.packages("ggplot2")
install.packages("reshape2")
install.packages("formattable")

#### ------- SETTING WORKING DIRECTORY AND IMPORT OF HELPER FUNCTIONS (START HERE) ------ ######

library(devtools)
library(dStruct)
library(ggplot2)
library(reshape2)
library(formattable)

wd <- setwd(".")
setwd(wd)

setwd('..')
getwd()

#### ------- LOADING DATA ------ ######
setwd("Reference_sequences")
setwd('../Data/Control_rRNA_dataset')

r1 <- read.table("Fun12_EtOH_160113_1_18S_18S_chemmod_out.txt", col.names=c("gene","position","nucleotide","readsmapped","dropoffs","hybs"))
r2 <- read.table("Fun12_EtOH_160113_2_18S_18S_chemmod_out.txt",  col.names=c("gene","position","nucleotide","readsmapped","dropoffs","hybs"))
r3 <- read.table("Fun12_EtOH_190712_18S_18S_chemmod_out.txt",  col.names=c("gene","position","nucleotide","readsmapped","dropoffs","hybs"))

r4 <- read.table("Fun12_DMS_160113_1_18S_18S_chemmod_out.txt", col.names=c("gene","position","nucleotide","readsmapped","dropoffs","hybs"))
r5 <- read.table("Fun12_DMS_160113_2_18S_18S_chemmod_out.txt",  col.names=c("gene","position","nucleotide","readsmapped","dropoffs","hybs"))
r6 <- read.table("Fun12_DMS_190712_18S_18S_chemmod_out.txt",  col.names=c("gene","position","nucleotide","readsmapped","dropoffs","hybs"))


#mergedcountswt18Snew

mergedcountswt18Snew <- subset(r1, select = -c(nucleotide,dropoffs,hybs, readsmapped ) )
mergedcountswt18Snew$position <- as.integer(as.character(mergedcountswt18Snew$position))
mergedcountswt18Snew$ETOH1 <- as.integer(as.character(r1$readsmapped))
mergedcountswt18Snew$ET0H2 <- as.integer(as.character(r2$readsmapped))
#mergedcountswt18Snew$ETOH3 <- as.integer(as.character(r3$readsmapped))

mergedcountswt18Snew$DMS1 <- as.integer(as.character(r4$readsmapped))
mergedcountswt18Snew$DMS2 <- as.integer(as.character(r5$readsmapped))
#mergedcountswt18Snew$DMS3 <- as.integer(as.character(r6$readsmapped))

mergedcountswt18Snew <- mergedcountswt18Snew[-1,]
mergedcountswt18Snew[nrow(mergedcountswt18Snew)+1,] <- 0

row.names(mergedcountswt18Snew) <- NULL
mergedcountswt18Snew[1800,"gene"] = "18S"
mergedcountswt18Snew[1800,"position"] = 1800


#mergedstartswt18Snew
mergedstartswt18Snew <-  subset(r1, select = -c(nucleotide,dropoffs,hybs, readsmapped ) )
mergedstartswt18Snew$position <- as.integer(as.character(mergedstartswt18Snew$position))


mergedstartswt18Snew$ETOH1 <- as.integer(as.character(r1$dropoffs))
mergedstartswt18Snew$ET0H2 <- as.integer(as.character(r2$dropoffs))
#mergedstartswt18Snew$ETOH3 <- as.integer(as.character(r3$dropoffs))

mergedstartswt18Snew$DMS1 <- as.integer(as.character(r4$dropoffs))
mergedstartswt18Snew$DMS2 <- as.integer(as.character(r5$dropoffs))
#mergedstartswt18Snew$DMS3 <- as.integer(as.character(r6$dropoffs))

mergedstartswt18Snew <- mergedstartswt18Snew[-1,]
mergedstartswt18Snew[nrow(mergedstartswt18Snew)+1,] <- 0

row.names(mergedstartswt18Snew) <- NULL
mergedstartswt18Snew[1800,"gene"] = "18S"
mergedstartswt18Snew[1800,"position"] = 1800

write.table(mergedcountswt18Snew,sep="\t",quote=FALSE,file="mature_rRNA_18SFun12_counts_EtOH1.txt",col.names = c("gene","position", "ETOH1", "ET0H2", "DMS1","DMS2"), row.names = TRUE)
write.table(mergedstartswt18Snew,sep="\t",quote=FALSE,file="mature_rRNA_18SFun12_starts_EtOH1.txt",col.names = c("gene","position", "ETOH1", "ET0H2", "DMS1","DMS2"), row.names = TRUE)

write.table(mergedcountswt18Snew,sep="\t",quote=FALSE,file="mature_rRNA_18SFun12_counts_EtOH2.txt",col.names = c("gene","position", "ETOH1", "ET0H2", "DMS1","DMS2"), row.names = TRUE)
write.table(mergedstartswt18Snew,sep="\t",quote=FALSE,file="mature_rRNA_18SFun12_starts_EtOH2.txt",col.names = c("gene","position", "ETOH1", "ET0H2",  "DMS1","DMS2"), row.names = TRUE)




a_cnts <- read.table("mature_rRNA_18SFun12_starts_EtOH1.txt",
                     sep = "\t", header= T)
a_cov <- read.table("mature_rRNA_18SFun12_counts_EtOH1.txt",
                    sep = "\t", header= T)

b_cnts <- read.table("mature_rRNA_18SFun12_starts_EtOH2.txt",
                     sep = "\t", header= T)
b_cov <- read.table("mature_rRNA_18SFun12_counts_EtOH2.txt",
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

setwd('..')
setwd('..')
setwd("Analysis/dStruct/18S_Fun12")

#------  WRITE TO FILE AND INSPECT SHAPE REACTIVITIES (OPTIONAL) --------

result <- dStruct(reac, reps_A = 2, reps_B = 1, min_length = 5) #Change the search length here

res <- subset(result, FDR < 0.15) #Change the FDR level here.

write.table(res,sep="\t",quote=FALSE,file='output_dStruct_18SFun12_res_table_5nt_0.15FDR.txt', row.names = FALSE)


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

