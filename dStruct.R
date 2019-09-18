#---------------------------------------
#dStruct INSTALLATION STEPS 

#install.packages("devtools")
#library(devtools)
#install_github("AviranLab/dStruct")
#library(dStruct)   (To check that installation was successful)
#---------------------------------------

setwd("/Users/maran/Desktop/diff_BUM_HMM_Project/Github/diff_BUM_HMM///")  
getwd()
library(dStruct)
library(ggplot2)
library(reshape2)
#---------------
#Reading data
cond <- c("delta5", "Erb1")
a_cnts <- read.table("Data/35S_control_delta5_merged_reads.sgr",
                       sep = "\t", header= T)
a_cov <- read.table("Data/35S_control_delta5_merged_dropoffcounts.sgr",
                    sep = "\t", header= T)

b_cnts <- read.table("Data/35S_control_Erb1_merged_reads.sgr",
                     sep = "\t", header= T)
b_cov <- read.table("Data/35S_control_Erb1_merged_dropoffcounts.sgr",
                    sep = "\t", header= T)






#----------------
#Calculating reactivities
a_rates <- a_cnts[, 3:6]/a_cov[, 3:6]
b_rates <- b_cnts[, 3:6]/b_cov[, 3:6]
a_raw_reac <- (a_rates[, 3:4] - a_rates[, 1:2])/(1 - a_rates[, 1:2])
b_raw_reac <- (b_rates[, 3:4] - b_rates[, 1:2])/(1 - b_rates[, 1:2])

reac <- cbind(a_raw_reac, b_raw_reac)
reac[reac<0] <- 0
reac <- as.data.frame(reac)
colnames(reac) <- c("delta5_1", "delta5_2", "Erb1_1", "Erb1_2")
reac <- apply(reac, 2, two.eight.normalize)

#----------------------
#Differential analysis
result <- dStruct(reac, reps_A = 2, reps_B = 2, 
                  min_length = 1) #Change the search length here.
wres <- subset(result, 
              FDR < 0.05) #Change the FDR level here.


#--------------
#Plotting results
df <- melt(data.frame(reac, n = 1:nrow(reac)), id.vars = "n")
for (i in 1:nrow(wres)) {
  ggsave(paste0(wres$Start[i], "_", wres$Stop[i], ".pdf"), 
         print(ggplot(df, aes(x= n, y = value)) + geom_bar(stat = "identity") +
          xlab("Nucleotide") + ylab("Normalized reactivity") +
    facet_grid(variable~.) + 
      coord_cartesian(ylim = c(0, 3), xlim = c(wres$Start[i], wres$Stop[i]))),
    width = 7, height = 7, units = "in")
}

