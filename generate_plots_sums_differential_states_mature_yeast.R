working_directory <-"/Users/maran/Desktop/diff_BUM_HMM_Project/Github/diff_BUM_HMM/"
setwd(working_directory)  


table1_incell <- read.delim("5S_negative_control_test_identical_conditions_diff_BUM_HMM_analysed.txt", stringsAsFactors=FALSE, col.names= c("UU","UM","MU","MM"))

differentiallymod <- table1_incell[,2] + table1_incell[,3]

png("diff_bum_hmm_output_5S_sum_of_diff_states.png")
plot(differentiallymod, xlab = 'Nucleotide position',
     ylab = 'Probability of modification (UM+MU)',
     main = 'diffBUMHMM output - 5S',
     ylim = c(0,1))
dev.off()



table2_incell <- read.delim("5.8S_negative_control_test_identical_conditions_diff_BUM_HMM_analysed.txt", stringsAsFactors=FALSE, col.names= c("UU","UM","MU","MM"))

differentiallymod_2 <- table2_incell[,2] + table2_incell[,3]

png("diff_bum_hmm_output_5_8S_sum_of_diff_states.png")
plot(differentiallymod_2, xlab = 'Nucleotide position',
     ylab = 'Probability of modification (UM+MU)',
     main = 'diffBUMHMM output - 5.8S',
     ylim = c(0,1))
dev.off()



table3_incell <- read.delim("18S_negative_control_test_identical_conditions_diff_BUM_HMM_analysed.txt", stringsAsFactors=FALSE, col.names= c("UU","UM","MU","MM"))

differentiallymod_3 <- table3_incell[,2] + table3_incell[,3]

png("diff_bum_hmm_output_18S_sum_of_diff_states.png")
plot(differentiallymod_3, xlab = 'Nucleotide position',
     ylab = 'Probability of modification (UM+MU)',
     main = 'diffBUMHMM output - 18S',
     ylim = c(0,1))
dev.off()


table4_incell <- read.delim("25S_negative_control_test_identical_conditions_diff_BUM_HMM_analysed.txt", stringsAsFactors=FALSE, col.names= c("UU","UM","MU","MM"))

differentiallymod_4 <- table4_incell[,2] + table4_incell[,3]

png("diff_bum_hmm_output_25S_sum_of_diff_states.png")
plot(differentiallymod_4, xlab = 'Nucleotide position',
     ylab = 'Probability of modification (UM+MU)',
     main = 'diffBUMHMM output - 25S',
     ylim = c(0,1))
dev.off()


