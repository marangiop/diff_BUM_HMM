source("https://bioconductor.org/biocLite.R")
#biocLite("BUMHMM")
#biocLite("Biostrings")

suppressPackageStartupMessages({ library(BUMHMM)
                                 library(Biostrings)
                                 library(SummarizedExperiment) })

getwd()
setwd("/Volumes/Granneman_Lab_Backup/Dropbox/diffBUM_HMM/")

mergedcounts <- read.table("Data/35S_control_delta5_merged_reads.sgr",comment.char="#",col.names=c("chromosome","position","35S_DMSO_1","35S_DMSO_2","35S_1M7_1","35S_1M7_2"))
mergedstarts <- read.table("Data/35S_control_delta5_merged_dropoffcounts.sgr",comment.char="#",col.names=c("chromosome","position","35S_DMSO_1","35S_DMSO_2","35S_1M7_1","35S_1M7_2"))

#head(mergedcounts)

mergedcounts <- mergedcounts[3:6]
mergedstarts <- mergedstarts[3:6]

head(mergedcounts)

mergeddors <- mergedstarts/mergedcounts
mergeddors <- replace(mergeddors,is.na(mergeddors),0)

head(mergeddors)

refsequence <- "rDNA.seq"
seq <- gsub("[\r\n\"]", "",readChar(refsequence, file.info(refsequence)$size))
dna <- DNAString(seq)
noreplicates <- 2

se <- SummarizedExperiment(
  list(
    coverage=as.matrix(mergedcounts),
    dropoff_count=as.matrix(mergedstarts),
    dropoff_rate=as.matrix(mergeddors)
  ), colData=DataFrame(
    replicate=rep(c("control", "treatment"), each=noreplicates)
  ), rowData=DataFrame(
    nucl=Views(dna, successiveIRanges(rep(1, nchar(dna))))
  ))
colnames(se) <- c('C1', 'C2', 'T1', 'T2')

se

controls <- se[, se$replicate == "control"]
treatments <- se[, se$replicate == "treatment"]
head(assay(controls, 'coverage'))
head(assay(treatments, 'coverage'))

Nc <- Nt <- noreplicates
t <- 1
nuclSelection <- selectNuclPos(se, Nc, Nt, t)
List(nuclSelection)
t(combn(Nc, 2))

length(nuclSelection$analysedC[[1]])
length(nuclSelection$analysedCT[[1]])

## Medians of original drop-off rates in each replicate
apply(assay(se, 'dropoff_rate'), 2, median)

## Scale drop-off rates
assay(se, "dropoff_rate") <- scaleDOR(se, nuclSelection, Nc, Nt)

## Medians of scaled drop-off rates in each replicate
apply(assay(se, 'dropoff_rate'), 2, median)

## ------------------------------------------------------------------------
stretches <- computeStretches(se, t)

## ------------------------------------------------------------------------
head(stretches)
assay(se, 'dropoff_count')[1748,]

## ------------------------------------------------------------------------
varStab <- stabiliseVariance(se, nuclSelection, Nc, Nt)
LDR_C <- varStab$LDR_C
LDR_CT <- varStab$LDR_CT

hist(LDR_C, breaks = 30, main = 'Null distribution of LDRs')

## ------------------------------------------------------------------------
nuclNum <- 3
patterns <- nuclPerm(nuclNum)
patterns

## ------------------------------------------------------------------------
## Extract the DNA sequence
sequence <- subject(rowData(se)$nucl)
sequence
nuclPosition <- findPatternPos(patterns, sequence, '+')
patterns[[1]]
head(nuclPosition[[1]])

## ------------------------------------------------------------------------
nuclPosition <- list()
nuclPosition[[1]] <- 1:nchar(sequence)

## Start of the stretch
nuclPosition[[1]][1]
## End of the stretch
nuclPosition[[1]][length(nuclPosition[[1]])]

## ------------------------------------------------------------------------
posteriors <- computeProbs(LDR_C, LDR_CT, Nc, Nt, '+', nuclPosition,
                           nuclSelection$analysedC, nuclSelection$analysedCT,
                           stretches)

## ------------------------------------------------------------------------
head(posteriors)

## ------------------------------------------------------------------------
shifted_posteriors <- matrix(, nrow=dim(posteriors)[1], ncol=1)
shifted_posteriors[1:(length(shifted_posteriors) - 1)] <-
  posteriors[2:dim(posteriors)[1], 2]

## ------------------------------------------------------------------------
plot(shifted_posteriors, xlab = 'Nucleotide position',
     ylab = 'Probability of modification',
     main = 'BUMHMM output for 35S DMS data set')

## ----eval=FALSE----------------------------------------------------------
## ## Call the function with the additonal tolerance parameter
## posteriors <- computeProbs(LDR_C, LDR_CT, Nc, Nt, '+', nuclPosition,
##                            nuclSelection$analysedC, nuclSelection$analysedCT,
##                            stretches, 0.001)

## ------------------------------------------------------------------------

shifted_posteriors <- replace(shifted_posteriors,is.na(shifted_posteriors),-999)
write.table(shifted_posteriors,sep="\t",quote=FALSE,file="35S_delta5_posteriors.txt")

