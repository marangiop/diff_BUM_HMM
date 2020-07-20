
calculateLDRs <- function(mergedcounts,mergedstarts,noreplicates,refsequence) {
    ### calculates and returns LDRs for the dataset
    if (noreplicates == 2) {
        mergedcounts <- mergedcounts[3:6]
        mergedstarts <- mergedstarts[3:6]
    } else if (noreplicates == 3){
        mergedcounts <- mergedcounts[3:8]
        mergedstarts <- mergedstarts[3:8]
    } else {
        mergedcounts <- mergedcounts[3:10]
        mergedstarts <- mergedstarts[3:10]
    }   
    
    
    
    #calculate the drop off rates for each nucleotide position, drop off rates for treatment should be higher than control
    mergeddors <- mergedstarts / mergedcounts
    mergeddors <- replace(mergeddors, is.na(mergeddors), 0)
    
    #introduce the reference DNA seqeunce, remove \r\n\ characters that can be found
    #in size info of the file (although I have no idea why)
    #store the reference sequence as a DNAString class object, under the name dna
    #here its needed to specify the number of replicates we have for each control and treatment samples
    
    setwd("Reference_sequences")
    
    seq <- gsub("[\r\n\"]", "", readChar(refsequence, file.info(refsequence)$size))
    dna <- DNAString(seq)
    
    setwd('..')
    
    
    
    #construct a list of matrices using the container SummarizedExperiment (se), containing DOC, coverage, and drop-off rate values.
    #each column represents the samples: 1 column for each control or treatment replicate
    #while each row represents features of interest we want to view: the nucleotide position along the reference DNA sequence
    #the rep function gives a handle to extract data from the matrices
    
    se <- SummarizedExperiment(
        list(
            coverage = as.matrix(mergedcounts),
            dropoff_count = as.matrix(mergedstarts),
            dropoff_rate = as.matrix(mergeddors)
        ),
        colData = DataFrame(replicate = rep(c(
            "control", "treatment"
        ), each = noreplicates)),
        rowData = DataFrame(nucl = Views(dna, successiveIRanges(rep(
            1, nchar(dna)
        ))))
    )
    
    col_name_vector = c()
    for (i in 1:noreplicates) {
        name=paste("C", i, sep="")
        col_name_vector=c(col_name_vector,name)
    }
    for (i in 1:noreplicates) {
        name=paste("T", i, sep="")
        col_name_vector=c(col_name_vector,name)
    }
    
    
    
    #colnames(se) <- c('C1', 'C2', 'T1', 'T2')
    colnames(se) <- col_name_vector
    
    
    ##this is a test to see if we can get data for the specified experimental replicates (controls/treatments)
    #for one of the matrices (coverage) we have contained in se
    controls <- se[, se$replicate == "control"]
    treatments <- se[, se$replicate == "treatment"]
    head(assay(controls, 'coverage'))
    head(assay(treatments, 'coverage'))
    
    ##This section selects nucleotide positions in each experimental replicate from which
    #LDRs will be calculated, using the selectNuclPos function. It takes in the coverage
    #and DOC info from se, the numbers of control (Nc) and treatment (Nt) replicates, and
    #user-specified coverage threshold. Nucleotides with coverage < t are eliminated.
    Nc <- Nt <- noreplicates
    t <- 1
    nuclSelection <- selectNuclPos(se, Nc, Nt, t)
    
    ##this is just a calculator that returns the combinations we can make
    #of n elements when taken m at a time. Here it just means we want to find the number
    #of comparisons we can make when there Nc replicates. It returns the indices of the
    #replicates for each possible combination.
    t(combn(Nc, 2))
    
    ##selectNuclPos returns lists that hold the positional information for the nucleotide
    #pairs that can be used to compute log ratios, in analysedC and analysed CT
    length(nuclSelection$analysedC[[1]])
    length(nuclSelection$analysedCT[[1]])
    
    ## Finds the medians of original drop-off rates in each experimental replicate
    apply(assay(se, 'dropoff_rate'), 2, median)
    
    ## This is a normalization strategy, to ensure the distribution of drop-off rates are
    #similar between replicates. Here we scale drop-off rates of the nucleotides selected for pairwise
    #comparisons to have a common median value.
    assay(se, "dropoff_rate") <- scaleDOR(se, nuclSelection, Nc, Nt)
    
    ## just a print out to confirm if medians of scaled drop-off rates between
    #replicates have become more similar
    apply(assay(se, 'dropoff_rate'), 2, median)
    
    ## ------------------------------------------------------------------------
    ##The computeStretches function finds uninterrupted stretches of nucleotides to be
    #used for calculation of posterior probabilites. It returns an IRanges object
    #called "stretches", a matrix which contains the index ranges and width
    #of the uninterrupted stretches.
    
    stretches <- computeStretches(se, t)
    
    ## ------------------------------------------------------------------------
    #corrects the log-ratios for any dependency on coverage, to correct for coverage
    #bias. This is also the part where the log drop-off ratios are calculated.
    varStab <- stabiliseVariance(se, nuclSelection, Nc, Nt)
    
    ## ------------------------------------------------------------------------
    #considering that immediate neighbours of a nucleotide can affect accessibility,
    #we want to remove the effect of sequence on nucleotides' susceptibility
    #to chemical modification.This is done by computing different null distributions
    #for user-defined sequence patterns. This step is omitted for the 18S.
    nuclNum <- 3
    patterns <- nuclPerm(nuclNum)
    patterns
    
    ## ------------------------------------------------------------------------
    ## Extract the DNA sequence
    sequence <- subject(rowData(se)$nucl)
    sequence
    nuclPosition <- findPatternPos(patterns, sequence, '+')
    
    ## ------------------------------------------------------------------------
    #Skip to here for 18S. Here we are using all the nucleotide positions for constructing
    #a single null distribution for quantifying drop-off rate variability. We specify
    #nuclPosition list to contain 1 element ([1]), corresponding to the single stretch
    #of reference DNA sequence
    nuclPosition <- list()
    nuclPosition[[1]] <- 1:nchar(sequence)
    
    ## Start of the stretch
    nuclPosition[[1]][1]
    ## End of the stretch, gives the length of the stretch of nucleotides
    nuclPosition[[1]][length(nuclPosition[[1]])]
    
    return(list("LDR_C" = varStab$LDR_C, "LDR_CT" = varStab$LDR_CT, "stretches" = stretches, "nuclSelection" = nuclSelection, "nuclPosition" =  nuclPosition, "stretches" =  stretches))
}