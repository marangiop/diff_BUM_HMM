##Annotations and notes for the computeProbs pipeline.
##The computepVals function takes in the LDR values, number of replicates, nucleotide
#positions used to derive LDRs. 

computePvals <- function(LDR_C, LDR_CT, Nc, Nt, strand, nuclPosition, analysedC,analysedCT,optimise=NULL) {
  
  ##safety checks to make sure the input data is uniform and in the correct format, 
  ##satisfying the min threshold for bumhmm
  if ((Nc < 2) | (Nt < 2)) {
    stop('Number of control and treatment replicates must be at least 2.')
  }
  else if (any(sapply(analysedC, function(x) length(x) == 0))) {
    stop('All lists of positions selected for pair-wise comparisons should
         be non-empty.')
  }
  else if (any(sapply(analysedCT, function(x) length(x) == 0))) {
    stop('All lists of positions selected for pair-wise comparisons should
         be non-empty.')
  }
  if ((strand != '+') & (strand != '-')) {
    stop('Strand should be either plus or minus, specified with a sign.')
  }
  if (length(nuclPosition) == 0) {
    stop('The list of considered nucleotide positions should be non-empty.
         If patterns are not used, a vector with all positions should be
         provided in the first element of the list.')
  }
  
  #length of analysedC and CT is the number of matrices positional for each comparison 
  # = to the number of replicate comparisons = number of columns in the LDR matrices
  length(analysedC) 
  if (dim(LDR_C)[2] != length(analysedC)) {
    stop('The matrix of control-control LDRs should have as many columns
         as there are control-control comparisons.')
  }
  if (dim(LDR_CT)[2] != length(analysedCT)) {
    stop('The matrix of treatment-control LDRs should have as many columns
         as there are treatment-control comparisons.')
  }
  if (!is.null(optimise) & !is.numeric(optimise)) {
    stop('Please provide a tolerance if shape parameters are to be optimised
         with EM algorithm.')
  }
  else {
    
    ## Number of nucleotides in the sequence = number of rows in LDR_CT
    nNucl <- dim(LDR_CT)[1]
    
    ## Construct a quantiles list
    quantiles <- list()
    
    ## Intervals at which to compute quantiles (quantiles are the cut points dividing the probability distribution): 
    ## the sequence here generates a sequence of numbers, within a 
    ## specified range and interval. Here we combine two sequences tgt
    ## in the same vector. arguments are seq(start, end, interval).
    intervals <- c(seq(0.01, 0.89, 0.01), seq(0.9, 0.999, 0.001))
    
    
    ## Number of intervals = 189, therefore creating 190 groups to split
    # the LDRs
    intervalNum <- length(intervals)
    
    
    message('Computing the quantiles of null distributions...')
    
    ## Precompute quantiles of LDR distributions for different patterns
    #for each of the matrices in nuclPosition
    for (i in 1:length(nuclPosition)) {
      
      ## Make an array named subset of distribution. Arrays can
      # store data in more than two dimensions. array(2,3,4) has
      #matrices with 2 rows and 3 columns respectively. For each column 
      #(LDRs from each possible replicate comparison) in LDR_C:                        
      subDistr <- array()
      for (j in 1:dim(LDR_C)[2]) {
        
        ## observedPos: getting the indices of positions where LDRs are to be observed
        # unlist = extract the all of the positional information stored in nuclPosition$analysedC. 
        #intersect = overlap the information across the different matrices (different comparisons)
        observedPos <- intersect(unlist(nuclPosition[i]),  
                                 analysedC[[j]])           
        ## Extract a subset of the null distribution from LDR_C, corresponding to
        ## the positions selected and recorded in nuclPosition, and iteratively update the subDistr array
        # each element in subDistr 
        #why need to remove the first LDR if there is only 1 comparison?
        subDistr <- c(LDR_C[observedPos, j], subDistr)    
        if (j == 1) {
          subDistr <- subDistr[1:(length(subDistr)-1)]
        }
      }
      
      ## Compute the quantiles where each C-C LDR falls in within the null distribution:
      # The quantile function in the stats package is used. 
      # Now we know at what quantile each discrete C-C LDR value lies in the null distribution
      quantiles[[i]] = as.numeric(stats::quantile(subDistr, intervals))
      
    }
    
    ## Annotate with the names of patterns ????
    names(quantiles) <- names(nuclPosition)
    
    
    ## draw up a matrix for p-values
    empPvals <- matrix(, nrow=nNucl, ncol=Nc * Nt)
    
    
    ## Compute empirical p-values
    
    message('Computing empirical p-values...')
    
    ## For each pattern
    for (i in 1:length(nuclPosition))  {
      
      ## Collect positions common between all treatment-control comparisons
      all_analysedCT <- Reduce(union, analysedCT)
      
      ## Collect all nucleotide positions that correspond to the pattern in
      # all_analysedCT and store in patternPos
      patternPos <- unlist(nuclPosition[i])
      patternPos <- intersect(patternPos, all_analysedCT)
      
      ## Extract treatment-control LDRs at these positions from all
      ## comparisons and store in subDist (note this is different from subDistr)
      subDist <- LDR_CT[patternPos, ]
      
      
      #if the subDist is an empty array, assign these dimensions to subDist
      if (is.null(dim(subDist))) {
        dim(subDist) <- c(length(subDist), 1)
      }
      
      ## Compare each LDR to the quantiles of the null distribution
      ## for that pattern
      pValPattern <- matrix(, nrow=length(patternPos), ncol=Nc * Nt)
      
      ## Select the closest quantile to the LDR:
      # apply(x, margin (don't understand this), function)
      # This applies a specified function to the T-C LDRs in subDist (x).
      # Here we define the function to calculate the quantile closest to the LDR =
      # getting the quantile whose absolute difference with the LDR is the smallest.
      pValPattern <- apply(subDist, c(1, 2),
                           function(x)
                             if (is.na(x)) {NA}
                           else {
                             which.min(abs(unlist(quantiles[[i]]) - x))
                           })
      
      ## Compute the p-value as 1 - quantile
      pValPattern <- apply(pValPattern, c(1, 2),
                           function(x)
                             1 - intervals[x])
      
      ## Store the empirical p-values for these positions
      empPvals[patternPos, ] <- pValPattern
      empPvals <- t(empPvals)
      
      ## Removing NAs and replacing them with ones. CHECK IF THIS SHOULD BE ZERO!
      #empPvals[is.na(empPvals)] = 1
    }
    return (empPvals)
  }
}