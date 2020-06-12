
hmmFwbw_differential <- function(pVals, trans, initialProb, alpha, beta){
    
    ## Number of experiments (treatment-control comparisons)
    nexp <- length(pVals[[1]][, 1])
    
    # calculates number of nucleotides
    nBins <- length(pVals[[1]][1, ]) 
    
    # calculates number of states
    nStates <- length(trans[, 1]) 
    
    ## Log-likelihood of observations given state
    
    ##generates a matrix where rows are states and columns are nucleotides
    obsLike <- matrix(1, ncol = nBins, nrow = nStates) 
    
    fwdMessage <- matrix(0, ncol = nBins, nrow = nStates)
    bwdMessage <- matrix(0, ncol = nBins, nrow = nStates)
    
    ## Non independence assumption
    ## Calculation of likelihoods (sum over replicates of likelihoods of eachexperiment)
    
    #iterates over a sequence of values from 1 to the number of experimental comparisons
    for (index in 1:nexp) {  
        #iterates over a sequence of values from 1 to the number of nucleotides
        for (index2 in 1:nBins) { 
            # Likelihod for UU 
            obsLike[1, index2] <- obsLike[1, index2] * stats::dbeta(pVals[[1]][index, index2], shape1 = 1, shape2 = 1) *
                stats::dbeta(pVals[[2]][index, index2], shape1 = 1, shape2 = 1)
            # Likelihod for UM
            obsLike[2, index2] <- obsLike[2, index2] * stats::dbeta(pVals[[1]][index, index2], shape1 = 1, shape2 = 1) *
                stats::dbeta(pVals[[2]][index, index2], shape1 = alpha, shape2 = beta)
            # Likelihod for MU
            obsLike[3, index2] <- obsLike[3, index2] * stats::dbeta(pVals[[1]][index, index2], shape1 = alpha, shape2 = beta) *
                stats::dbeta(pVals[[2]][index, index2], shape1 = 1, shape2 = 1)
            # Likelihod for MM
            obsLike[4, index2] <- obsLike[4, index2] * stats::dbeta(pVals[[1]][index, index2], shape1 = alpha, shape2 = beta) *
                stats::dbeta(pVals[[2]][index, index2], shape1 = alpha, shape2 = beta)
        }
    }
    
    ## Calculation of the forward messages
    fwdMessage[, 1] <- initialProb * obsLike[, 1] 
    fwdMessage[, 1] <- fwdMessage[, 1] / sum(fwdMessage[, 1]) 
    
    for (index in 2:nBins) {
        fwdMessage[, index] <- (trans %*% fwdMessage[, index - 1]) *
            obsLike[, index]    
        fwdMessage[, index] <- fwdMessage[, index] / sum(fwdMessage[, index]) 
    }
    
    ## Calculation of the backward message
    bwdMessage[, nBins] <- 1
    
    for (index in (nBins - 1):1) {
        bwdMessage[, index] <- (trans %*% bwdMessage[, index + 1]) * obsLike[, index]
        bwdMessage[, index] <- bwdMessage[, index] / sum(bwdMessage[, index])
    }
    
    ## Calculation of posteriors
    posterior <- fwdMessage * bwdMessage
    posterior <- posterior / (matrix(1, nrow = length(posterior[, 1]))
                              %*% colSums(posterior))
    return(t(posterior))
    }