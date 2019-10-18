
hmmFwbw_differential_two_betas <- function(pValues){
      
  ## Settings for HMM
  ## Set the transition matrix: try amending this for performance 
  trans <- matrix(c(0.9025, 0.19, 0.19, 0.04,
                    0.0475, 0.76, 0.01, 0.16,
                    0.0475, 0.01, 0.76, 0.16,
                    0.0025, 0.04, 0.04, 0.64), nrow = 4, ncol = 4, byrow = TRUE)
  
  ## Set the values for Beta shape parameters in the emission mixture
  ## model
  alpha_P1 = 1
  beta_P1 = 10
  
  alpha_P2 = 1
  beta_P2 = 10
  
  ## Set initial probability to 0.25 in each state
  initialProb = c(0.25, 0.25, 0.25, 0.25)
  ## Number of experiments (treatment-control comparisons)
  nexp <- length(pValues[[1]][, 1])
  
  # calculates number of nucleotides
  nBins <- length(pValues[[1]][1, ]) 
  
  # calculates number of states
  nStates <- length(trans[, 1]) 
  
  ## Log-likelihood of observations given state
  obsLike <- matrix(1, ncol = nBins, nrow = nStates) ##generates a matrix where rows are states and columns are nucleotides
  
  fwdMessage <- matrix(0, ncol = nBins, nrow = nStates)
  bwdMessage <- matrix(0, ncol = nBins, nrow = nStates)
  
  ## Non independence assumption
  ## Calculation of likelihoods (sum over replicates of likelihoods of each experiment)
  #iterates over a sequence of values from 1 to the number of experimental comparisons
  for (index in 1:nexp ) { 
    #iterates over a sequence of values from 1 to the number of nucleotides
    for (index2 in 1:nBins) { 
      # Likelihod for UU
      if (is.na(pValues[[1]][index, index2]) || is.na(pValues[[2]][index, index2])) {
        obsLike[, index2] <- obsLike[, index2]
      } else { 
        obsLike[1, index2] <- obsLike[1, index2] * stats::dbeta(pValues[[1]][index, index2], shape1 = 1, shape2 = 1) * 
          stats::dbeta(pValues[[2]][index, index2], shape1 = 1, shape2 = 1)
        # Likelihod for MU
        obsLike[2, index2] <- obsLike[2, index2] * stats::dbeta(pValues[[1]][index, index2], shape1 = 1, shape2 = 1) *
          stats::dbeta(pValues[[2]][index, index2], shape1 = alpha_P2, shape2 = beta_P2)
        # Likelihod for UM
        obsLike[3, index2] <- obsLike[3, index2] * stats::dbeta(pValues[[1]][index, index2], shape1 = alpha_P1, shape2 = beta_P1) *
          stats::dbeta(pValues[[2]][index, index2], shape1 = 1, shape2 = 1)
        # Likelihod for MM
        obsLike[4, index2] <- obsLike[4, index2] * stats::dbeta(pValues[[1]][index, index2], shape1 = alpha_P1, shape2 = beta_P1) *
          stats::dbeta(pValues[[2]][index, index2], shape1 = alpha_P2, shape2 = beta_P2)
      }
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