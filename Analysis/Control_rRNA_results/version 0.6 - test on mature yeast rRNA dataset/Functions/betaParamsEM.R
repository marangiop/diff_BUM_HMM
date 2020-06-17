## This function is for internal use by computeProbs.R

betaParamsEM <- function(posteriors, pVals, alpha, beta,
                                     tolerance, strand, stretches,
                                     trans, in_prob) {

    ## Initialise with parameter values
    all_alphas <- alpha
    all_betas <- beta

    ## Dummy previous estimate
    prev_alpha <- alpha + 1
    prev_beta <- beta + 1

    max_iterations <- 10
    counter <- 1

    param_estimations <- list()
    function_values <- list()

    ## Comopute constant values
    logP <- colSums(log(pVals))
    log1minusP <- colSums(log(1 - pVals))

    ## While current estimate is changing a lot compared to previous and
    ## the maximum number of iterations is not reached
    while (((abs(alpha - prev_alpha) > tolerance) |
            (abs(beta - prev_beta) > tolerance)) &
           (counter <= max_iterations)) {

        ## Compute new estimates with M-step
        parameters <- betaParamsMStep(posteriors, pVals, alpha, beta, tolerance,
                                      logP, log1minusP)

        ## Current estimates are now the previous ones
        prev_alpha <- alpha
        prev_beta <- beta

        ## New estimates
        alpha <- parameters$param[1]
        beta <- parameters$param[2]

        ## Store estimates and function values
        param_estimations[[counter]] <- parameters$points
        function_values[[counter]] <- parameters$function_values

        ## Update lists
        all_alphas <- c(all_alphas, alpha)
        all_betas <- c(all_betas, beta)

        ## If new estimates differ from the previous ones
        if (!((alpha == prev_alpha) & (beta == prev_beta))) {

            ## Compute new posterior probabilities with E-step
            for (i in 1:length(stretches)) {

                stretchStart <- start(stretches)[i]
                stretchEnd <- end(stretches)[i]

                # Run the HMM inference for the stretches
                if (strand == '+') {
                    current_pVals <- pVals[, stretchStart:stretchEnd]
                    posterior <- hmmFwbw(current_pVals, trans, in_prob, alpha,
                                          beta)
                    posteriors[stretchStart:stretchEnd, ] <- t(posterior)
                }
                if (strand == '-') {
                    current_pVals <- pVals[, stretchEnd:stretchStart]
                    posterior <- hmmFwbw(current_pVals, trans, in_prob, alpha,
                                          beta)
                    posteriors[stretchStart:stretchEnd, ] <- t(posterior)
                }
            }
        }

        counter <- counter + 1
    }

    return(list("posteriors" = posteriors,
                "parameters" = c(alpha, beta),
                "param_estimations" = param_estimations,
                "function_values" = function_values))
}
