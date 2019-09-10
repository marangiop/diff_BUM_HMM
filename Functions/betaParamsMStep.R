## This is a function for internal use by betaParamsEM.R

betaParamsMStep <- function(posteriors, pVals, alpha, beta, tolerance,
                                     logP, log1minusP) {

    ## Move into log space
    zeta <- log(alpha)
    ita <- log(beta)

    ## Current estimate
    current_p <- c(alpha, beta)

    ## Number of iterations
    max_iterations <- 20

    ## Number of experimental replicates
    nexp <- length(pVals[, 1])

    points <- matrix(, nrow=max_iterations + 1, ncol=2)
    points[1, ] <- current_p

    ## Compute the current function value
    function_current <- sum(posteriors[, 2] * ((alpha - 1) * logP + (beta - 1) *
                                                 log1minusP -
                                                 nexp * lbeta(alpha, beta)),
                            na.rm=TRUE)

    function_values <- matrix(, nrow=max_iterations + 1, ncol=1)
    function_values[1] <- function_current

    ## Compute constants from posteriors and p-values
    d <- sum(posteriors[, 2] * logP, na.rm=TRUE)
    c <- nexp * sum(posteriors[, 2], na.rm=TRUE)
    f <- sum(posteriors[, 2] * log1minusP, na.rm=TRUE)

    counter <- 1

    new_p <- c(alpha + 1, beta + 1)

    gamma <- 1
    stopping_flag <- FALSE

    ## While estimates continue changing, the maximum number of iterations is
    ## not reached and function value is still decreasing
    while (all(abs(current_p - new_p) > tolerance) &
           (counter <= max_iterations) &
           (!stopping_flag)) {

        if (counter > 1) {
            current_p <- new_p
        }

        gamma <- 1

        ## First derivatives
        alpha_derivative <- c * exp(zeta) * (digamma(exp(zeta) + exp(ita)) -
                            digamma(exp(zeta))) + d * exp(zeta)
        beta_derivative <- c * exp(ita) * (digamma(exp(zeta) + exp(ita)) -
                           digamma(exp(ita))) + f * exp(ita)

        ## Gradient matrix
        gradient <- c(alpha_derivative, beta_derivative)

        ## Second derivatives
        dalpha2_da <- c * exp(zeta) * (digamma(exp(zeta) + exp(ita)) -
                      digamma(exp(zeta)) + exp(zeta) *
                      (trigamma(exp(zeta) + exp(ita)) -
                      trigamma(exp(zeta)))) + d * exp(zeta)
        dalpha2_db <- c * exp(zeta) * trigamma(exp(zeta) + exp(ita)) * exp(ita)
        dbeta2_da <- c * exp(ita) * trigamma(exp(zeta) + exp(ita)) * exp(zeta)
        dbeta2_db <- c * exp(ita) * (digamma(exp(zeta) + exp(ita)) -
                     digamma(exp(ita)) + exp(ita) *
                     (trigamma(exp(zeta) + exp(ita)) - trigamma(exp(ita)))) +
                     f * exp(ita)

        ## Hessian matrix
        hessian <- matrix(c(dalpha2_da, dalpha2_db, dbeta2_da, dbeta2_db),
                          nrow=2, ncol=2, byrow=TRUE)

        ## New guess
        new_p <- current_p - gamma * solve(hessian) %*% gradient
        function_old <- function_current

        if (any(new_p < 0)) {
            function_current <- -Inf
        } else {
            function_current <- sum(posteriors[, 2] * ((new_p[1] - 1) * logP +
                                (new_p[2] - 1) * log1minusP - nexp *
                                lbeta(new_p[1], new_p[2])), na.rm=TRUE)
        }

        ## While function value is decreasing and step size is not too small
        while ((function_current < function_old) & (gamma > 1/256)) {

            ## Half the step size and compute new point
            gamma <- gamma / 2
            new_p <- current_p - gamma * solve(hessian) %*% gradient

            if (any(new_p < 0)) {
                function_current <- -Inf
            } else {
                function_current <- sum(posteriors[, 2] * ((new_p[1] - 1) *
                                    logP + (new_p[2] - 1) * log1minusP - nexp *
                                    lbeta(new_p[1], new_p[2])), na.rm=TRUE)
            }
        }

        ## If new function value is greater than the previous
        if (function_current > function_old) {

            ## Store and update current guess
            points[counter + 1, ] <- new_p
            function_values[counter + 1] <- function_current

            ## Update transformation variables
            zeta <- log(new_p[1])
            ita <- log(new_p[2])

            counter <- counter + 1
        } else {
            ## Keep old guess and stop
            stopping_flag <- TRUE
        }
    }

    return(list('param' = points[counter, ],
                'points' = points,
                'function_values' = function_values))
}
