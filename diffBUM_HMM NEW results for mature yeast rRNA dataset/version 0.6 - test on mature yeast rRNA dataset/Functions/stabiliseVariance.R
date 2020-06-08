.stabiliseVariance <- function(covFile, dorFile, analysedC, analysedCT, Nc, Nt)
{

    if ((Nc < 2) | (Nt < 2)) {
        stop('Number of control and treatment replicates must be at least 2.')
    }
    else if ((length(analysedC) == 0) | (length(analysedCT) == 0)) {
        stop('All lists of positions selected for pair-wise comparisons should
             be non-empty.')
    }
    else if (any(sapply(analysedC, function(x) length(x) == 0))) {
        stop('All lists of positions selected for pair-wise comparisons should
             be non-empty.')
    }
    else if (any(sapply(analysedCT, function(x) length(x) == 0))) {
        stop('All lists of positions selected for pair-wise comparisons should
             be non-empty.')
    }
    else if (any(is.na(covFile)) | any(is.na(dorFile))) {
        stop('The coverage and drop-off count matrices should not have NA
             entries.')
    }
    else {
        ## Number of nucleotides
        nNucl <- dim(covFile)[1]
        ## Number of comparisons between controls
        nCompC <- length(analysedC)

        covVector <- array()
        ldRatVector <- array()

        ## Enumerate all pairs of control replicates
        index <- t(combn(Nc, 2))

        for (i in 1:length(analysedC)) {

            ## Extract coverage and drop-off rates of nucleotides in
            ## control-control comparisons
            coverageC <- covFile[analysedC[[i]], index[i, ]]
            dorC <- dorFile[analysedC[[i]], index[i, ]]

            ## Compute average coverage between the pair of control replicates
            ## at the selected nucleotides
            avgCov <- rowMeans(coverageC)
            dim(avgCov) <- c(length(analysedC[[i]]), 1)

            ## Put them in a vector
            covVector <- c(avgCov, covVector)
            if (i == 1) {
                covVector <- covVector[1:(length(covVector) - 1)]
            }

            ## Compute LDRs for the pair of replicates
            logRat <- log(dorC[, 1]) - log(dorC[, 2])
            dim(logRat) <- c(length(analysedC[[i]]), 1)

            ## Put them in a vector
            ldRatVector <- c(logRat, ldRatVector)
            if (i == 1) {
                ldRatVector <- ldRatVector[1:(length(ldRatVector) - 1)]
            }
        }

        ## Sort by average coverage and order coverages and log-ratios
        orderedInd <- sort(covVector, index.return = TRUE)$ix
        orderedCov <- covVector[orderedInd]
        orderedLd <- ldRatVector[orderedInd]

        ## Split all nucleotides in bins of spanning average coverage in
        ## increments of 100
        bins <- list()

        ## The first nucleotide (sorted by coverage) is placed in the first bin
        c <- 1
        bins[c] <- 1

        ## Store the current average coverage
        current_avgCov <- orderedCov[1]

        ## If the greatest average coverage is not 100 greater than the first,
        ## set the increment to be their difference divided by 10
        if ((orderedCov[length(orderedCov)] - orderedCov[1]) > 100) {
            cov_diff <- 100
        } else {
            cov_diff <- (orderedCov[length(orderedCov)] - orderedCov[1]) / 10
        }

        ## Once a nucleotide is found with coverage greater than the current
        ## coverage + increment, place it in the next bin
        for (i in 2:length(orderedCov)) {
            if (orderedCov[i] >= current_avgCov + cov_diff) {
                c <- c + 1
                bins[c] <- i
                current_avgCov <- orderedCov[i]
            }
        }

        bins <- unlist(bins)

        quantiles <- matrix(, length(bins), 1)
        binCoverage <- matrix(, length(bins), 1)

        ## For each bin, compute the 95th quantile of LDRs with subtracted mean,
        ## and find the average coverage in that bin
        for (i in 1:length(bins)) {
            if (i < length(bins)) {
                binnedLd <- orderedLd[bins[i] : (bins[i+1] - 1)]
                q <- stats::quantile(binnedLd, c(0.95))
                quantiles[i] <- as.numeric(q) - mean(binnedLd)
                binCoverage[i] <- mean(orderedCov[bins[i] : (bins[i+1] - 1)])
            } else {
                binnedLd <- orderedLd[bins[length(bins)] : length(orderedCov)]
                q <- stats::quantile(binnedLd, c(0.95))
                quantiles[i] <- as.numeric(q) - mean(binnedLd)
                binCoverage[i] <- mean(orderedCov[bins[length(bins)] :
                                                    length(orderedCov)])
            }
        }

        ## Find non-linear least-squares estimates for the parameters of a model
        ## for the distribution of LDRs

        if (all(quantiles == 0)) {
            stop("Unable to fit the model for correcting the coverage bias.")
        } else {
            parameters <- stats::nls(quantiles ~ (1/sqrt(binCoverage)) * k + b,
                              start=c(k = 1, b = 0))
            k <- as.numeric(parameters$m$getPars()[1])
            b <- as.numeric(parameters$m$getPars()[2])

            quantileModel <- (1 / sqrt(covVector)) * k + b

            ## Scale LDRs of observed nucleotides according to the model
            ldRatVector <- ldRatVector / quantileModel

            LDR_C <- matrix(, nrow=nNucl, ncol=length(analysedC))

            start <- 0

            ## Arrange these back into a matrix with columns for each comparison
            for (i in 1:length(analysedC)) {
                LDR_C[analysedC[[i]], i] <- ldRatVector[(start + 1):
                                            (start + length(analysedC[[i]]))]
                start <- start + length(analysedC[[i]])
            }

            ## Enumerate all pairs of control replicates
            indexT <- t(matrix(c(rep((Nc+1):(Nc+Nt), each=Nc), rep(1:Nc, Nt)),
                               2, byrow=TRUE))

            ## Number of treatment-control comparisons
            nCompCT <- length(analysedCT)

            LDR_CT <- matrix(, nrow=nNucl, ncol=nCompCT)

            for (i in 1:length(analysedCT)) {

                ## Extract coverage and drop-off rates of nucleotides in
                ## treatment-control comparisons
                coverageCT <- covFile[analysedCT[[i]], indexT[i, ]]
                dorCT <- dorFile[analysedCT[[i]], indexT[i, ]]

                ## Compute average coverage between the pair of replicates at
                ## the selected nucleotides
                avgCov <- rowMeans(coverageCT)
                dim(avgCov) <- c(length(analysedCT[[i]]), 1)

                ## Compute LDRs for the pair of replicates
                logRat <- log(dorCT[, 1]) - log(dorCT[, 2])
                dim(logRat) <- c(length(analysedCT[[i]]), 1)

                # Scale LDRs by the same transformation
                f <- (1 / sqrt(avgCov)) * k + b
                logRat <- logRat / f
                LDR_CT[analysedCT[[i]], i] <- logRat
            }

            return(list("LDR_C" = LDR_C, "LDR_CT" = LDR_CT))

        }
    }
}

stabiliseVariance <- function(se, nuclSelection, Nc, Nt) {
    .stabiliseVariance(assay(se, "coverage"), assay(se, "dropoff_rate"),
                       nuclSelection$analysedC, nuclSelection$analysedCT,
                       Nc, Nt)
}
