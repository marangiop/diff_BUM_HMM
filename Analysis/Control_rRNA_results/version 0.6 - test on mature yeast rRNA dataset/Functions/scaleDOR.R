## This is a hidden helper function
.poolNucl <- function(c, comparisons, indices) {

  ## Pool nucleotide positions selected for comparisons including replicate c
  rep <- which(comparisons == c, arr.ind=TRUE)[, 1]
  return(Reduce(union, indices[rep]))
}

## This is a hidden function
.scaleDOR <- function(dorFile, nuclSelection, Nc, Nt) {

    if ((Nc < 2) | (Nt < 2)) {
        stop('Number of control and treatment replicates must be at least 2.')
    }
    else if (length(nuclSelection) != 2) {
        stop('Nucleotide selection should have two elements.')
    }
    else if (any(sapply(nuclSelection, function(x) length(x)) == 0)) {
        stop('All lists of positions selected for pair-wise comparisons should be non-empty.')
    }
    else if (any(is.na(dorFile))) {
        stop('Drop-off rate matrix should not have any NA entries.')
    }
    else {
        ## Enumerate control-control comparisons
        index <- t(combn(Nc, 2))

        ## Enumerate treatment-control comparisons
        indexT <- t(matrix(c(rep((Nc+1):(Nc+Nt), each=Nc), rep(1:Nc, Nt)), 2,
                           byrow=TRUE))

        ## Get all positions in control replicates selected for pair-wise
        ## comparisons
        allSelectedNuclC <- list()
        for (i in 1:Nc) {
            allSelectedNuclC[[i]] <- union(
            .poolNucl(i, index, nuclSelection$analysedC),
            .poolNucl(i, indexT, nuclSelection$analysedCT))
        }

        ## Get all positions in treatment replicates selected for pair-wise
        ## comparisons
        allSelectedNuclT <- list()
        for (i in (Nc+1):(Nc+Nt)) {
            allSelectedNuclT[[i - Nc]] <- .poolNucl(i, indexT,
                                                 nuclSelection$analysedCT)
        }

        ## Drop-off rates of selected nucleotides in all control replicates
        ## pooled in one array
        dorC <- array()
        for (i in 1:length(allSelectedNuclC)) {
            dorC <- c(dorFile[allSelectedNuclC[[i]], i], dorC)

            ## Get rid of empty entry at the end
            if (i == 1) {
                dorC <- dorC[1:(length(dorC)-1)]
            }
        }

        ## Compute the median of these drop-off rates
        dorC <- stats::median(dorC)

        ## Compute factors by which the drop-off rates in each replicate must be
        ## scaled to have the same median as above
        factors <- list()

        ## Process control replicates
        for (i in 1:Nc) {
            factors[[i]] <- dorC / stats::median(dorFile[allSelectedNuclC[[i]],
                                                         i])
        }
        ## Process treatment replicates
        for (i in (Nc+1):(Nc + Nt)) {
            factors[[i]] <- dorC / stats::median(dorFile[allSelectedNuclT[[i - Nc]],
                                                  i])
        }

        ## Scale all drop-off rates in each replicate by the corresponding
        ## factor
        for (i in 1:(Nc + Nt)) {
            dorFile[, i] <- dorFile[, i] * factors[[i]]
        }

        return(dorFile)
    }
}

## This is the function visible to a user
scaleDOR <- function(se, nuclSelection, Nc, Nt) {
    .scaleDOR(assay(se, "dropoff_rate"), nuclSelection, Nc, Nt)
}
