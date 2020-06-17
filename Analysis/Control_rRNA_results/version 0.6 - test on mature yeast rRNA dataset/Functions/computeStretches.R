computeStretches <- function(se, t) {

    if (t < 0) {
      stop('The minumum coverage threshold must be non-negative.')
    }

    ## Find positions with coverage >= t in all replicates
    allowedCoverage <- rowSums(assay(se, "coverage") >= t) == ncol(se)

    ## Find positions with drop-off count > 0 in treatment replicates
    tReps <- (se$replicate == "treatment")
    tDOC <- assay(se, "dropoff_count")[, tReps, drop=FALSE] > 0

    ## Find positions with coverage >= t in all replicates and
    ## drop-off count > 0 in at least one treatment replicate
    stretches <- IRanges(rowSums(allowedCoverage & tDOC) > 0)

    ## Pick those stretches that are at least 2 positions wide
    stretches <- stretches[which(stretches@width >= 2)]

    return(stretches)
}
