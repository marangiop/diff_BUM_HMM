## This is a hidden function
.selectNuclPos <- function(covFile, docFile, Nc, Nt, t) {

    if ((Nc < 2) | (Nt < 2)) {
        stop('The number of experimental replicates must be at least 2.')
    } else if (t < 0) {
        stop('The minumum coverage threshold must be non-negative.')
    }
    else if (any(is.na(covFile)) | any(is.na(docFile))) {
        stop('The coverage and drop-off count matrices should not have NA
              entries.')
    }
    else {

        ## Find nucleotides with significant coverage in each control replicate
        observedC <- list()
        for (col in 1:Nc) {
            observedC[[col]] <- which(covFile[, col] >= t)
        }

        ## Find nucleotides with significant coverage in each treatment
        ## replicate
        observedT <- list()
        i <- 1
        for (col in (Nc+1):(Nc+Nt)) {
            observedT[[i]] <- which(covFile[, col] >= t)
            i <- i+1
        }

        ## Enumerate all pairs of control replicates
        index <- t(combn(Nc, 2))

        ## Find nucleotides with significant coverage in all pairs of control
        ## experiments
        obsC <- list()
        for (i in 1:dim(index)[1]) {
            obsC[[i]] <- intersect(observedC[[index[i,1]]],
                                   observedC[[index[i,2]]])
        }

        ## Enumerate all pairs of treatment-control replicates
        indexT = t(matrix(c(rep((Nc+1):(Nc+Nt), each=Nc), rep(1:Nc, Nt)), 2,
                          byrow=TRUE))

        ## Find nucleotides with significant coverage in all pairs of
        ## treatment-control experiments
        obsCT <- list()
        for (i in 1:dim(indexT)[1]) {
            obsCT[[i]] <- intersect(observedT[[indexT[i,1] - Nc]],
                                    observedC[[indexT[i,2]]])
        }

        ## Include in the analysis only those nucleotides selected for
        ## control-control comparisons that have drop-off count > 0
        analysedC <- list()
        for (i in 1:length(obsC)) {
            analysedC[[i]] <- obsC[[i]][which((docFile[obsC[[i]],
                                      index[i,1]] > 0)
                                      & (docFile[obsC[[i]], index[i,2]] > 0))]
        }

        ## Include in the analysis only those nucleotides selected for
        ## treatment-control comparisons that have drop-off count > 0 in both
        ## replicates
        analysedCT <- list()
        for (i in 1:length(obsCT)) {
            analysedCT[[i]] <- obsCT[[i]][which((docFile[obsCT[[i]],
                                          indexT[i,1]] > 0) &
                                          (docFile[obsCT[[i]],
                                          indexT[i,2]] > 0))]
        }

        return(list("analysedC" = analysedC,
                    "analysedCT" = analysedCT))
    }
}

## The is the function visible to a user
selectNuclPos <- function(se, Nc, Nt, t) {
  res <- .selectNuclPos(assay(se, "coverage"),
                        assay(se, "dropoff_count"), Nc, Nt, t)
  return(res)
}
