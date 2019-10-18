nuclPerm <- function(n) {

    if (n <= 0) {
        stop('The length of patterns provided is not a positive number.')
    }
    else {
        ## Create all permutations of 4 nucleobases of length n
        perm <- gtools::permutations(4, n, c('A','T','G','C'),
                                     repeats.allowed=TRUE)

        patterns <- ''
        for (i in 1:dim(perm)[1]) {
            patterns[i] <- Reduce(paste0, perm[i, ])
        }
    }

    return(patterns)

}
