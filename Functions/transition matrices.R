#paolo's original transition matrix
transpaolo <- matrix(c(0.8, 0.1, 0.1, 0.1,
                  0.05, 0.6, 0.15, 0.05,
                  0.05, 0.15, 0.6, 0.05,
                  0.1, 0.15, 0.15, 0.8), nrow = 4, ncol = 4, byrow = TRUE)

#applying alina's trans matrix to diffBUM_HMM, this is the one used
trans1 <- matrix(c(0.9025, 0.19, 0.19, 0.04,
                  0.0475, 0.76, 0.01, 0.16,
                  0.0475, 0.01, 0.76, 0.16,
                  0.0025, 0.04, 0.04, 0.64), nrow = 4, ncol = 4, byrow = TRUE)
