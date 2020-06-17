getwd()
#this cannot work because we are overwriting slices of matrices which is why each element
#in the list ends up 2x4. Need to first generate a list of empty matrices then fill them in
Pv1 <- matrix(ncol = 1,nrow = length(empPvals_1[,1]))
Pv2 <- matrix(ncol = 1,nrow = length(empPvals_2[,1]))
pvaluesstretch<-list(Pv1, Pv2)
for (i in 1:length(stretches)) {
  if (i>1 & i<=length(stretches)) {
    ## Extract start and end of a current stretch
    stretchStart <- start(stretches)[i]
    stretchEnd <- end(stretches)[i]
    Pv1 <-cbind(Pv1, matrix(nrow = length(empPvals_1[,1]), ncol = (stretchStart - end(stretches[i-1])-1)))
    Pv2 <-cbind(Pv2, matrix(nrow = length(empPvals_2[,1]), ncol = (stretchStart - end(stretches[i-1])-1)))
    Pv1 <- cbind(Pv1,empPvals_1[,stretchStart:stretchEnd])
    Pv2 <- cbind(Pv2,empPvals_2[,stretchStart:stretchEnd])
    pvaluesstretch <-list(Pv1, Pv2)
    next()
  } else {
    ## Extract start and end of a current stretch
    stretchStart <- start(stretches)[i]
    stretchEnd <- end(stretches)[i]
    Pv1 <- cbind(Pv1,empPvals_1[,stretchStart:stretchEnd])
    Pv2 <- cbind(Pv2,empPvals_2[,stretchStart:stretchEnd])
    pvaluesstretch <- list(Pv1,Pv2)
    next()
  }
  return(pvaluesstretch)
}
pvaluesstretch

pvaluesstretch[[1]][,100:200]
pvaluesstretch
pVals <- list(empPvals_1, empPvals_2)
pVals[[1]][, 100:200]
length(pVals[[1]][,1])
pVals
