stratify_by_class <-
function(yclass,kfold){
  
  y<-yclass
  ytable <- table(y)
  dist_by_class <- ytable/sum(ytable)
  n<-length(y)
  siz <- n - floor(n/kfold)
  classize <- roundvector(dist_by_class * siz, siz)
  if (any(ytable < kfold))
    warning("One or several classes are smaller than the number of folds. \n")
  indlist <- sapply(names(ytable), function(z) which(y ==z), simplify = FALSE)
  
  templist <- vector(mode = "list", length = length(indlist))
  suppressMessages(library(CMA))
  for (i in 1:length(indlist)) {
    outp <- do.call(GenerateLearningsets, args = list(n = ytable[i],
                                                      method = "CV", niter = niter, fold = kfold))@learnmatrix
    templist[[i]] <- t(apply(outp, 1, function(z) ifelse(z ==
                                                           0, 0, indlist[[i]][z])))
  }
  topass <- lapply(templist, function(z) z[1:fold,
                                           , drop = FALSE])
  swaporder <- rowswaps(topass)
  nrep <- 1
  while (nrep < niter) {
    swaporder <- rbind(swaporder, swaporder[1:fold,
                                            , drop = FALSE] + fold * nrep)
    nrep <- nrep + 1
  }
  for (i in 1:length(templist)) templist[[i]] <- templist[[i]][swaporder[,
                                                                         i], ]
  learnmatrix <- templist[[1]]
  for (i in 2:length(indlist)) learnmatrix <- cbind(learnmatrix,
                                                    templist[[i]])
}
