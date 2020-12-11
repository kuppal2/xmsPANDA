QRILCimpute <-
function(dataSet.mvs, tune.sigma = 1)
{
  #samples in cols; features in rows
  
  nFeatures = dim(dataSet.mvs)[1]
  nSamples = dim(dataSet.mvs)[2]
  dataSet.imputed = dataSet.mvs
  QR.obj = list()
  
  cl<-makeCluster(detectCores()*0.5)
  clusterExport(cl,"quantile")
  clusterExport(cl,"lm")
  clusterEvalQ(cl,library(tmvtnorm))
  
  results<-parLapply(cl,1:nSamples,function(i,dataSet.mvs){
    #for (i in 1:nSamples) {
    curr.sample = dataSet.mvs[, i]
    pNAs = length(which(is.na(curr.sample)))/length(curr.sample)
    upper.q = 0.95
    q.normal = qnorm(seq(pNAs, upper.q, (upper.q - pNAs)/(upper.q * 10000)), mean = 0, sd = 1)
    q.curr.sample = quantile(curr.sample, probs = seq(0,upper.q, 1e-04), na.rm = T)
    temp.QR = lm(q.curr.sample ~ q.normal)
    QR.obj[[i]] = temp.QR
    mean.CDD = temp.QR$coefficients[1]
    sd.CDD = as.numeric(temp.QR$coefficients[2])
    data.to.imp = rtmvnorm(n = nFeatures, mean = mean.CDD,
                           sigma = sd.CDD * tune.sigma, upper = qnorm(pNAs,
                                                                      mean = mean.CDD, sd = sd.CDD), algorithm = c("gibbs"))
    curr.sample.imputed = curr.sample
    curr.sample.imputed[which(is.na(curr.sample))] = data.to.imp[which(is.na(curr.sample))]
    #dataSet.imputed[, i] = curr.sample.imputed
    #}
    return(curr.sample.imputed)
  },dataSet.mvs=dataSet.mvs)
  
  stopCluster(cl)
  
  
  ##list(dataSet.imputed, QR.obj)
  
  return(results)
}
