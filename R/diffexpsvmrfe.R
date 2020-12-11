diffexpsvmrfe <-
function(x,y,svmkernel="radial"){
  
  #Checking for the variables
  stopifnot(!is.null(x) == TRUE, !is.null(y) == TRUE)
  
  n = ncol(x)
  survivingFeaturesIndexes = seq_len(n)
  featureRankedList = vector(length=n)
  rankedFeatureIndex = n
  
  while(length(survivingFeaturesIndexes)>0){
    #train the support vector machine
    svmModel = svm(x[, survivingFeaturesIndexes], y, type="nu-classification", kernel=svmkernel)
    
    #compute the weight vector
    w = t(svmModel$coefs)%*%svmModel$SV
    
    #compute ranking criteria
    rankingCriteria = w * w
    
    #rank the features
    ranking = sort(rankingCriteria, index.return = TRUE)$ix
    
    #update feature ranked list
    featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]
    rankedFeatureIndex = rankedFeatureIndex - 1
    
    if(length(survivingFeaturesIndexes)==n){
      
      featureWeights=rankingCriteria
      
    }
    
    #eliminate the feature with smallest ranking criterion
    (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
    
    
  }
  
  return (list(featureRankedList=featureRankedList,featureWeights=featureWeights))
}
