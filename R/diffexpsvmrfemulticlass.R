diffexpsvmrfemulticlass <-
function(x,y,svmkernel="radial"){
  n = ncol(x)
  survivingFeaturesIndexes = seq(n)
  featureRankedList = vector(length=n)
  rankedFeatureIndex = n
  while(length(survivingFeaturesIndexes)>0){
    #train the support vector machine
    svmModel = svm(x[, survivingFeaturesIndexes],
                   y,
                   scale= FALSE,
                   type="nu-classification", kernel=svmkernel)
    #compute the weight vector
    multiclassWeights = svm.getweights(svmModel)
    #compute ranking criteria
    multiclassWeights = multiclassWeights * multiclassWeights
    rankingCriteria = 0
    for(i in 1:ncol(multiclassWeights))rankingCriteria[i] = mean(multiclassWeights[,i])
    #rank the features
    (ranking = sort(rankingCriteria, index.return = TRUE)$ix)
    
    if(length(survivingFeaturesIndexes)==n){
      
      featureWeights=rankingCriteria
      
    }
    
    ## New update feature ranked list
    #featureRankedList[rev((s-r+1):s)] <-survivingFeaturesIndexes[ranking[1:r]]
    featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]
    rankedFeatureIndex = rankedFeatureIndex - 1
    ## New to remove perc.rem
    #  rankedFeatureIndex <- rankedFeatureIndex - r
    
    #eliminate the feature with smallest ranking criterion
    #survivingFeaturesIndexes <-
    #   survivingFeaturesIndexes[-ranking[1:r]]
    
    (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
  }
  #return(featureRankedList)
  return (list(featureRankedList=featureRankedList,featureWeights=featureWeights))
}
