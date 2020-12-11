confusion_matrix2 <-
function(
  prediction,
  target,
  positive = NULL,
  prevalence = NULL,
  return_table = FALSE,
  dnn = c('Predicted', 'Target'),
  longer = FALSE,
  ...
) {
  
  
  init <- data.frame(prediction, target) %>%
    dplyr::mutate_if(is.logical, as.numeric) %>%
    dplyr::mutate_all(as.factor)
  
  
  
  if (any(levels(init$target) != levels(init$prediction))) {
    warning(
      "Levels are not the same for target and prediction.
    \nRefactoring prediction to match. Some statistics may not be available."
    )
    
    init <- init %>%
      dplyr::mutate(prediction = factor(prediction, levels = levels(target)))
  }
  
  prediction <- init$prediction
  target   <- init$target
  
  # changed focus to be on target levels; prediction can have a single class
  # without failure.
  classLevels <- levels(target)
  numLevels   <- length(classLevels)
  
  
  if(numLevels == 2 & is.null(positive))  positive <- levels(target)[1]
  
  # create confusion matrix
  
  conf_mat <- table(prediction, target, dnn = dnn)
  
  return(conf_mat)
}
