calc_stats <-
function(tabble, prevalence = NULL, positive, ...) {
  # checks
  # using original all.equal checks will fail
  if (!identical(nrow(tabble), ncol(tabble)))
    stop("the table must have nrow = ncol")
  
  # this doesn't really check order
  if (!identical(rownames(tabble), colnames(tabble)))
    stop("the table must the same groups in the same order")
  
  tabble_init <- tabble
  
  print(tabble)
  
  # Calculate Sensitivity ---------------------------------------------------
  
  if (nrow(tabble_init) > 2) {
    tmp <- tabble_init
    tabble <- matrix(NA, 2, 2)
    
    colnames(tabble) <- rownames(tabble) <- c("pos", "neg")
    posCol <- which(colnames(tmp) %in% positive)
    negCol <- which(!(colnames(tmp) %in% positive))
    
    tabble[1, 1] <- sum(tmp[posCol, posCol])
    tabble[1, 2] <- sum(tmp[posCol, negCol])
    tabble[2, 1] <- sum(tmp[negCol, posCol])
    tabble[2, 2] <- sum(tmp[negCol, negCol])
    tabble <- as.table(tabble)
    
    
    pos <- "pos"
    neg <- "neg"
    
    rm(tmp)
  } else {
    pos <- positive
    neg <- rownames(tabble_init)[rownames(tabble_init) != positive]
  }
  
  numer <- sum(tabble[pos, pos])
  denom <- sum(tabble[, pos])
  sens  <- ifelse(denom > 0, numer/denom, NA)
  
  detection_rate <- sum(tabble[pos, pos])/sum(tabble)
  detection_prevalence <- sum(tabble[pos, ])/sum(tabble)
  
  
  # Calculate Specificity ---------------------------------------------------
  
  numer <- sum(tabble[neg, neg])
  denom <- sum(tabble[, neg])
  spec  <- ifelse(denom > 0, numer/denom, NA)
  
  
  # Calculate Prevalence ----------------------------------------------------
  
  if (is.null(prevalence))
    prevalence <- sum(tabble_init[, positive]) / sum(tabble_init)
  
  
  # Calculate PPV/NPV -------------------------------------------------------
  
  ppv <-
    (sens * prevalence) /
    ((sens * prevalence) + ((1 - spec) *(1 - prevalence)))
  
  npv <-
    (spec * (1 - prevalence)) /
    (((1 - sens) * prevalence) + ((spec) * (1 - prevalence)))
  
  
  # Calculate F1 ------------------------------------------------------------
  
  f1 <- 2/(1/sens + 1/ppv)
  
  
  # Calculate d-prime/AUC ---------------------------------------------------
  
  # check for inability to calculate
  if (any(rowSums(tabble) == 0)) {
    d_prime <- NA
    auc <- NA
  }
  else {
    d_prime <- qnorm(sens) - qnorm(1-spec)  # primary calculation
    
    # check if sens/spec 1/0 and fudge with warning
    if (is.infinite(d_prime)) {
      warning('Encountered infinite values for d_prime,
    fudge factor introduced to correct.')
      sens_   <- abs(sens - .000001)
      spec_   <- abs(spec - .000001)
      d_prime <- qnorm(sens_) - qnorm(1 - spec_)
      
      xmax <- max(4, d_prime + 3)
      x <- seq(-3, xmax, 0.05)
      
      vpx <- stats::pnorm(x + stats::qnorm(sens_))
      fpx <- stats::pnorm(x - stats::qnorm(spec_))
    }
    else {
      xmax <- max(4, d_prime + 3)
      x <- seq(-3, xmax, 0.05)
      
      vpx <- stats::pnorm(x + stats::qnorm(sens))
      fpx <- stats::pnorm(x - stats::qnorm(spec))
    }
    
    fpx.diff <- diff(fpx)
    lower.sum <- sum(fpx.diff * vpx[-1])
    upper.sum <- sum(fpx.diff * vpx[-length(vpx)])
    auc <- (lower.sum + upper.sum)/2
    auc <- ifelse(auc < .5, 1 - auc, auc)
    # shortcut auc = pnorm(tab$`D Prime`/sqrt(2))
  }
  
  
  # Return result -----------------------------------------------------------
  
  dplyr::tibble(
    `Sensitivity/Recall/TPR` = sens,
    `Specificity/TNR` = spec,
    `PPV/Precision` = ppv,
    `NPV` = npv,
    `F1/Dice` = f1,
    `Prevalence` = prevalence,
    `Detection Rate` = detection_rate,
    `Detection Prevalence` = detection_prevalence,
    `Balanced Accuracy` = (sens + spec)/2,
    `FDR` = 1 - ppv,
    `FOR`  = 1 - npv,
    `FPR/Fallout`  = 1 - spec,
    `FNR`  = 1 - sens,
    `D Prime` = d_prime,
    `AUC` = auc
  )
}
