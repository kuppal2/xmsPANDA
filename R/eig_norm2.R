eig_norm2 <-
function(rv,plot.trend=FALSE) {
  # UNPUT:
  #   rv - return value from the eig_norm1
  #   if user wants to change the number of bias trends that will be eliminated h.c in rv should
  #   be updates to the desired number
  #
  # OUTPUT:
  #   normalized - matrix of normalized abundances with 2 columns of protein and peptdie names
  #   norm_m - matrix of normalized abundances, no extra columns
  #   eigentrends - found in raw data, bias trendsup to h.c
  #   rescrange - rescaling range for the addition of the while noise to avoid overfitting
  #   norm.svd - trends in normalized data, if one wanted to plot at later time.
  #   exPeps - excluded variables - excluded due to exception in fitting a linear model
  
  m = rv$pres # yuliya: use pres matrix, as we cannot deal with m anyways, need to narrow it down to 'complete' variables
  treatment = rv$treatment
  my.svd = rv$my.svd
  pres = rv$pres
  n.treatment = rv$n.treatment
  n.u.treatment = rv$n.u.treatment
  numFact = dim(rv$treatment)[2]
  print(paste('Unique number of treatment combinations:', n.u.treatment) )
  h.c = rv$h.c
  present = rv$present
  toplot1 = rv$toplot1
  # vector of indicators of variables that threw exeptions
  exPeps = vector(mode = "numeric", length = nrow(pres))
  
  print("Normalizing...")
  treatment = data.frame(treatment) # does this need to be done?
  if(n.u.treatment > 1) {
    lm.fm = makeLMFormula(treatment, 'ftemp')
    mtmp = model.matrix(lm.fm$lm.formula, data=treatment, eval(parse(text=lm.fm$lm.params)))  #contrasts=list(bl="contr.sum", it="contr.sum",Pi="contr.sum", tp="contr.sum"))
  } else {  # have 1 treatment group
    mtmp = treatment # as.numeric(t(treatment))
  }
  
  
  # above needed to know how many values will get back for some matrices
  # create some variables:
  betahat = matrix(NA,nrow=dim(mtmp)[2],ncol=nrow(pres))
  newR = array(NA, c(nrow(pres), ncol(pres))) #, n.treatment))
  norm_m = array(NA, c(nrow(pres), ncol(pres))) # , n.treatment))
  numsamp = dim(pres)[2]
  numpep = dim(pres)[1]
  betahat_n = matrix(NA,nrow=dim(mtmp)[2],ncol=nrow(pres))
  rm(mtmp)
  
  
  V0 = my.svd$v[,1:h.c,drop=F]   # residual eigenvariables
  
  
  if(n.u.treatment == 1) { # got 1 treatment group
    for (ii in 1:nrow(pres)) {
      if(ii%%250 == 0) { print(paste('Processing variable ',ii))  }
      pep = pres[ii, ]
      pos = !is.na(pep)
      peptemp = as.matrix(pep[pos]) # take only the observed values
      resm = rep(NA, numsamp)
      resm[pos] = as.numeric(pep[pos])
      bias = array(NA, numsamp)
      bias[pos] = resm[pos] %*% V0[pos,] %*% t(V0[pos,])
      norm_m[ii, ] = as.numeric(pep - bias)
    }
    
  } else { # got 2+ treatment groups
    for (ii in 1:nrow(pres)) {
      if(ii %% 100 == 0) { print(paste('Processing variable ',ii))  }
      pep = pres[ii, ]
      
      pos = !is.na(pep)
      peptemp = as.matrix(pep[pos]) # take only the observed values, may not be needed in R? but this works
      
      
      ftemp = treatment[pos,]
      ftemp = data.frame(ftemp)
      #### use try, not entirely sure if need for modt, need it for solve lm?!
      options(warn = -1)
      lm.fm = makeLMFormula(ftemp, 'ftemp') # using general function that can accomodate for 1+ number of factors
      modt = try(model.matrix(lm.fm$lm.formula, data=ftemp, eval(parse(text=lm.fm$lm.params))), silent=TRUE)
      options(warn = 0)
      
      if(!inherits(modt, "try-error")) { # do nothing if could not make model matrix
        options(warn = -1)
        # if we are able to solve this, we are able to estimate bias
        bhat =  try(solve(t(modt) %*% modt) %*% t(modt) %*% peptemp)
        options(warn = 0)
        if(!inherits(bhat, "try-error")) {
          betahat[,ii] = bhat
          ceffects = modt %*% bhat  # these are the group effects, from estimated coefficients betahat
          
          resm = rep(NA, numsamp) # really a vector only, not m
          resm[pos] = as.numeric(pep[pos] - ceffects)
          bias = array(NA, numsamp)
          bias[pos] = resm[pos] %*% V0[pos,] %*% t(V0[pos,])
          norm_m[ii, ] = as.numeric(pep - bias)
          
          # yuliya:  but newR should be computed on Normalized data
          resm_n = rep(NA, numsamp)
          bhat_n =  solve(t(modt) %*% modt) %*% t(modt) %*% norm_m[ii, pos]
          betahat_n[,ii] = bhat_n
          ceffects_n = modt %*% bhat_n
          resm_n[pos] = norm_m[ii,pos] - ceffects
          newR[ii, ] = resm_n
        } else {
          print(paste('got exception 2 at variable:', ii, 'should not get here...'))
          exPeps[ii] = 2 # should not get 2 here ever...
        }
      } else {
        print(paste('got exception at variable:', ii))
        exPeps[ii] = 1 # keep track of variables that threw exeptions, check why...
      }
    }
  } # end else - got 2+ treatment groups
  
  
  #####################################################################################
  # rescaling has been eliminated form the code after discussion that bias
  # adds variation and we remove it, so no need to rescale after as we removed what was introduced
  y_rescaled = norm_m # for 1 group normalization only, we do not rescale
  # add column names to y-rescaled, now X1, X2,...
  colnames(y_rescaled) = colnames(pres) # these have same number of cols
  rownames(y_rescaled) = rownames(pres)
  y_resc = data.frame(present, y_rescaled)
  rownames(y_resc) = rownames(pres)  # rownames(rv$normalized)
  final = y_resc # row names are assumed to be UNIQUE, variable IDs are unique
  
  # rows with all observations present
  complete_all = y_rescaled[rowSums(is.na(y_rescaled))==0,,drop=F]
  
  #  x11() # make R open new figure window
  par(mfcol=c(3,2))
  par(mar = c(2,2,2,2))
  # center each variable around zero (subtract its mean across samples)
  # note: we are not changing matrix itself, only centerig what we pass to svd
  complete_all_center = t(scale(t(complete_all), center = TRUE, scale = FALSE))
  toplot3 = svd(complete_all_center)
  
  if(plot.trend==TRUE){
    plot.eigentrends(toplot1, "Raw Data")
    plot.eigentrends(toplot3, "Normalized Data")
  }
  # print("Done with normalization!!!")
  colnames(V0) =  paste("Trend", 1:ncol(V0), sep="_")
  
  maxrange = NULL # no rescaling # data.matrix(maxrange)
  return(list(normalized=final, norm_m=y_rescaled, eigentrends=V0, rescrange=maxrange,
              norm.svd=toplot3, exPeps=exPeps))
}
