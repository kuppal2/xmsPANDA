makeLMFormula <-
function(eff, var_name='') {
  # eff - effects used in contrasts
  # var_name - for singe factor use var-name that is passed in as variable names, otherwise it has no colnmae
  #           only used for a single factor
  if(is.factor(eff))
  {
    ndims = 1
    cols1 = var_name # ftemp in EigenMS
  }
  else
  {
    ndims = dim(eff)[2]
    cols1 = colnames(eff)
  }
  lhs = cols1[1]
  lm.fm = NULL
  # check if can have a list if only have 1 factor...
<<<<<<< HEAD
  
  params = paste('contrasts=list(', cols1[1], '=contr.sum', sep=)
  
=======

  params = paste('contrasts=list(', cols1[1], '=contr.sum', sep=)
 
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
  if (ndims > 1) { # removed ndims[2] here, now ndims holds only 1 dimention...
    for (ii in 2:length(cols1))
    {
      lhs = paste(lhs, "+", cols1[ii])  # bl="contr.sum",
<<<<<<< HEAD
      params = paste(params, ',', cols1[ii], '=contr.sum', sep='')
    }
=======
        params = paste(params, ',', cols1[ii], '=contr.sum', sep='')
      }
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
  }
  params = paste(params,")")
  lm.formula = as.formula(paste('~', lhs))
  lm.fm$lm.formula = lm.formula
  lm.fm$lm.params = params
  return(lm.fm)
}
